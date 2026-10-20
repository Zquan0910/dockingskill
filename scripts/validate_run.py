"""Check a real acceptance run against saved structures, command, chemistry and image data."""
import argparse
import json
from pathlib import Path
import numpy as np
from PIL import Image
from rdkit import Chem
from preparation import coordinates, save_json
from workflow import d_path


def validate(job):
    job = d_path(job)
    read = lambda name: json.loads((job / name).read_text(encoding='utf-8'))
    result = read('result.json')
    assert result['status'] == 'SUCCESS'
    assert result['render_status'] == 'SUCCESS'
    grid = read('grid.json')
    xyz = coordinates(job / 'receptor' / 'receptor.pdbqt')
    center, size = np.array(grid['center']), np.array(grid['size'])
    assert np.all(xyz >= center-size/2-1e-9) and np.all(xyz <= center+size/2+1e-9)
    args = read('vina_run.json')['command']
    assert args[args.index('--exhaustiveness')+1] == '64'
    meta = read('structure_selection.json')
    assert meta['alignment']['coverage'] >= .95
    prepared = read('receptor/receptor_metadata.json')
    assert not prepared['missing_prepared_residues']
    assert not prepared['missing_prepared_heavy_atoms']
    assert all(m['present'] for m in prepared['metals'])
    ligand = read('ligand/ligand_metadata.json')
    assert ligand['optimised_conformers'] >= 1
    original = next(m for m in Chem.SDMolSupplier(str(job / 'ligand/ligand_prepared.sdf')) if m)
    exported = next(m for m in Chem.SDMolSupplier(str(job / f'pose_rank{result["best_model"]}.sdf')) if m)
    assert Chem.MolToSmiles(original) == Chem.MolToSmiles(exported)
    contacts = read('figures/interactions.json')['contacts']
    hbonds = [c for c in contacts if c['kind'] == 'hbond']
    for bond in hbonds:
        measured = np.linalg.norm(np.array(bond['PROTCOO'])-np.array(bond['LIGCOO']))
        assert abs(measured-float(bond['DIST_D-A'])) <= .03
    image = np.array(Image.open(job / 'figures/overview.png').convert('RGB'), dtype=float)
    green = (image[:,:,1]-image[:,:,0] > 5) & (image[:,:,1]-image[:,:,2] > 5)
    yy, xx = np.where(green)
    assert len(xx) > 4, 'No visible ligand in overview'
    projection = read('figures/projection.json')
    points = np.array(projection['overview_ligand_pixels'])
    low, high = points.min(axis=0)-20, points.max(axis=0)+20
    assert xx.min() >= low[0] and xx.max() <= high[0] and yy.min() >= low[1] and yy.max() <= high[1], 'Overview marker misses ligand'
    layout = read('figures/layout_qa.json')
    assert layout['annotation_engine'] == 'PyMOL'
    assert projection['annotation_engine'] == 'PyMOL'
    assert projection['label_background_transparency'] == 1
    from pymol import cmd
    cmd.reinitialize()
    cmd.load(str(job / 'figures/scene.pse'))
    native_text = []
    cmd.iterate('annotation_*', 'rows.append(label)', space={'rows': native_text})
    expected_text = [r['text'] for r in projection['labels'] + projection['distances']]
    assert sorted(native_text) == sorted(expected_text), 'Native PSE labels do not match the exported annotations'
    assert cmd.get_setting_float('label_bg_transparency') == 1
    panel = np.array(Image.open(job / 'figures/panel.png').convert('RGB'))
    assert layout['black_frames'] == []
    assert layout['marker_style'] == 'dashed'
    assert layout['connector_color'] == '#e6b800'
    left, right = layout['layout_regions']
    assert left[2]-left[0] == right[2]-right[0]
    for placement, source in zip(layout['image_placements'], ('overview.png', 'detail.png')):
        with Image.open(job / 'figures' / source) as native:
            _, image_size, _ = placement
            assert abs(image_size[0]/image_size[1] - native.width/native.height) < .003
    assert np.count_nonzero(np.all(panel == [230,184,0], axis=2)) > 1000
    assert (job / 'figures/detail.png').read_bytes() == (job / 'figures/detail_annotated.png').read_bytes()
    labels = layout['label_rectangles']
    for i, a in enumerate(labels):
        assert 0 <= a[0] < a[2] <= 1600 and 0 <= a[1] < a[3] <= 900
        for b in labels[i+1:]:
            assert a[2] < b[0] or a[0] > b[2] or a[3] < b[1] or a[1] > b[3]
    checks = {'passed': True, 'atom_count': len(xyz), 'coverage': meta['alignment']['coverage'],
              'identity': meta['alignment']['identity'], 'score': result['best_binding_energy_kcal_mol'],
              'grid_size': size.tolist(), 'exhaustiveness': 64, 'plip_contacts': len(contacts),
              'hydrogen_bonds_distance_checked': len(hbonds), 'overview_marker_contains_visible_ligand': True,
              'native_label_footprints_nonoverlapping': True, 'borderless_yellow_dashed_layout': True,
              'native_pse_labels_verified': len(native_text),
              'had_3d_ligand_optimised': ligand['had_3d']}
    save_json(job / 'verification.json', checks)
    print(json.dumps(checks, indent=2))
    return checks


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('job', type=Path)
    validate(parser.parse_args().job)
