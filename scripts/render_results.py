"""PLIP analysis and reference-style overview/detail panels from actual docked poses."""
from __future__ import annotations
import argparse
import csv
import json
import math
import shutil
import subprocess
import sys
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont
import numpy as np

from workflow import ROOT, d_path
from preparation import save_json
from export_pose import export_pose


def analyse(job, model):
    from plip.structure.preparation import PDBComplex
    from plip.exchange.report import BindingSiteReport, StructureReport
    from plip.basic import config
    folder = job / 'figures'
    folder.mkdir(exist_ok=True)
    config.NOFIX = True
    config.NOHYDRO = False
    complex_ = PDBComplex()
    complex_.output_path = str(folder) + '/'
    complex_.load_pdb(str(job / f'complex_rank{model}.pdb'))
    pose = json.loads((job / 'pose_export.json').read_text())
    targets = [l for l in complex_.ligands if l.hetid == 'LIG' and l.chain == pose['ligand_chain']]
    if len(targets) != 1:
        raise RuntimeError('PLIP did not recognise the exported ligand uniquely.')
    complex_.characterize_complex(targets[0])
    key = f'LIG:{pose["ligand_chain"]}:1'
    if key not in complex_.interaction_sets:
        raise RuntimeError('PLIP interaction set missing.')
    report = BindingSiteReport(complex_.interaction_sets[key])
    StructureReport(complex_).write_xml()
    contacts = []
    for kind in ('hydrophobic', 'hbond', 'waterbridge', 'saltbridge', 'pistacking', 'pication', 'halogen', 'metal'):
        for values in getattr(report, kind + '_info'):
            row = dict(zip(getattr(report, kind + '_features'), values))
            row['kind'] = kind
            if 'PROTCOO' in row and 'LIGCOO' in row:
                row['coordinate_distance'] = float(np.linalg.norm(np.array(row['PROTCOO']) - np.array(row['LIGCOO'])))
            contacts.append(row)
    # JSON roundtrip normalises tuple/numpy values before any further processing.
    contacts = json.loads(json.dumps(contacts, default=lambda x: x.item() if hasattr(x, 'item') else list(x)))
    save_json(folder / 'interactions.json', {'engine': 'PLIP', 'version': config.__version__, 'site': key, 'contacts': contacts,
              'distance_unit': 'angstrom', 'displayed_lines': 'hydrogen-bond donor-acceptor distances only'})
    with (folder / 'interactions.csv').open('w', newline='', encoding='utf-8-sig') as handle:
        fields = sorted({k for r in contacts for k in r}) or ['kind']
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(contacts)
    return contacts, pose


def project(coords, view, width, height, fov):
    points = (np.array(coords) - np.array(view[12:15])) @ np.array(view[:9]).reshape(3, 3)
    scale = height / (2 * abs(view[11]) * math.tan(math.radians(float(fov)) / 2))
    return np.column_stack([width / 2 + (points[:, 0] + view[9]) * scale,
                            height / 2 - (points[:, 1] + view[10]) * scale]).tolist()


def pymol_preflight(folder):
    """Check the same local interpreter used by the native renderer before analysis."""
    probe = subprocess.run([sys.executable, '-B', '-c',
                            'import json,pymol; from pymol import cmd; '
                            'print(json.dumps(dict(version=cmd.get_version()[0],module=pymol.__file__)))'],
                           capture_output=True, text=True, timeout=60)
    if probe.returncode:
        raise RuntimeError('Local PyMOL is unavailable in ' + sys.executable + ': ' + probe.stderr)
    info = json.loads(probe.stdout.strip().splitlines()[-1])
    info['interpreter'] = sys.executable
    save_json(folder / 'pymol_environment.json', info)
    return info


def native_labels(cmd, entries, view):
    """Place transparent PyMOL labels in screen space, then map anchors to the scene."""
    occupied = []
    rotation = np.array(view[:9]).reshape(3, 3)
    scale = 900 / (2 * abs(view[11]) * math.tan(math.radians(float(cmd.get('field_of_view'))) / 2))
    model = cmd.get_model('(ligand or protein) and rep sticks and not hydro')
    xyz = [a.coord for a in model.atom]
    for bond in model.bond:
        a, b = [np.array(model.atom[index].coord) for index in bond.index]
        xyz.extend([(a*(1-t)+b*t).tolist() for t in (.25, .5, .75)])
    obstacles = np.array(project(xyz, view, 1600, 900, cmd.get('field_of_view'))) if xyz else np.empty((0, 2))
    cmd.set('label_color', 'black')
    cmd.set('label_font_id', 5)
    cmd.set('label_size', 23)
    cmd.set('label_shadow_mode', 0)
    cmd.set('label_outline_color', -1)
    cmd.set('label_bg_transparency', 1)
    cmd.set('label_connector', 0)
    cmd.set('float_labels', 1)
    cmd.set('label_position', [0, 0, 0])
    for i, entry in enumerate(entries):
        x, y = project([entry['xyz']], view, 1600, 900, cmd.get('field_of_view'))[0]
        w, h = max(55, len(entry['text']) * 18), 42
        chosen = None
        for radius in (35, 65, 95, 130, 170, 220):
            for angle in (90, 270, 0, 180, 45, 135, 225, 315):
                px = min(1590-w, max(10, x + radius * math.cos(math.radians(angle)) - w/2))
                py = min(890-h, max(10, y + radius * math.sin(math.radians(angle)) - h/2))
                rect = (px, py, px+w, py+h)
                if np.any((obstacles[:,0] > px-12) & (obstacles[:,0] < px+w+12) &
                          (obstacles[:,1] > py-12) & (obstacles[:,1] < py+h+12)):
                    continue
                if all(rect[2] < r[0] or rect[0] > r[2] or rect[3] < r[1] or rect[1] > r[3] for r in occupied):
                    chosen = rect
                    break
            if chosen:
                break
        if chosen is None:
            raise ValueError('Cannot place native PyMOL labels without overlap.')
        occupied.append(chosen)
        # Native labels start at their anchor; reserve a conservative text footprint.
        anchor = [chosen[0]+5, chosen[1]+h/2+7]
        camera_z = ((np.array(entry['xyz'])-np.array(view[12:15])) @ rotation)[2]
        delta = np.array([(anchor[0]-x)/scale, -(anchor[1]-y)/scale, 12-camera_z]) @ rotation.T
        name = f'annotation_{i}'
        cmd.pseudoatom(name, pos=(np.array(entry['xyz']) + delta).tolist())
        cmd.hide('everything', name)
        cmd.label(name, repr(entry['text']))
        entry.update(pixel=[x, y], label_anchor=anchor, label_object=name, rectangle=list(chosen))
    return occupied


def render_scene(job, model):
    import pymol
    from pymol import cmd
    folder = job / 'figures'
    contacts = json.loads((folder / 'interactions.json').read_text())['contacts']
    pose = json.loads((job / 'pose_export.json').read_text())
    cmd.reinitialize()
    cmd.bg_color('white')
    cmd.set('orthoscopic', 1)
    cmd.set('ray_opaque_background', 1)
    cmd.set('ray_shadows', 0)
    cmd.set('specular', 0.45)
    cmd.set('shininess', 45)
    cmd.set('reflect', 0.35)
    cmd.set('direct', 0.7)
    cmd.set('stick_quality', 24)
    cmd.set('depth_cue', 0)
    cmd.set('ray_trace_fog', 0)
    cmd.set('antialias', 2)
    cmd.set('ambient', 0.30)
    cmd.set('ray_trace_mode', 0)
    cmd.load(str(job / f'complex_rank{model}.pdb'), 'complex')
    cmd.select('ligand', f'resn LIG and chain {pose["ligand_chain"]}')
    cmd.select('protein', 'polymer.protein')
    cmd.hide('everything')
    cmd.set_color('protein_cyan', [0.72, 0.86, 0.87])
    cmd.set_color('ligand_green', [0.14, 0.62, 0.24])
    cmd.set_color('residue_red', [0.73, 0.14, 0.12])
    cmd.color('protein_cyan', 'protein')
    cmd.show('surface', 'protein')
    cmd.set('transparency', 0.50, 'protein')
    cmd.show('sticks', 'ligand and not hydro')
    cmd.color('ligand_green', 'ligand and elem C')
    cmd.color('red', 'ligand and elem O')
    cmd.color('blue', 'ligand and elem N')
    cmd.set('stick_radius', 0.32)
    cmd.viewport(1050, 900)
    cmd.orient('protein')
    cmd.zoom('protein or ligand', 5)
    overview_view = list(cmd.get_view())
    overview_xyz = cmd.get_coords('ligand and not hydro').tolist()
    cmd.png(str(folder / 'overview.png'), width=1050, height=900, dpi=300, ray=1)
    overview_view = list(cmd.get_view())
    overview_pixels = project(overview_xyz, overview_view, 1050, 900, cmd.get('field_of_view'))
    cmd.hide('surface')
    cmd.set('stick_radius', 0.24)
    cmd.show('cartoon', 'protein')
    cmd.set('cartoon_transparency', 0.65)
    priority = sorted(contacts, key=lambda r: (r['kind'] != 'hbond', float(r.get('coordinate_distance', 99))))
    display_bonds = [r for r in priority if r['kind'] == 'hbond'][:4]
    residues, shown = [], []
    for r in display_bonds or priority[:4]:
        key = (str(r['RESCHAIN']), str(r['RESNR']), str(r['RESTYPE']))
        if key not in residues and len(residues) < 4:
            residues.append(key)
    for chain, resi, name in residues:
        selection = f'protein and chain {chain} and resi {resi}'
        cmd.show('sticks', selection + ' and not hydro')
        cmd.color('residue_red', selection + ' and elem C')
        cmd.color('red', selection + ' and elem O')
        cmd.color('blue', selection + ' and elem N')
        coords = cmd.get_coords(selection + ' and name CA')
        if coords is not None:
            shown.append({'text': f'{name}-{resi}', 'xyz': coords[0].tolist()})
    distances = []
    for i, r in enumerate(display_bonds):
        if r['kind'] != 'hbond' or (str(r['RESCHAIN']), str(r['RESNR']), str(r['RESTYPE'])) not in residues:
            continue
        a, b = r['PROTCOO'], r['LIGCOO']
        measured = float(np.linalg.norm(np.array(a) - np.array(b)))
        if abs(measured - float(r['DIST_D-A'])) > 0.03:
            raise ValueError('PLIP hydrogen-bond coordinates disagree with the reported distance.')
        cmd.pseudoatom(f'pa{i}', pos=a)
        cmd.pseudoatom(f'pb{i}', pos=b)
        cmd.distance(f'hbond{i}', f'pa{i}', f'pb{i}')
        cmd.hide('everything', f'pa{i} or pb{i}')
        cmd.hide('labels', f'hbond{i}')
        distances.append({'text': f'{measured:.2f} A', 'xyz': ((np.array(a) + np.array(b)) / 2).tolist()})
    cmd.set('dash_color', 'yellow')
    cmd.set('dash_width', 2.5)
    cmd.set('dash_gap', 0.25)
    cmd.select('focus', 'ligand or (byres (protein within 4.0 of ligand))')
    cmd.viewport(1600, 900)
    cmd.orient('focus')
    # Try three principal views; select the one with most ligand spread in the image plane.
    base = cmd.get_view()
    choices = []
    for axis, angle in [('x', 0), ('x', 90), ('y', 90)]:
        cmd.set_view(base)
        cmd.turn(axis, angle)
        view = cmd.get_view()
        pts = (np.array(cmd.get_coords('ligand and not hydro')) - np.array(view[12:15])) @ np.array(view[:9]).reshape(3, 3)
        choices.append((np.ptp(pts[:, 0]) * np.ptp(pts[:, 1]), view))
    cmd.set_view(max(choices, key=lambda t: t[0])[1])
    cmd.zoom('focus', 4.5)
    cmd.clip('slab', 35)
    view = list(cmd.get_view())
    native_labels(cmd, distances + shown, view)
    cmd.png(str(folder / 'detail.png'), width=1600, height=900, dpi=300, ray=1)
    view = list(cmd.get_view())
    for label in shown + distances:
        label['pixel'] = project([label['xyz']], view, 1600, 900, float(cmd.get('field_of_view')))[0]
    save_json(folder / 'projection.json', {'pymol_version': cmd.get_version()[0], 'overview_ligand_pixels': overview_pixels,
              'labels': shown, 'distances': distances, 'view': view,
              'annotation_engine': 'PyMOL', 'label_background_transparency': 1})
    cmd.save(str(folder / 'scene.pse'))


def font(size):
    for path in ['D:/GeneProteinDocking/env/Lib/site-packages/rdkit/Data/Fonts/Amadeus.ttf',
                 'C:/Windows/Fonts/msyh.ttc', 'C:/Windows/Fonts/arial.ttf']:
        if Path(path).is_file() and 'Amadeus' not in path:
            return ImageFont.truetype(path, size)
    return ImageFont.load_default(size=size)


def dashed(draw, start, end, fill, width=2):
    a, b = np.array(start), np.array(end)
    length = float(np.linalg.norm(b - a))
    for offset in range(0, int(length), 16):
        p, q = a + (b-a) * offset / length, a + (b-a) * min(offset+9, length) / length
        draw.line([tuple(p), tuple(q)], fill=fill, width=width)


def compose(job, result):
    folder = job / 'figures'
    overview = Image.open(folder / 'overview.png').convert('RGB')
    detail = Image.open(folder / 'detail.png').convert('RGB')
    projection = json.loads((folder / 'projection.json').read_text())
    occupied = [entry['rectangle'] for entry in projection['distances'] + projection['labels']]
    # Structure annotations are already ray-traced by PyMOL, never painted here.
    shutil.copy2(folder / 'detail.png', folder / 'detail_annotated.png')
    panel = Image.new('RGB', (2800, 1080), 'white')
    # Equal-width layout regions; preserve each native image's aspect ratio.
    left_region, right_region = (40, 110, 1290, 1010), (1510, 110, 2760, 1010)
    placements = []
    for source, region in ((overview, left_region), (detail, right_region)):
        x1, y1, x2, y2 = region
        scale = min((x2-x1)/source.width, (y2-y1)/source.height)
        size = (round(source.width*scale), round(source.height*scale))
        origin = (x1+(x2-x1-size[0])//2, y1+(y2-y1-size[1])//2)
        panel.paste(source.resize(size, Image.Resampling.LANCZOS), origin)
        placements.append((origin, size, scale))
    draw = ImageDraw.Draw(panel)
    title = f"{result.get('gene_name') or result.get('protein_target')}  {result.get('compound_name')}"
    size = 38
    while draw.textbbox((0,0), title, font=font(size))[2] > 2650 and size > 14:
        size -= 1
    draw.text((55, 24), title, font=font(size), fill='#151515')
    points = np.array(projection['overview_ligand_pixels'])
    low, high = points.min(axis=0), points.max(axis=0)
    pixels = np.array(overview, dtype=float)
    green = (pixels[:,:,1]-pixels[:,:,0] > 5) & (pixels[:,:,1]-pixels[:,:,2] > 5)
    yy, xx = np.where(green)
    if len(xx) < 5 or xx.min() < low[0]-20 or xx.max() > high[0]+20 or yy.min() < low[1]-20 or yy.max() > high[1]+20:
        raise ValueError('Overview ligand is not visible inside the projected zoom marker.')
    origin, _, scale = placements[0]
    box = [max(10, low[0]-30)*scale+origin[0], max(10, low[1]-30)*scale+origin[1],
           min(overview.width-10, high[0]+30)*scale+origin[0],
           min(overview.height-10, high[1]+30)*scale+origin[1]]
    if box[2] <= box[0] or box[3] <= box[1]:
        raise ValueError('Ligand projects outside overview; refusing an incorrect zoom marker.')
    yellow = '#e6b800'
    corners = [(box[0],box[1]), (box[2],box[1]), (box[2],box[3]), (box[0],box[3])]
    for start, end in zip(corners, corners[1:] + corners[:1]):
        dashed(draw, start, end, yellow, 4)
    (dx, dy), (dw, dh), _ = placements[1]
    connectors = [((box[2],box[1]), (dx,dy)), ((box[2],box[3]), (dx,dy+dh))]
    for start, end in connectors:
        dashed(draw, start, end, yellow, 4)
    draw.text((1530, 1020), f"Vina score: {result['best_binding_energy_kcal_mol']} kcal/mol", font=font(29), fill='#303030')
    panel.save(folder / 'panel.png', dpi=(300, 300))
    save_json(folder / 'layout_qa.json', {'label_rectangles': occupied, 'label_footprints_do_not_overlap': True,
              'annotation_engine': 'PyMOL', 'label_background_transparency': 1, 'black_frames': [],
              'layout_regions': [left_region, right_region], 'image_placements': placements,
              'connector_color': yellow, 'connectors': connectors, 'marker_style': 'dashed',
              'overview_marker': box, 'overview_ligand_pixels': points.tolist(), 'panel_size': list(panel.size)})
    return folder / 'panel.png'


def assemble_summary(panels, out_dir, rows_per_page=None):
    out_dir = d_path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    if rows_per_page is None:
        from dock_pipeline import read_config
        rows_per_page = int(read_config().get('render', {}).get('rows_per_page', 4))
    if not 1 <= rows_per_page <= 4:
        raise ValueError('Rows per page must be 1..4.')
    pages, names = [], []
    for offset in range(0, len(panels), rows_per_page):
        chunk = panels[offset:offset+rows_per_page]
        page = Image.new('RGB', (2800, 1080 * len(chunk)), 'white')
        for i, path in enumerate(chunk):
            panel = Image.open(path).convert('RGB')
            page.paste(panel, (0, i * 1080))
            ImageDraw.Draw(page).text((8, i*1080+26), chr(65 + (offset+i) % 26), font=font(33), fill='#111111')
        name = out_dir / f'summary_{len(pages)+1:02d}.png'
        page.save(name, dpi=(300,300))
        names.append(str(name))
        pages.append(page)
    if pages:
        pdf = out_dir / 'summary.pdf'
        pages[0].save(pdf, save_all=True, append_images=pages[1:], resolution=300)
        names.append(str(pdf))
    return names


def render_one(job):
    job = d_path(job)
    result = json.loads((job / 'result.json').read_text(encoding='utf-8'))
    if result['status'] != 'SUCCESS':
        raise ValueError('Only successful docking results can be rendered.')
    model = result.get('best_model', 1)
    folder = job / 'figures'
    folder.mkdir(exist_ok=True)
    pymol_preflight(folder)
    export_pose(job, model)
    analyse(job, model)
    folder = job / 'figures'
    # A standalone copy and relative launcher reproduce the same scene after moving a run.
    shutil.copy2(Path(__file__).with_name('render_results.py'), folder / 'render_source.py')
    replay = folder / 'reproduce.py'
    replay.write_text('from pathlib import Path\nimport sys\n'
                     + 'sys.path.insert(0, ' + repr(str(Path(__file__).parent)) + ')\n'
                     + 'from render_source import render_scene\n'
                     + f'render_scene(Path(__file__).resolve().parents[1], {model})\n', encoding='utf-8')
    with (folder / 'pymol.log').open('w', encoding='utf-8') as log:
        proc = subprocess.run([sys.executable, '-B', str(replay)], cwd=folder, stdout=log, stderr=subprocess.STDOUT, timeout=900)
    if proc.returncode:
        raise RuntimeError('PyMOL rendering failed; inspect figures/pymol.log.')
    for name in ('overview.png', 'detail.png'):
        im = np.asarray(Image.open(folder / name).convert('RGB'))
        if float(np.mean(np.min(im, axis=2) < 235)) < 0.003:
            raise RuntimeError(f'Blank render: {name}')
    panel = compose(job, result)
    outputs = [str(panel), str(folder / 'scene.pse'), str(folder / 'interactions.json')]
    outputs += assemble_summary([panel], folder)
    return outputs


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('paths', nargs='+', type=Path, help='Job directories or batch results CSV')
    parser.add_argument('--out-dir', type=Path)
    args = parser.parse_args()
    jobs = []
    for path in args.paths:
        if path.suffix.lower() == '.csv':
            with path.open(encoding='utf-8-sig', newline='') as handle:
                jobs += [Path(r['work_dir']) for r in csv.DictReader(handle) if r.get('status') == 'SUCCESS']
        else: jobs.append(path)
    panels, failed = [], []
    for job in jobs:
        job = d_path(job)
        result = json.loads((job / 'result.json').read_text(encoding='utf-8'))
        try:
            outputs = render_one(job)
            result.update(render_status='SUCCESS', figures=outputs)
            panels.append(job / 'figures' / 'panel.png')
        except Exception as exc:
            import traceback
            (job / 'render_error.log').write_text(traceback.format_exc(), encoding='utf-8')
            result.update(render_status='FAILED', render_error=str(exc))
            failed.append(str(job))
        save_json(job / 'result.json', result)
    out = args.out_dir or ROOT / 'figures' / __import__('uuid').uuid4().hex[:12]
    outputs = assemble_summary(panels, out) if panels else []
    save_json(d_path(out) / 'render_manifest.json', {'outputs': outputs, 'failed_jobs': failed})
    print(json.dumps({'outputs': outputs, 'failed_jobs': failed}, indent=2))
    return 1 if failed or not panels else 0


if __name__ == '__main__':
    raise SystemExit(main())
