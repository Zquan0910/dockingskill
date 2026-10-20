"""Structure selection, traceable preparation, and whole-receptor grid validation."""
from __future__ import annotations

import io
import json
import math
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
from Bio import Align, PDB
from Bio.SeqUtils import seq1
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.MolStandardize import rdMolStandardize

METALS = {'ZN', 'FE', 'MG', 'MN', 'CA', 'CU', 'CO', 'NI'}


def save_json(path, data):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(json.dumps(data, ensure_ascii=False, indent=2), encoding='utf-8')


def load_model(path):
    parser = PDB.MMCIFParser(QUIET=True) if Path(path).suffix.lower() in {'.cif', '.mmcif'} else PDB.PDBParser(QUIET=True)
    return next(parser.get_structure('receptor', str(path)).get_models())


def protein_residues(chain):
    return [r for r in chain if r.id[0] == ' ' and PDB.is_aa(r, standard=False)]


def chain_alignment(chain, target):
    residues = protein_residues(chain)
    sequence = ''.join(seq1(r.resname, custom_map={'MSE': 'M'}) for r in residues)
    if not sequence:
        return {'chain': chain.id, 'coverage': 0, 'identity': 0, 'mapping': []}
    aligner = Align.PairwiseAligner(mode='global', match_score=2, mismatch_score=-1, open_gap_score=-5, extend_gap_score=-0.5)
    alignment = aligner.align(target, sequence)[0]
    mapping, matched = [], 0
    for (a, b), (c, d) in zip(*alignment.aligned):
        for i, j in zip(range(a, b), range(c, d)):
            matched += target[i] == sequence[j]
            mapping.append({'uniprot': int(i + 1), 'pdb_number': residues[j].id[1], 'insertion': residues[j].id[2].strip()})
    return {'chain': chain.id, 'coverage': len(mapping) / len(target),
            'identity': matched / max(1, len(mapping)), 'mapping': mapping,
            'observed_residues': len(residues), 'sequence_length': len(target)}


def select_structure(gene, job, cache, minimum=0.95, local=None, chain=None):
    import dock_pipeline as sources
    notes = []
    if local:
        path = Path(local)
        if not path.is_file():
            path = sources.download_pdb(str(local), cache) if len(str(local)) == 4 else path
        if not path.is_file():
            raise ValueError(f'Structure does not exist: {local}')
        source, pdb_id = 'local_file', path.stem
    elif gene:
        path = None
        for pdb_id in gene.pdb_refs or sources.search_rcsb(gene.accession, cache):
            try:
                candidate = sources.download_pdb(pdb_id, cache)
                matches = [chain_alignment(c, gene.sequence) for c in load_model(candidate)]
                best = max(matches, key=lambda x: (x['coverage'], x['identity']))
                notes.append({'pdb_id': pdb_id, **{k: v for k, v in best.items() if k != 'mapping'}})
                if best['coverage'] >= minimum and best['identity'] >= 0.95:
                    path, source = candidate, 'rcsb'
                    break
            except Exception as exc:
                notes.append({'pdb_id': pdb_id, 'error': str(exc)})
        if path is None:
            path = sources.download_alphafold(gene.accession, cache)
            source, pdb_id = 'alphafold', None
    else:
        raise ValueError('Provide a gene or local protein structure.')
    model = load_model(path)
    matches = [chain_alignment(c, gene.sequence) for c in model] if gene else []
    if chain:
        if chain not in model or not protein_residues(model[chain]):
            raise ValueError(f'Protein chain {chain} is absent.')
        chosen = chain
    elif gene:
        chosen = max(matches, key=lambda x: (x['coverage'], x['identity']))['chain']
    else:
        chains = [c.id for c in model if protein_residues(c)]
        if len(chains) != 1:
            raise ValueError(f'Local structure has {len(chains)} protein chains; specify --chain.')
        chosen = chains[0]
    match = next((m for m in matches if m['chain'] == chosen), None)
    if match and match['identity'] < 0.95:
        raise ValueError('Selected chain does not match the requested protein (identity <95%).')
    if not local and source == 'alphafold' and match['coverage'] < minimum:
        raise ValueError('AlphaFold model does not cover the complete requested target; explicit local input required.')
    raw = job / ('receptor_original' + path.suffix)
    shutil.copy2(path, raw)
    confidence = [float(a.bfactor) for r in protein_residues(model[chosen]) for a in r if a.name == 'CA']
    meta = {'source': source, 'pdb_id': pdb_id, 'accession': gene.accession if gene else None,
            'raw_structure': str(raw), 'target_chain': chosen, 'model': model.id,
            'alignment': match, 'candidates': notes, 'minimum_coverage': minimum,
            'missing_target_positions': sorted(set(range(1, gene.length + 1)) - {m['uniprot'] for m in match['mapping']}) if match else None,
            'warnings': ['No missing-fragment modelling; no pH microstate enumeration.']}
    if source == 'alphafold':
        meta['plddt'] = {'mean': float(np.mean(confidence)), 'fraction_below_70': float(np.mean(np.array(confidence) < 70))}
    save_json(job / 'structure_selection.json', meta)
    return model, chosen, meta


def prepare_receptor(model, chain, job, meta):
    folder = job / 'receptor'
    folder.mkdir()
    protein = model[chain]
    unsupported = [r.resname for r in protein if r.id[0] != ' ' and PDB.is_aa(r, standard=False)]
    if unsupported:
        raise ValueError(f'Nonstandard polymer residues require explicit curation: {unsupported}')
    # Select one physical conformer before writing; all downstream coordinates share it.
    expected = []
    for residue in protein_residues(protein):
        if residue.is_disordered() == 2:
            raise ValueError(f'Disordered residue identity requires explicit curation: {residue.id}')
        expected.append((chain, residue.id[1], residue.id[2]))
        for atom in residue:
            if atom.is_disordered():
                children = atom.child_dict
                selected = 'A' if 'A' in children else max(children, key=lambda k: children[k].occupancy or 0)
                atom.disordered_select(selected)
    xyz = np.array([a.coord for r in protein_residues(protein) for a in r])
    retained_metals = []
    for c in model:
        for r in c:
            if r.resname.strip().upper() in METALS:
                for a in r:
                    if np.min(np.linalg.norm(xyz - a.coord, axis=1)) <= 4.0:
                        retained_metals.append((c.id, r.id, a.name))
    class SelectTarget(PDB.Select):
        def accept_residue(self, residue):
            return (residue.parent.id == chain and residue.id[0] == ' ' and PDB.is_aa(residue, standard=False)) or any(residue.parent.id == c and residue.id == rid for c, rid, _ in retained_metals)
        def accept_atom(self, atom):
            if atom.is_disordered():
                return atom.altloc == atom.parent[atom.name].selected_child.altloc
            return True
    clean = folder / 'receptor_clean.pdb'
    writer = PDB.PDBIO()
    writer.set_structure(model)
    writer.save(str(clean), SelectTarget())
    text = clean.read_text()
    # Blank altloc after resolving it, keeping residue/insertion numbering intact.
    clean.write_text('\n'.join(l[:16] + ' ' + l[17:] if l.startswith(('ATOM', 'HETATM')) else l for l in text.splitlines()) + '\n', encoding='ascii')
    pdbqt = folder / 'receptor.pdbqt'
    cmd = [sys.executable, '-B', '-m', 'meeko.cli.mk_prepare_receptor', '--read_pdb', str(clean), '-p', str(pdbqt), '--default_altloc', 'A']
    proc = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
    (folder / 'meeko_receptor.log').write_text(proc.stdout + '\n' + proc.stderr, encoding='utf-8')
    if proc.returncode or not pdbqt.is_file():
        raise ValueError('Meeko receptor preparation failed; inspect receptor/meeko_receptor.log. No permissive residue deletion was enabled.')
    prepared = pdbqt.read_text().splitlines()
    observed = {(l[21], int(l[22:26]), l[26]) for l in prepared if l.startswith(('ATOM', 'HETATM'))}
    missing = [list(r) for r in expected if r not in observed]
    metal_checks = []
    original_metals = [a for c in model for r in c for a in r if (c.id, r.id, a.name) in retained_metals]
    qt_xyz = coordinates(pdbqt)
    from scipy.spatial import cKDTree
    tree = cKDTree(qt_xyz)
    heavy_atoms = [a for r in protein_residues(protein) for a in r if a.element not in {'H', 'D'}]
    missing_atoms = [f'{a.parent.resname}:{a.parent.id}:{a.name}' for a in heavy_atoms
                     if tree.query(a.coord)[0] > 0.02]
    for atom in original_metals:
        present = bool(np.any(np.linalg.norm(qt_xyz - atom.coord, axis=1) < 0.01))
        metal_checks.append({'element': atom.element, 'xyz': atom.coord.tolist(), 'present': present})
    details = {**meta, 'expected_residues': len(expected), 'missing_prepared_residues': missing,
               'missing_prepared_heavy_atoms': missing_atoms, 'metals': metal_checks,
               'command': cmd, 'preparation': 'Meeko templates; no missing-residue modelling or pH enumeration'}
    save_json(folder / 'receptor_metadata.json', details)
    if missing or missing_atoms or any(not m['present'] for m in metal_checks):
        raise ValueError('Prepared receptor lost residues or retained metals; see receptor_metadata.json.')
    return pdbqt


def prepare_ligand(compound, ligand_input, job, seed):
    folder = job / 'ligand'
    folder.mkdir()
    mol = None
    original = None
    # Explicit input takes precedence; PubChem identity must not silently override it.
    source = ligand_input or (compound.sdf_path if compound and compound.sdf_path else None)
    if source:
        try:
            is_file = Path(str(source)).is_file()
        except OSError:
            is_file = False
        if is_file:
            original = folder / 'ligand_original.sdf'
            shutil.copy2(source, original)
            mol = next((m for m in Chem.SDMolSupplier(str(original), removeHs=False) if m is not None), None)
        else:
            mol = Chem.MolFromSmiles(str(source))
    if mol is None and compound and not ligand_input:
        mol = Chem.MolFromSmiles(compound.isomeric_smiles or compound.canonical_smiles)
    if mol is None:
        raise ValueError('Cannot parse ligand structure.')
    before = Chem.MolToSmiles(Chem.RemoveHs(mol), isomericSmiles=True)
    mol = rdMolStandardize.LargestFragmentChooser().choose(mol)
    after = Chem.MolToSmiles(Chem.RemoveHs(mol), isomericSmiles=True)
    mol = Chem.AddHs(mol, addCoords=True)
    had_3d = bool(mol.GetNumConformers() and mol.GetConformer().Is3D())
    if had_3d:
        ids = [mol.GetConformer().GetId()]
    else:
        mol.RemoveAllConformers()
        params = AllChem.ETKDGv3()
        params.randomSeed = seed
        ids = list(AllChem.EmbedMultipleConfs(mol, numConfs=20, params=params))
    if not ids:
        raise ValueError('Ligand conformer generation failed.')
    props = AllChem.MMFFGetMoleculeProperties(mol)
    energies = []
    for cid in ids:
        ff = AllChem.MMFFGetMoleculeForceField(mol, props, confId=int(cid)) if props else AllChem.UFFGetMoleculeForceField(mol, confId=int(cid))
        if ff is None:
            raise ValueError('No supported forcefield for ligand optimisation.')
        convergence = ff.Minimize(maxIts=2000)
        energies.append((float(ff.CalcEnergy()), int(cid), convergence))
    best = min(energies)
    chosen = Chem.Conformer(mol.GetConformer(best[1]))
    mol.RemoveAllConformers()
    mol.AddConformer(chosen, assignId=True)
    sdf = folder / 'ligand_prepared.sdf'
    with Chem.SDWriter(str(sdf)) as writer:
        writer.write(mol)
    pdbqt = folder / 'ligand.pdbqt'
    cmd = [sys.executable, '-B', '-m', 'meeko.cli.mk_prepare_ligand', '-i', str(sdf), '-o', str(pdbqt)]
    proc = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
    (folder / 'meeko_ligand.log').write_text(proc.stdout + '\n' + proc.stderr, encoding='utf-8')
    save_json(folder / 'ligand_metadata.json', {'query': compound.query if compound else str(ligand_input),
              'pubchem_cid': compound.cid if compound else None, 'inchikey': compound.inchikey if compound else None,
              'smiles_before': before, 'smiles_after': after, 'formal_charge': Chem.GetFormalCharge(mol),
              'stereocenters': Chem.FindMolChiralCenters(mol, includeUnassigned=True), 'had_3d': had_3d,
              'forcefield': 'MMFF94' if props else 'UFF', 'optimised_conformers': len(energies),
              'selected_energy': best[0], 'converged': best[2] == 0, 'seed': seed,
              'ph_microstates_enumerated': False, 'command': cmd})
    if proc.returncode or not pdbqt.is_file():
        raise ValueError('Meeko ligand preparation failed; see meeko_ligand.log.')
    return pdbqt


def coordinates(path):
    xyz = [[float(l[30:38]), float(l[38:46]), float(l[46:54])] for l in Path(path).read_text().splitlines() if l.startswith(('ATOM', 'HETATM'))]
    if not xyz or not np.isfinite(xyz).all():
        raise ValueError('Receptor has no valid finite atom coordinates.')
    return np.array(xyz)


def local_pocket(model, chain, gene, meta, job):
    import dock_pipeline as sources
    protein_xyz = np.array([a.coord for r in protein_residues(model[chain]) for a in r])
    class LocalContext(PDB.Select):
        def accept_residue(self, residue):
            if residue.id[0] == ' ':
                return residue.parent.id == chain and PDB.is_aa(residue, standard=False)
            return any(np.min(np.linalg.norm(protein_xyz - atom.coord, axis=1)) <= 4.5 for atom in residue)
    context = job / 'local_pocket_context.pdb'
    writer = PDB.PDBIO()
    writer.set_structure(model)
    writer.save(str(context), LocalContext())
    bound = sources.choose_cocrystal_pocket(context)
    if bound:
        return (*bound[:2], 'co_crystal_ligand')
    if gene and meta.get('alignment'):
        positions = set(sources.feature_positions(gene.features))
        ids = {(m['pdb_number'], m['insertion']) for m in meta['alignment']['mapping'] if m['uniprot'] in positions}
        atoms = [{'x': float(a.coord[0]), 'y': float(a.coord[1]), 'z': float(a.coord[2])}
                 for r in protein_residues(model[chain]) if (r.id[1], r.id[2].strip()) in ids for a in r]
        if atoms:
            return (*sources.box_from_atoms(atoms), 'mapped_uniprot_functional_site')
    metal = sources.choose_metal_site_pocket(context)
    return metal or sources.run_p2rank(job / 'receptor' / 'receptor_clean.pdb', job)


def whole_box(path, padding=5.0, center=None, size=None):
    xyz = coordinates(path)
    if not math.isfinite(padding) or padding < 0:
        raise ValueError('Padding must be finite and nonnegative.')
    low, high = xyz.min(axis=0), xyz.max(axis=0)
    if (center is None) != (size is None):
        raise ValueError('Specify both --center and --size.')
    center = np.round(np.array(center if center is not None else (low + high) / 2, dtype=float), 3)
    size = np.array(size if size is not None else 2 * np.max(np.abs(xyz - center), axis=0) + 2 * padding, dtype=float)
    center, size = np.round(center, 3), np.ceil(size * 1000) / 1000
    if not np.isfinite(center).all() or not np.isfinite(size).all() or np.any(size <= 0):
        raise ValueError('Grid must have finite coordinates and positive dimensions.')
    if np.any(xyz < center - size / 2 - 1e-9) or np.any(xyz > center + size / 2 + 1e-9):
        raise ValueError('Grid does not enclose every prepared receptor atom.')
    return center.tolist(), size.tolist(), {'atom_count': len(xyz), 'min': low.tolist(), 'max': high.tolist(), 'all_atoms_inside': True}
