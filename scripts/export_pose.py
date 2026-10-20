"""Export Meeko atom-mapped poses without guessing elements or chemical bonds."""
from __future__ import annotations
import argparse
import io
import json
from pathlib import Path
from rdkit import Chem
from meeko import PDBQTMolecule, RDKitMolCreate


def pose_molecule(job_dir, model_idx=1):
    if model_idx < 1:
        raise ValueError('Pose index must be positive.')
    docked = PDBQTMolecule.from_file(str(Path(job_dir) / 'vina_out.pdbqt'), skip_typing=True)
    sdf, failures = RDKitMolCreate.write_sd_string(docked)
    if failures:
        raise ValueError(f'Meeko could not reconstruct atom mapping: {failures}')
    molecules = list(Chem.ForwardSDMolSupplier(io.BytesIO(sdf.encode()), removeHs=False))
    if model_idx > len(molecules) or molecules[model_idx - 1] is None:
        raise ValueError(f'Pose {model_idx} not present.')
    mol = molecules[model_idx - 1]
    original = next(m for m in Chem.SDMolSupplier(str(Path(job_dir) / 'ligand' / 'ligand_prepared.sdf'), removeHs=False) if m)
    if Chem.MolToSmiles(Chem.RemoveHs(mol), isomericSmiles=True) != Chem.MolToSmiles(Chem.RemoveHs(original), isomericSmiles=True):
        raise ValueError('Export changed ligand identity, bond orders, or stereochemistry.')
    return mol


def export_pose(job_dir, model_idx=1, out_pdb=None, make_complex=True):
    from workflow import d_path
    job_dir = d_path(job_dir)
    mol = pose_molecule(job_dir, model_idx)
    out = d_path(out_pdb or job_dir / f'pose_rank{model_idx}.pdb')
    out.parent.mkdir(parents=True, exist_ok=True)
    with Chem.SDWriter(str(job_dir / f'pose_rank{model_idx}.sdf')) as writer:
        writer.write(mol)
    receptor = job_dir / 'receptor' / 'receptor_clean.pdb'
    rec_lines = [l for l in receptor.read_text().splitlines() if l.startswith(('ATOM', 'HETATM', 'TER'))]
    chains = {l[21] for l in rec_lines if l.startswith(('ATOM', 'HETATM'))}
    chain = next(c for c in 'ZYXWVUTSRQPONMLKJIHGFEDCBA0123456789' if c not in chains)
    maximum = max(int(l[6:11]) for l in rec_lines if l[6:11].strip().isdigit())
    if maximum + mol.GetNumAtoms() > 99999:
        raise ValueError('PDB serial range exceeded; mmCIF complex support required.')
    for i, atom in enumerate(mol.GetAtoms(), 1):
        info = Chem.AtomPDBResidueInfo()
        info.SetName(f'{atom.GetSymbol()}{i}'.rjust(4))
        info.SetResidueName('LIG')
        info.SetResidueNumber(1)
        info.SetChainId(chain)
        info.SetIsHeteroAtom(True)
        atom.SetMonomerInfo(info)
    ligand = []
    for line in Chem.MolToPDBBlock(mol).splitlines():
        if line.startswith(('ATOM', 'HETATM')):
            ligand.append(f'HETATM{maximum + int(line[6:11]):5d}' + line[11:])
        elif line.startswith('CONECT'):
            ids = [int(line[k:k+5]) for k in range(6, len(line), 5) if line[k:k+5].strip()]
            ligand.append('CONECT' + ''.join(f'{maximum + n:5d}' for n in ids))
    out.write_text('\n'.join(ligand) + '\nEND\n', encoding='ascii')
    if make_complex:
        (job_dir / f'complex_rank{model_idx}.pdb').write_text('\n'.join(rec_lines + ligand) + '\nEND\n', encoding='ascii')
    (job_dir / 'pose_export.json').write_text(json.dumps({'model': model_idx, 'ligand_chain': chain, 'residue': 'LIG', 'residue_number': 1, 'atom_mapping': 'Meeko REMARK SMILES/SMILES IDX/H PARENT'}, indent=2))
    return 0


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('job_dir', type=Path)
    parser.add_argument('--model', type=int, default=None)
    parser.add_argument('--out', type=Path)
    parser.add_argument('--no-complex', action='store_true')
    args = parser.parse_args()
    result_path = args.job_dir / 'result.json'
    best = json.loads(result_path.read_text(encoding='utf-8')).get('best_model', 1) if result_path.is_file() else 1
    return export_pose(args.job_dir, args.model or best, args.out, not args.no_complex)


if __name__ == '__main__':
    raise SystemExit(main())
