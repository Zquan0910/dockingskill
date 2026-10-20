"""Validated local docking workflow. Retrieval/local-pocket helpers live in dock_pipeline."""
from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.metadata
import json
import math
import os
import shutil
import subprocess
import sys
import traceback
import uuid
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path('D:/ChatGPT/geneprotein-dock-runs')


def d_path(path):
    p = Path(path).resolve()
    if p.drive.upper() != 'D:':
        raise ValueError(f'Output must be on D: {p}')
    return p


def runtime():
    temp = ROOT / 'tmp'
    temp.mkdir(parents=True, exist_ok=True)
    for key in ('TEMP', 'TMP', 'TMPDIR'):
        os.environ[key] = str(temp)
    os.environ['PYTHONDONTWRITEBYTECODE'] = '1'
    os.environ['MPLCONFIGDIR'] = str(ROOT / 'mplconfig')
    os.environ['PATH'] = 'D:/GeneProteinDocking/env/Library/bin;' + os.environ.get('PATH', '')


runtime()
import dock_pipeline as sources
from preparation import save_json, select_structure, prepare_receptor, prepare_ligand, whole_box, local_pocket


def settings(overrides):
    defaults = dict(exhaustiveness=64, num_modes=10, energy_range=5, seed=20260715,
                    cpu=0, timeout=7200, mode='whole', padding=5.0)
    config = sources.read_config()
    defaults.update({k: v for k, v in config.get('docking', {}).items() if k in defaults})
    defaults.update({k: v for k, v in overrides.items() if k in defaults and v is not None})
    for key in ('exhaustiveness', 'num_modes', 'seed', 'cpu', 'timeout'):
        defaults[key] = int(defaults[key])
    for key in ('energy_range', 'padding'):
        defaults[key] = float(defaults[key])
    if defaults['exhaustiveness'] < 1 or defaults['num_modes'] < 1 or defaults['timeout'] < 1 or defaults['cpu'] < 0:
        raise ValueError('Invalid Vina integer parameter.')
    if any(not math.isfinite(defaults[k]) or defaults[k] < 0 for k in ('energy_range', 'padding')):
        raise ValueError('Energy range and padding must be finite and nonnegative.')
    if defaults['mode'] not in ('whole', 'local'):
        raise ValueError('Mode must be whole or local.')
    return defaults


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def fingerprint(receptor, ligand, center, size, params, vina):
    versions = {n: importlib.metadata.version(n) for n in ('meeko', 'rdkit', 'biopython')}
    scripts = {name: sha(Path(__file__).parent / name) for name in ('workflow.py', 'preparation.py', 'dock_pipeline.py')}
    data = {'receptor_sha256': sha(receptor), 'ligand_sha256': sha(ligand), 'center': center,
            'size': size, 'parameters': params, 'vina_sha256': sha(vina), 'versions': versions, 'scripts': scripts}
    return hashlib.sha256(json.dumps(data, sort_keys=True).encode()).hexdigest(), data


def run_vina(vina, receptor, ligand, center, size, job, params):
    out = job / 'vina_out.pdbqt'
    args = [str(vina), '--receptor', str(receptor), '--ligand', str(ligand)]
    for name, vals in [('center', center), ('size', size)]:
        for axis, val in zip('xyz', vals):
            args += [f'--{name}_{axis}', f'{val:.3f}']
    for key in ('exhaustiveness', 'num_modes', 'energy_range', 'seed', 'cpu'):
        args += ['--' + key, str(params[key])]
    args += ['--out', str(out)]
    save_json(job / 'vina_run.json', {'command': args, 'timeout': params['timeout'], 'parameters': params})
    print(f'Vina: grid={size}; exhaustiveness={params["exhaustiveness"]}; job={job}', flush=True)
    with (job / 'vina_stdout.txt').open('w', encoding='utf-8') as stdout, (job / 'vina_stderr.txt').open('w', encoding='utf-8') as stderr:
        proc = subprocess.Popen(args, cwd=job, stdout=stdout, stderr=stderr)
        try:
            code = proc.wait(timeout=params['timeout'])
        except subprocess.TimeoutExpired:
            proc.kill()
            proc.wait()
            raise RuntimeError('Vina timed out. Logs retained; no box/search reduction attempted.')
    if code or not out.is_file():
        raise RuntimeError(f'Vina failed ({code}); see vina_stderr.txt.')
    energies = []
    for line in out.read_text().splitlines():
        if line.startswith('REMARK VINA RESULT:'):
            energies.append(float(line.split(':', 1)[1].split()[0]))
    if not energies or not all(math.isfinite(e) for e in energies):
        raise RuntimeError('Vina produced no valid scored poses.')
    return energies


def render_job(result):
    from render_results import render_one
    try:
        outputs = render_one(Path(result['work_dir']))
        result.update(render_status='SUCCESS', figures=outputs)
    except Exception as exc:
        result.update(render_status='FAILED', render_error=str(exc))
        (Path(result['work_dir']) / 'render_error.log').write_text(traceback.format_exc(), encoding='utf-8')
    save_json(Path(result['work_dir']) / 'result.json', result)


def run_dock_job(job_id=None, gene_name=None, compound_name=None, protein_path=None, ligand_input=None,
                 custom_center=None, custom_size=None, work_dir=None, cache_dir=None,
                 exhaustiveness=None, chain=None, render=None, resume=True, prepare_only=False, **overrides):
    work = d_path(work_dir or ROOT / 'work')
    cache = d_path(cache_dir or ROOT / 'cache')
    identifier = sources.slug(job_id or f'{gene_name or Path(str(protein_path)).stem}_{compound_name or "ligand"}')
    stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S')
    job = work / f'{identifier}_{stamp}_{uuid.uuid4().hex[:8]}'
    job.mkdir(parents=True)
    result = {'job_id': job_id or identifier, 'gene_name': gene_name or '', 'compound_name': overrides.get('compound_label') or compound_name or str(ligand_input or ''),
              'work_dir': str(job), 'status': 'PREPARING', 'created_at': stamp}
    save_json(job / 'result.json', result)
    try:
        params = settings(dict(overrides, exhaustiveness=exhaustiveness))
        gene = sources.resolve_gene(gene_name, cache) if gene_name else None
        compound = sources.resolve_compound(compound_name, cache) if compound_name and not ligand_input else None
        model, selected, meta = select_structure(gene, job, cache,
                float(sources.read_config().get('structure', {}).get('minimum_coverage', 0.95)), protein_path, chain)
        receptor = prepare_receptor(model, selected, job, meta)
        ligand = prepare_ligand(compound, ligand_input, job, params['seed'])
        if params['mode'] == 'whole':
            center, size, bounds = whole_box(receptor, params['padding'], custom_center, custom_size)
            pocket = 'whole_receptor'
        else:
            if (custom_center is None) != (custom_size is None):
                raise ValueError('Specify both center and size for a manual box.')
            local = None
            if custom_center is not None:
                local = (custom_center, custom_size, 'user_specified')
            else:
                local = local_pocket(model, selected, gene, meta, job)
            if not local:
                raise ValueError('No local pocket; provide a manual local box.')
            center, size, pocket = local
            if not all(math.isfinite(v) for v in (*center, *size)) or any(v <= 0 for v in size):
                raise ValueError('Invalid local grid.')
            center, size, bounds = list(center), list(size), {'all_atoms_inside': False, 'mode': 'explicit_local'}
        save_json(job / 'grid.json', {'center': center, 'size': size, 'padding_per_side': params['padding'], 'source': pocket, **bounds})
        vina = sources.detect_vina()
        key, details = fingerprint(receptor, ligand, center, size, params, vina)
        save_json(job / 'fingerprint.json', {'sha256': key, **details})
        result.update(protein_target=meta.get('accession') or meta.get('pdb_id'), pubchem_cid=compound.cid if compound else '',
                      structure_source=meta['source'], target_chain=selected, pocket_source=pocket, grid_center=center,
                      grid_size=size, parameters=params, fingerprint=key, render_status='PENDING')
        do_render = render if render is not None else sources.read_config().get('render', {}).get('enabled', True)
        if prepare_only:
            result.update(status='PREPARED', render_status='NOT_RUN')
            save_json(job / 'result.json', result)
            return result
        if resume:
            for old in sorted(work.glob('*/result.json'), reverse=True):
                if old.parent == job:
                    continue
                try:
                    previous = json.loads(old.read_text(encoding='utf-8'))
                    pose = old.parent / 'vina_out.pdbqt'
                    if previous.get('fingerprint') == key and previous.get('status') == 'SUCCESS' and pose.is_file() and previous.get('pose_sha256') == sha(pose):
                        for name in ('vina_out.pdbqt', 'vina_run.json', 'vina_stdout.txt', 'vina_stderr.txt'):
                            if (old.parent / name).is_file():
                                shutil.copy2(old.parent / name, job / name)
                        for name in ('vina_scores_kcal_mol', 'best_model', 'best_binding_energy_kcal_mol', 'pose_sha256'):
                            result[name] = previous[name]
                        result.update(status='SUCCESS', reused=True, reused_from=previous['work_dir'], render_status='NOT_RUN')
                        save_json(job / 'result.json', result)
                        if do_render:
                            render_job(result)
                        return result
                except (OSError, ValueError):
                    continue
        energies = run_vina(vina, receptor, ligand, center, size, job, params)
        result.update(status='SUCCESS', vina_scores_kcal_mol=energies, best_model=energies.index(min(energies)) + 1,
                      best_binding_energy_kcal_mol=f'{min(energies):.2f}', pose_sha256=sha(job / 'vina_out.pdbqt'))
        save_json(job / 'result.json', result)
        if do_render:
            render_job(result)
        else:
            result['render_status'] = 'NOT_RUN'
            save_json(job / 'result.json', result)
    except Exception as exc:
        result.update(status='FAILED', message=str(exc), render_status='NOT_RUN')
        (job / 'error.log').write_text(traceback.format_exc(), encoding='utf-8')
        save_json(job / 'result.json', result)
    return result


def write_results(rows, path):
    fields = ['job_id', 'gene_name', 'protein_target', 'compound_name', 'pubchem_cid', 'pocket_source',
              'best_binding_energy_kcal_mol', 'status', 'render_status', 'message', 'work_dir', 'fingerprint']
    path = d_path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', encoding='utf-8-sig', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)


def run_batch_csv(csv_path, out_csv=None, **options):
    csv_path = Path(csv_path).resolve()
    with csv_path.open(encoding='utf-8-sig', newline='') as handle:
        jobs = list(csv.DictReader(handle))
    batch_dir = ROOT / 'batches' / (datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S') + '_' + uuid.uuid4().hex[:8])
    batch_dir.mkdir(parents=True)
    rows = []
    for i, item in enumerate(jobs, 1):
        compound = f"CID:{item['pubchem_cid']}" if item.get('pubchem_cid') else item.get('compound_name') or item.get('compound')
        protein = item.get('protein_path') or item.get('protein')
        ligand = item.get('ligand_input') or item.get('ligand')
        for key, value in [('protein', protein), ('ligand', ligand)]:
            if value and (csv_path.parent / value).is_file():
                if key == 'protein': protein = str(csv_path.parent / value)
                else: ligand = str(csv_path.parent / value)
        job_options = dict(options, compound_label=item.get('compound_name') or item.get('compound') or compound)
        if item.get('chain'):
            job_options['chain'] = item['chain']
        row = run_dock_job(job_id=item.get('job_id') or f'JOB_{i:04}', gene_name=item.get('gene_name') or item.get('gene'),
                           compound_name=compound, protein_path=protein, ligand_input=ligand, **job_options)
        rows.append(row)
        write_results(rows, out_csv or batch_dir / 'results.csv')
        print(f'{i}/{len(jobs)}: {row["status"]} {row.get("best_binding_energy_kcal_mol", "")}', flush=True)
    if rows and all(r.get('render_status') == 'SUCCESS' for r in rows):
        from render_results import assemble_summary
        assemble_summary([Path(r['work_dir']) / 'figures' / 'panel.png' for r in rows], batch_dir)
    return 0 if rows and all(r['status'] == 'SUCCESS' and r.get('render_status') != 'FAILED' for r in rows) else 1


def main(argv=None):
    parser = argparse.ArgumentParser(description='Whole-protein Vina docking and reference-style PyMOL rendering')
    parser.add_argument('gene', nargs='?')
    parser.add_argument('compound', nargs='?')
    parser.add_argument('--gene', dest='gene_opt')
    parser.add_argument('--compound', dest='compound_opt')
    parser.add_argument('--protein', help='PDB/mmCIF path or PDB ID; overrides automatic structure selection')
    parser.add_argument('--ligand', help='SMILES or SDF path; overrides named compound structure')
    parser.add_argument('--chain')
    parser.add_argument('--batch', type=Path)
    parser.add_argument('--job-id')
    parser.add_argument('--out-dir', type=Path)
    parser.add_argument('--center', nargs=3, type=float)
    parser.add_argument('--size', nargs=3, type=float)
    parser.add_argument('--mode', choices=['whole', 'local'])
    parser.add_argument('--padding', type=float)
    for key in ('exhaustiveness', 'num_modes', 'seed', 'cpu', 'timeout'):
        parser.add_argument('--' + key.replace('_', '-'), type=int, default=None)
    parser.add_argument('--energy-range', type=float)
    parser.add_argument('--no-render', action='store_true')
    parser.add_argument('--no-resume', action='store_true')
    parser.add_argument('--prepare-only', action='store_true')
    args = parser.parse_args(argv)
    options = {k: getattr(args, k) for k in ('exhaustiveness', 'num_modes', 'seed', 'cpu', 'timeout', 'energy_range', 'mode', 'padding', 'chain', 'prepare_only')}
    options.update(work_dir=args.out_dir, custom_center=args.center, custom_size=args.size,
                   render=False if args.no_render else None, resume=not args.no_resume)
    if args.batch:
        return run_batch_csv(args.batch, **options)
    gene, compound = args.gene_opt or args.gene, args.compound_opt or args.compound
    if not (gene or args.protein) or not (compound or args.ligand):
        parser.error('Provide a gene/protein and a compound/ligand.')
    result = run_dock_job(args.job_id, gene, compound, args.protein, args.ligand, **options)
    write_results([result], Path(result['work_dir']) / 'results.csv')
    print(json.dumps(result, ensure_ascii=False, indent=2), flush=True)
    return 0 if result['status'] in ('SUCCESS', 'PREPARED') and result.get('render_status') != 'FAILED' else 1


if __name__ == '__main__':
    raise SystemExit(main())
