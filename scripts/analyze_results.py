"""Group Vina scores by target; never assign universal affinity thresholds."""
from __future__ import annotations
import argparse
import csv
import json
import math
from pathlib import Path
from workflow import ROOT, d_path, write_results


def load_results(path):
    path = Path(path)
    if path.is_dir():
        files = sorted(path.glob('*/result.json'))
        if not files:
            files = sorted(path.glob('work/*/result.json'))
        return [json.loads(p.read_text(encoding='utf-8')) for p in files]
    if not path.exists():
        return []
    if path.suffix.lower() == '.json':
        data = json.loads(path.read_text(encoding='utf-8'))
        return [data] if isinstance(data, dict) else data
    with path.open(encoding='utf-8-sig', newline='') as handle:
        return list(csv.DictReader(handle))


def analyze(rows, export_csv=None):
    groups, failed, seen = {}, [], set()
    for row in rows:
        if row.get('status') == 'REUSED':
            continue
        try:
            score = float(row['best_binding_energy_kcal_mol'])
            valid = row.get('status') == 'SUCCESS' and math.isfinite(score)
        except (ValueError, KeyError, TypeError):
            valid = False
        if not valid:
            failed.append(row)
            continue
        key = row.get('fingerprint') or row.get('work_dir') or row.get('job_id')
        if key in seen:
            continue
        seen.add(key)
        target = row.get('gene_name') or row.get('protein_target') or 'UNKNOWN'
        groups.setdefault(target, []).append((score, row))
    ranked = []
    for target, group in groups.items():
        print(f'\nTarget: {target} | Vina score (kcal/mol)')
        for rank, (score, row) in enumerate(sorted(group, key=lambda x: x[0]), 1):
            print(f'{rank:3}  {score:8.2f}  {row.get("compound_name", "")}  render={row.get("render_status", "unknown")}')
            ranked.append(row)
    for row in failed:
        print(f'{row.get("job_id")}: {row.get("status")} {row.get("message", "")}')
    if not rows:
        print('No results found.')
    if export_csv:
        write_results(ranked, d_path(export_csv))
    return ranked


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('path', nargs='?', type=Path, default=ROOT / 'work')
    parser.add_argument('--export', type=Path)
    args = parser.parse_args()
    analyze(load_results(args.path), args.export)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
