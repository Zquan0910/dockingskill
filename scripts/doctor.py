from __future__ import annotations

import warnings
warnings.filterwarnings("ignore")
import argparse
import importlib
import json
import os
import shutil
import subprocess
import sys
import urllib.request
from pathlib import Path
from typing import Any

USER_AGENT = "DockingSkill-Doctor/1.0"

# Keep diagnostics independent of scientific imports so missing libraries can be reported.
TEMP = Path('D:/ChatGPT/geneprotein-dock-runs/tmp')
TEMP.mkdir(parents=True, exist_ok=True)
for key in ('TEMP', 'TMP', 'TMPDIR'):
    os.environ[key] = str(TEMP)
os.environ['PYTHONDONTWRITEBYTECODE'] = '1'
os.environ['MPLCONFIGDIR'] = str(TEMP.parent / 'mplconfig')
os.environ['PATH'] = 'D:/GeneProteinDocking/env/Library/bin;' + os.environ.get('PATH', '')


def check_python() -> dict[str, Any]:
    return {
        "name": "python",
        "version": sys.version.split()[0],
        "executable": sys.executable,
        "ok": sys.version_info >= (3, 10),
    }


def check_module(name: str, hint: str = "") -> dict[str, Any]:
    try:
        mod = importlib.import_module(name)
        ver = getattr(mod, "__version__", "installed")
        return {"name": name, "ok": True, "version": str(ver)}
    except Exception as exc:
        return {"name": name, "ok": False, "error": str(exc), "hint": hint}


def check_vina() -> dict[str, Any]:
    from dock_pipeline import read_config
    candidates = [
        Path("D:/Download/vina.exe"),
        Path("D:/GeneProteinDocking/tools/vina/vina.exe"),
    ]
    try:
        configured = read_config().get('software', {}).get('vina_executable')
    except Exception as exc:
        return {'name': 'vina', 'ok': False, 'error': f'Cannot read config.yaml: {exc}'}
    if configured:
        candidates.insert(0, Path(configured))
    which = shutil.which("vina") or shutil.which("vina.exe")
    if which:
        candidates.append(Path(which))

    for exe in candidates:
        if exe.exists():
            try:
                proc = subprocess.run([str(exe), "--version"], capture_output=True, text=True, timeout=10)
                output = (proc.stdout + proc.stderr).strip()
                if proc.returncode == 0 or "AutoDock Vina" in output:
                    return {"name": "vina", "ok": True, "path": str(exe), "version": output.splitlines()[0]}
            except Exception:
                continue
    return {"name": "vina", "ok": False, "error": "vina.exe not found or not executable"}


def check_p2rank() -> dict[str, Any]:
    bat_path = Path("D:/GeneProteinDocking/tools/p2rank_run.bat")
    if bat_path.exists():
        try:
            proc = subprocess.run([str(bat_path), "--version"], capture_output=True, text=True, timeout=15)
            out = (proc.stdout + proc.stderr).strip()
            if proc.returncode == 0 or "P2Rank" in out:
                return {"name": "p2rank", "ok": True, "path": str(bat_path), "version": out.splitlines()[0]}
        except Exception as exc:
            return {"name": "p2rank", "ok": False, "error": str(exc)}
    return {"name": "p2rank", "ok": False, "error": "p2rank launcher not found at D:/GeneProteinDocking/tools/p2rank_run.bat"}


def check_network(name: str, url: str, method: str = "GET", payload: dict | None = None) -> dict[str, Any]:
    try:
        data = json.dumps(payload).encode("utf-8") if payload else None
        headers = {"User-Agent": USER_AGENT}
        if payload:
            headers["Content-Type"] = "application/json"
        req = urllib.request.Request(url, data=data, headers=headers, method=method)
        with urllib.request.urlopen(req, timeout=15) as resp:
            return {"name": name, "ok": 200 <= resp.status < 300, "status": resp.status}
    except Exception as exc:
        return {"name": name, "ok": False, "error": str(exc)}


def run_diagnostics(local_only=False) -> dict[str, Any]:
    modules = [
        ("rdkit", "conda install -c conda-forge rdkit"),
        ("meeko", "pip install meeko"),
        ("gemmi", "conda install -c conda-forge gemmi"),
        ("Bio", "conda install -c conda-forge biopython"),
        ("pandas", "conda install pandas"),
        ("numpy", "conda install numpy"),
        ("scipy", "conda install scipy"),
        ("yaml", "conda install pyyaml"),
        ("openbabel", "conda install -c conda-forge openbabel"),
        ("pymol", "Install pymol-open-source into the D: docking environment"),
        ("plip", "Install PLIP into the D: docking environment"),
        ("lxml", "Install lxml into the D: docking environment"),
        ("PIL", "Install Pillow into the D: docking environment"),
    ]

    report: dict[str, Any] = {
        "python": check_python(),
        "tools": [check_vina(), check_p2rank()],
        "modules": [check_module(m, h) for m, h in modules],
        "network": [] if local_only else [
            check_network("UniProt", "https://rest.uniprot.org/uniprotkb/search?query=gene_exact:YES1%20AND%20organism_id:9606&format=json&size=1"),
            check_network("RCSB_REST", "https://data.rcsb.org/rest/v1/core/entry/4hhb"),
            check_network(
                "RCSB_POST",
                "https://search.rcsb.org/rcsbsearch/v2/query",
                method="POST",
                payload={
                    "query": {
                        "type": "terminal",
                        "service": "text",
                        "parameters": {
                            "attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession",
                            "operator": "exact_match",
                            "value": "P07947",
                        },
                    },
                    "return_type": "entry",
                    "request_options": {"paginate": {"start": 0, "rows": 1}},
                },
            ),
            check_network("PubChem", "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/quercetin/cids/JSON"),
            check_network("AlphaFold_DB", "https://alphafold.ebi.ac.uk/api/prediction/P07947"),
        ],
    }

    report["all_tools_ok"] = all(t["ok"] for t in report["tools"] if t['name'] != 'p2rank')
    report['network_skipped'] = local_only
    report['notes'] = ['P2Rank/Java are optional for explicit local-pocket mode. Rendering requires PyMOL and PLIP.']
    report["all_modules_ok"] = all(m["ok"] for m in report["modules"])
    report["all_network_ok"] = all(n["ok"] for n in report["network"])
    report["healthy"] = report["python"]["ok"] and report["all_tools_ok"] and report["all_modules_ok"]
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description='Check local docking and rendering tools before starting a job.')
    parser.add_argument('--local-only', action='store_true', help='Check tools/libraries without contacting databases')
    args = parser.parse_args()
    report = run_diagnostics(local_only=args.local_only)
    print("=" * 60)
    print("DockingSkill Doctor Diagnostics")
    print("=" * 60)
    py = report["python"]
    print(f"Python:    {'OK' if py['ok'] else 'FAIL'} (v{py['version']}) -> {py['executable']}")
    print("-" * 60)
    print("External Tools:")
    for t in report["tools"]:
        status = "OK" if t["ok"] else "FAIL"
        extra = t.get("version") or t.get("error", "")
        print(f"  [{status}] {t['name']:<10}: {extra}")
    print("-" * 60)
    print("Python Libraries:")
    for m in report["modules"]:
        status = "OK" if m["ok"] else "FAIL"
        extra = m.get("version") if m["ok"] else m.get("error", "")
        print(f"  [{status}] {m['name']:<10}: {extra}")
    print("-" * 60)
    print("Network Services:")
    if report['network_skipped']:
        print('  SKIPPED (--local-only); online structure retrieval still requires working network services.')
    for n in report["network"]:
        status = "OK" if n["ok"] else "FAIL"
        extra = f"HTTP {n.get('status')}" if n["ok"] else n.get("error", "")
        print(f"  [{status}] {n['name']:<12}: {extra}")
    print("=" * 60)
    if report["healthy"]:
        print("Local toolchain: HEALTHY. Ready for molecular docking.")
        if not report['network_skipped'] and not report['all_network_ok']:
            print('Network warning: some databases failed; local inputs may still be usable.')
        return 0
    else:
        print("System status: DEGRADED. Check failed items above.")
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
