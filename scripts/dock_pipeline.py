from __future__ import annotations

import warnings
warnings.filterwarnings("ignore")

import argparse
import csv
import hashlib
import json
import math
import os
import re
import shutil
import subprocess
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import dataclass
from pathlib import Path
from typing import Any

SCRIPT_DIR = Path(__file__).resolve().parent
SKILL_ROOT = SCRIPT_DIR.parent
DEFAULT_ROOT = Path("D:/ChatGPT/geneprotein-dock-runs")
CONFIG_PATH = SKILL_ROOT / "config.yaml"
USER_AGENT = "DockingSkill/2.0 (Windows; automated academic docking)"
DEFAULT_SEED = 20260715


class DockingError(Exception):
    code = "ERROR"
    user_message = "流程失败。"

    def __init__(self, message: str | None = None, code: str | None = None):
        super().__init__(message or self.user_message)
        if code:
            self.code = code
        self.user_message = message or self.user_message


class GeneNotResolved(DockingError):
    code = "GENE_NOT_RESOLVED"


class CompoundNotResolved(DockingError):
    code = "COMPOUND_NOT_RESOLVED"


class StructureNotFound(DockingError):
    code = "STRUCTURE_NOT_FOUND"


class PocketNotFound(DockingError):
    code = "POCKET_NOT_FOUND"


class PreparationFailed(DockingError):
    code = "PREPARATION_FAILED"


class VinaFailed(DockingError):
    code = "VINA_FAILED"


def slug(value: str) -> str:
    text = re.sub(r"[^A-Za-z0-9_.-]+", "_", value.strip())
    return text.strip("_") or hashlib.sha1(value.encode("utf-8")).hexdigest()[:12]


def read_config() -> dict[str, Any]:
    if not CONFIG_PATH.exists():
        return {}
    import yaml
    config = yaml.safe_load(CONFIG_PATH.read_text(encoding="utf-8")) or {}
    if not isinstance(config, dict):
        raise ValueError('config.yaml must contain a mapping.')
    return config


def http_get(url: str, cache_path: Path | None = None, accept: str = "application/json") -> bytes:
    if cache_path and cache_path.exists() and cache_path.stat().st_size > 0:
        return cache_path.read_bytes()
    last: Exception | None = None
    for attempt in range(3):
        req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT, "Accept": accept})
        try:
            with urllib.request.urlopen(req, timeout=45) as resp:
                if resp.status < 200 or resp.status >= 300:
                    raise DockingError(f"HTTP {resp.status}: {url}", "HTTP_ERROR")
                data = resp.read()
                if cache_path:
                    cache_path.parent.mkdir(parents=True, exist_ok=True)
                    cache_path.write_bytes(data)
                return data
        except Exception as exc:
            last = exc
            time.sleep(1.5 * (2**attempt))
    raise DockingError(f"网络请求失败：{url}；{last}", "NETWORK_ERROR")


def http_json(url: str, cache_path: Path | None = None) -> Any:
    return json.loads(http_get(url, cache_path, "application/json").decode("utf-8"))


def http_post_json(url: str, payload: dict[str, Any], cache_path: Path | None = None) -> Any:
    if cache_path and cache_path.exists() and cache_path.stat().st_size > 0:
        return json.loads(cache_path.read_text(encoding="utf-8"))
    data = json.dumps(payload).encode("utf-8")
    last: Exception | None = None
    for attempt in range(3):
        req = urllib.request.Request(
            url,
            data=data,
            headers={"User-Agent": USER_AGENT, "Accept": "application/json", "Content-Type": "application/json"},
            method="POST",
        )
        try:
            with urllib.request.urlopen(req, timeout=60) as resp:
                raw = resp.read()
                if not raw:
                    return {"result_set": []}
                if cache_path:
                    cache_path.parent.mkdir(parents=True, exist_ok=True)
                    cache_path.write_bytes(raw)
                return json.loads(raw.decode("utf-8"))
        except Exception as exc:
            last = exc
            time.sleep(1.5 * (2**attempt))
    raise DockingError(f"网络请求失败：{url}；{last}", "NETWORK_ERROR")


def detect_vina() -> Path:
    config = read_config()
    configured = config.get("software", {}).get("vina_executable")
    candidates: list[Path] = []
    if configured:
        if not Path(configured).is_file():
            raise VinaFailed(f'Configured Vina executable does not exist: {configured}')
        candidates.append(Path(configured))
    which = shutil.which("vina") or shutil.which("vina.exe")
    if which:
        candidates.append(Path(which))
    for p in [
        Path("D:/Download/vina.exe"),
        Path("D:/GeneProteinDocking/tools/vina/vina.exe"),
        Path("C:/Program Files/vina/vina.exe"),
        Path("C:/Program Files (x86)/vina/vina.exe"),
    ]:
        if p.exists():
            candidates.append(p)
    seen: set[Path] = set()
    for cand in candidates:
        try:
            exe = cand.resolve()
        except Exception:
            exe = cand
        if exe in seen or not exe.exists():
            continue
        seen.add(exe)
        try:
            proc = subprocess.run([str(exe), "--version"], text=True, capture_output=True, timeout=10)
            if proc.returncode == 0 or "AutoDock Vina" in (proc.stdout + proc.stderr):
                return exe
        except Exception:
            continue
    raise VinaFailed("未找到可用的 vina.exe。请检查安装路径并在 config.yaml 中配置。", "VINA_NOT_FOUND")


def detect_p2rank() -> str | None:
    config = read_config()
    configured = config.get("software", {}).get("p2rank_command")
    if configured and Path(configured.split()[0]).exists():
        return configured
    candidates = [
        "D:/GeneProteinDocking/tools/p2rank_run.bat",
        "D:/GeneProteinDocking/tools/p2rank_2.5/prank.bat",
    ]
    for cand in candidates:
        if Path(cand).exists():
            return str(Path(cand).resolve())
    return None


@dataclass
class GeneInfo:
    query: str
    symbol: str
    accession: str
    protein_name: str
    sequence: str
    length: int
    features: list[dict[str, Any]]
    pdb_refs: list[str]


def resolve_gene(gene_name: str, cache_dir: Path) -> GeneInfo:
    q = gene_name.strip()
    if not q:
        raise GeneNotResolved("基因名称为空。")
    encoded = urllib.parse.quote(f"(gene_exact:{q}) AND (organism_id:9606)")
    url = f"https://rest.uniprot.org/uniprotkb/search?query={encoded}&format=json&size=25"
    raw_path = cache_dir / "uniprot" / f"{slug(q)}.json"
    payload = http_json(url, raw_path)
    results = payload.get("results", [])
    exact: list[dict[str, Any]] = []
    for item in results:
        genes = item.get("genes", [])
        names: set[str] = set()
        for gene in genes:
            if gene.get("geneName", {}).get("value"):
                names.add(gene["geneName"]["value"].upper())
            for syn in gene.get("synonyms", []):
                if syn.get("value"):
                    names.add(syn["value"].upper())
        if q.upper() in names:
            exact.append(item)
    exact.sort(key=lambda item: (
        0 if (("reviewed" in item.get("entryType", "").lower()) and ("unreviewed" not in item.get("entryType", "").lower())) else 1,
        item.get("primaryAccession", ""),
    ))
    if not exact:
        raise GeneNotResolved(f"无法将基因 '{gene_name}' 解析为唯一的人源蛋白。")
    reviewed = [item for item in exact if 'unreviewed' not in item.get('entryType', '').lower()]
    if len(reviewed or exact) != 1:
        raise GeneNotResolved(f"Ambiguous gene {q}: specify an unambiguous human gene symbol.")
    best = (reviewed or exact)[0]
    accession = best.get("primaryAccession")
    sequence = best.get("sequence", {}).get("value", "")
    features = best.get("features", [])
    pdb_refs = [
        ref.get("id")
        for ref in best.get("uniProtKBCrossReferences", [])
        if ref.get("database") == "PDB" and ref.get("id")
    ]
    protein_name = best.get("proteinDescription", {}).get("recommendedName", {}).get("fullName", {}).get("value", "")
    parsed = {
        "query": q,
        "symbol": q.upper(),
        "accession": accession,
        "protein_name": protein_name,
        "length": len(sequence),
        "features": features,
        "pdb_refs": pdb_refs,
    }
    (cache_dir / "uniprot" / f"{slug(q)}.parsed.json").write_text(json.dumps(parsed, ensure_ascii=False, indent=2), encoding="utf-8")
    return GeneInfo(q, q.upper(), accession, protein_name, sequence, len(sequence), features, pdb_refs)


CHINESE_COMPOUNDS = {
    "槲皮素": "quercetin",
    "山奈酚": "kaempferol",
    "木犀草素": "luteolin",
    "芹菜素": "apigenin",
    "黄芩素": "baicalein",
    "黄芩苷": "baicalin",
    "柚皮素": "naringenin",
    "葛根素": "puerarin",
    "熊果酸": "ursolic acid",
    "齐墩果酸": "oleanolic acid",
    "咖啡酸苯乙酯": "caffeic acid phenethyl ester",
    "岩藻黄素": "fucoxanthin",
    "槲皮苷": "quercitrin",
}


@dataclass
class CompoundInfo:
    query: str
    name: str
    cid: int | None
    canonical_smiles: str
    isomeric_smiles: str
    inchikey: str
    molecular_formula: str
    molecular_weight: float | None
    charge: int | None
    sdf_path: Path | None


def resolve_compound(name_or_input: str, cache_dir: Path) -> CompoundInfo:
    query = name_or_input.strip()
    if not query:
        raise CompoundNotResolved("小分子输入为空。")

    if query.upper().startswith("CID:"):
        cid = int(query.split(":", 1)[1].strip())
        lookup_name = f"CID:{cid}"
    else:
        lookup_name = CHINESE_COMPOUNDS.get(query, query)
        encoded = urllib.parse.quote(lookup_name)
        try:
            cids = http_json(
                f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{encoded}/cids/JSON",
                cache_dir / "pubchem" / f"name_{slug(lookup_name)}_cids.json",
            )
        except Exception:
            raise CompoundNotResolved(f"PubChem 未能检索到化合物 '{query}'。")
        cid_list = cids.get("IdentifierList", {}).get("CID", [])
        unique = list(dict.fromkeys(cid_list))
        if not unique:
            raise CompoundNotResolved(f"PubChem 未能检索到化合物 '{query}'。")
        if len(unique) > 1:
            raise CompoundNotResolved(f"化合物名称 '{query}' 对应多个不同化学结构（CID: {unique[:3]}...），请明确指定 PubChem CID。")
        cid = int(unique[0])

    props = "MolecularFormula,MolecularWeight,CanonicalSMILES,IsomericSMILES,InChIKey,IUPACName,Charge"
    data = http_json(
        f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/{props}/JSON",
        cache_dir / "pubchem" / f"cid_{cid}_props.json",
    )
    rows = data.get("PropertyTable", {}).get("Properties", [])
    if not rows:
        raise CompoundNotResolved(f"未能获取 CID {cid} 的分子属性。")
    row = rows[0]
    smiles = row.get("IsomericSMILES") or row.get("SMILES") or row.get("CanonicalSMILES") or row.get("ConnectivitySMILES", "")

    sdf_path = cache_dir / "pubchem" / f"cid_{cid}_3d.sdf"
    if not sdf_path.exists() or sdf_path.stat().st_size == 0:
        sdf_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/SDF?record_type=3d"
        try:
            http_get(sdf_url, sdf_path, "chemical/x-mdl-sdfile")
        except Exception:
            sdf_path.unlink(missing_ok=True)
            sdf_path = None

    info = CompoundInfo(
        query=query,
        name=row.get("IUPACName") or lookup_name,
        cid=cid,
        canonical_smiles=row.get("CanonicalSMILES", smiles),
        isomeric_smiles=smiles,
        inchikey=row.get("InChIKey", ""),
        molecular_formula=row.get("MolecularFormula", ""),
        molecular_weight=float(row["MolecularWeight"]) if row.get("MolecularWeight") else None,
        charge=int(row["Charge"]) if row.get("Charge") is not None else None,
        sdf_path=sdf_path if (sdf_path and sdf_path.exists()) else None,
    )
    return info


@dataclass
class ReceptorChoice:
    source: str
    pdb_id: str | None
    accession: str
    structure_path: Path
    pocket_source: str
    center: tuple[float, float, float]
    size: tuple[float, float, float]
    chain_ids: list[str]


COMMON_NON_LIGANDS = {
    "HOH", "WAT", "DOD", "SO4", "PO4", "CL", "NA", "K", "MG", "MN", "CA", "ZN", "FE", "NI", "CO",
    "GOL", "EDO", "PEG", "MPD", "ACT", "ACE", "TRS", "MES", "HEP", "BME", "DMS", "DMSO", "FMT",
}


def search_rcsb(accession: str, cache_dir: Path) -> list[str]:
    query = {
        "query": {
            "type": "terminal",
            "service": "text",
            "parameters": {
                "attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession",
                "operator": "exact_match",
                "value": accession,
            },
        },
        "return_type": "entry",
        "request_options": {"paginate": {"start": 0, "rows": 50}, "sort": [{"sort_by": "score", "direction": "desc"}]},
    }
    data = http_post_json("https://search.rcsb.org/rcsbsearch/v2/query", query, cache_dir / "rcsb" / f"{accession}_search.json")
    return [x["identifier"] for x in data.get("result_set", []) if x.get("identifier")]


def download_pdb(pdb_id: str, cache_dir: Path) -> Path:
    path = cache_dir / "rcsb" / f"{pdb_id.lower()}.pdb"
    if path.exists() and path.stat().st_size > 0:
        return path
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    http_get(url, path, "chemical/x-pdb")
    return path


def download_alphafold(accession: str, cache_dir: Path) -> Path:
    api_url = f"https://alphafold.ebi.ac.uk/api/prediction/{accession}"
    try:
        preds = http_json(api_url, cache_dir / "alphafold" / f"{accession}_meta.json")
        if preds and isinstance(preds, list) and preds[0].get("pdbUrl"):
            pdb_url = preds[0]["pdbUrl"]
            filename = Path(urllib.parse.urlparse(pdb_url).path).name
            path = cache_dir / "alphafold" / filename
            http_get(pdb_url, path, "chemical/x-pdb")
            if path.exists() and path.stat().st_size > 0:
                return path
    except Exception:
        pass

    for version in (6, 5, 4, 3, 2, 1):
        path = cache_dir / "alphafold" / f"AF-{accession}-F1-model_v{version}.pdb"
        if path.exists() and path.stat().st_size > 0:
            return path
        url = f"https://alphafold.ebi.ac.uk/files/AF-{accession}-F1-model_v{version}.pdb"
        try:
            http_get(url, path, "chemical/x-pdb")
            if path.exists() and path.stat().st_size > 0:
                return path
        except Exception:
            path.unlink(missing_ok=True)
    raise StructureNotFound(f"未能为 Accession {accession} 下载到 AlphaFold 预测结构。")


def parse_pdb_atoms(path: Path) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    atoms: list[dict[str, Any]] = []
    hets: list[dict[str, Any]] = []
    for line in path.read_text(errors="ignore").splitlines():
        rec = line[:6].strip()
        if rec not in {"ATOM", "HETATM"}:
            continue
        try:
            atom = {
                "record": rec,
                "name": line[12:16].strip(),
                "altloc": line[16].strip(),
                "resn": line[17:20].strip(),
                "chain": line[21].strip() or "_",
                "resi": int(line[22:26].strip()),
                "x": float(line[30:38].strip()),
                "y": float(line[38:46].strip()),
                "z": float(line[46:54].strip()),
                "element": line[76:78].strip() if len(line) >= 78 else "",
            }
        except Exception:
            continue
        if rec == "ATOM":
            atoms.append(atom)
        else:
            hets.append(atom)
    return atoms, hets


def ligand_groups(hets: list[dict[str, Any]]) -> list[list[dict[str, Any]]]:
    groups: dict[tuple[str, str, int], list[dict[str, Any]]] = {}
    for h in hets:
        resn = h["resn"]
        if resn in COMMON_NON_LIGANDS:
            continue
        key = (h["chain"], resn, h["resi"])
        groups.setdefault(key, []).append(h)
    return [g for g in groups.values() if len(g) >= 6]


def distance(a: dict[str, Any], b: dict[str, Any]) -> float:
    return math.sqrt((a["x"] - b["x"]) ** 2 + (a["y"] - b["y"]) ** 2 + (a["z"] - b["z"]) ** 2)


def box_from_atoms(atoms: list[dict[str, Any]], pad: float = 10.0, min_size: float = 18.0, max_size: float = 30.0) -> tuple[tuple[float, float, float], tuple[float, float, float]]:
    xs = [a["x"] for a in atoms]
    ys = [a["y"] for a in atoms]
    zs = [a["z"] for a in atoms]
    center = ((min(xs) + max(xs)) / 2, (min(ys) + max(ys)) / 2, (min(zs) + max(zs)) / 2)
    size = tuple(max(min_size, min(max_size, (max(vals) - min(vals)) + pad)) for vals in (xs, ys, zs))
    return center, size  # type: ignore[return-value]


def choose_cocrystal_pocket(path: Path) -> tuple[tuple[float, float, float], tuple[float, float, float], str] | None:
    atoms, hets = parse_pdb_atoms(path)
    groups = ligand_groups(hets)
    if not groups:
        return None
    scored: list[tuple[int, int, list[dict[str, Any]]]] = []
    for group in groups:
        contacts = 0
        for lig_atom in group:
            if any(distance(lig_atom, atom) <= 4.5 for atom in atoms):
                contacts += 1
        if contacts >= 3:
            scored.append((contacts, len(group), group))
    if not scored:
        return None
    scored.sort(reverse=True, key=lambda x: (x[0], x[1]))
    ligand = scored[0][2]
    center, size = box_from_atoms(ligand, pad=10.0, min_size=18.0, max_size=30.0)
    return center, size, ligand[0]["resn"]


def feature_positions(features: list[dict[str, Any]]) -> list[int]:
    positions: list[int] = []
    for feature in features:
        ftype = feature.get("type", "")
        if ftype not in {"Active site", "Binding site", "Metal binding"}:
            continue
        loc = feature.get("location", {})
        start = loc.get("start", {}).get("value")
        end = loc.get("end", {}).get("value", start)
        if isinstance(start, int) and isinstance(end, int):
            positions.extend(range(start, end + 1))
    return sorted(set(positions))


def choose_functional_site_pocket(gene: GeneInfo, structure_path: Path) -> tuple[tuple[float, float, float], tuple[float, float, float], str] | None:
    positions = feature_positions(gene.features)
    if not positions:
        return None
    atoms, _ = parse_pdb_atoms(structure_path)
    site_atoms = [atom for atom in atoms if atom.get("resi") in positions]
    if not site_atoms:
        return None
    center, size = box_from_atoms(site_atoms, pad=10.0, min_size=18.0, max_size=30.0)
    return center, size, "uniprot_functional_site"


def choose_metal_site_pocket(structure_path: Path) -> tuple[tuple[float, float, float], tuple[float, float, float], str] | None:
    atoms, hets = parse_pdb_atoms(structure_path)
    metals = [h for h in hets if h.get("resn") in {"ZN", "FE", "MG", "MN", "CA"} or h.get("element") in {"ZN", "FE", "MG", "MN", "CA"}]
    if not metals:
        return None
    site_atoms = list(metals)
    for atom in atoms:
        if any(distance(atom, metal) <= 6.0 for metal in metals):
            site_atoms.append(atom)
    center, size = box_from_atoms(site_atoms, pad=12.0, min_size=18.0, max_size=30.0)
    return center, size, f"metal_site_{metals[0]['resn'].lower()}"


def run_p2rank(structure_path: Path, job_dir: Path) -> tuple[tuple[float, float, float], tuple[float, float, float], str] | None:
    cmd = detect_p2rank()
    if not cmd:
        return None
    out_dir = job_dir / "p2rank"
    if out_dir.exists():
        shutil.rmtree(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    try:
        proc = subprocess.run([cmd, "predict", "-f", str(structure_path), "-o", str(out_dir)], capture_output=True, text=True, timeout=900)
        (job_dir / "p2rank_stdout.txt").write_text(proc.stdout, encoding="utf-8", errors="ignore")
        (job_dir / "p2rank_stderr.txt").write_text(proc.stderr, encoding="utf-8", errors="ignore")
    except Exception:
        return None
    if proc.returncode != 0:
        return None
    csvs = sorted(out_dir.glob("*predictions.csv"), key=lambda p: p.stat().st_size, reverse=True)
    for csv_path in csvs:
        reader = csv.DictReader(csv_path.read_text(encoding="utf-8", errors="ignore").splitlines())
        rows = [{(k.strip() if k else k): (v.strip() if isinstance(v, str) else v) for k, v in row.items()} for row in reader]
        for idx, row in enumerate(rows[:3], 1):
            try:
                center = (float(row["center_x"]), float(row["center_y"]), float(row["center_z"]))
                size = (22.0, 22.0, 22.0)
                return center, size, f"p2rank_rank_{idx}"
            except Exception:
                continue
    return None


def run_dock_job(*args, **kwargs):
    from workflow import run_dock_job as execute
    return execute(*args, **kwargs)


def run_batch_csv(*args, **kwargs):
    from workflow import run_batch_csv as execute
    return execute(*args, **kwargs)


def main(argv=None):
    from workflow import main as execute
    return execute(argv)


if __name__ == "__main__":
    raise SystemExit(main())
