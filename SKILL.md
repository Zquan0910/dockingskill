---
name: dockingskill
description: Retrieve human gene-matched protein and named compound structures, prepare them with RDKit and Meeko, run local whole-protein AutoDock Vina docking with exhaustiveness 64, and generate reference-style PyMOL/PLIP figures. Supports explicit local-pocket docking, batches, provenance, pose export, and verified resume.
metadata:
  short-description: Whole-protein Vina docking and PyMOL figures
---

# DockingSkill

Resolve the skill root, then use its PowerShell launcher. All new outputs, caches,
temporary files and logs belong under D:/ChatGPT/geneprotein-dock-runs.
Never write to C:. Existing installations and source inputs may be read there.

This package lives at D:/ChatGPT/dockingskill and can be explicitly requested as
`$dockingskill` when loaded. A folder outside Codex's configured discovery roots
is not automatically registered: provide this SKILL.md path when requesting it.
Do not copy or install this package into a C: skill directory.
The existing runtime is D:/GeneProteinDocking/env/python.exe; Vina is configured
in config.yaml. See [environment.yml](environment.yml) for the verified dependency
snapshot and [troubleshooting](references/troubleshooting.md) for environment use.
Scripts resolve imports, config and resources from this package, not the old skill.

In commands below, `<skill-root>` means `D:\ChatGPT\dockingskill`.

## Default Workflow

Before any structure download, preparation or docking, check the entire local
toolchain with `doctor --local-only`. Continue only when the required tools and
libraries pass. Missing components must be listed and resolved before proceeding;
do not silently install packages or start a partial workflow. The PowerShell
launcher enforces this for dock/batch, including aliases. Direct Python/API users
must explicitly perform the same preflight first. P2Rank/Java are optional for
the default whole-protein workflow; local-pocket mode may require them.
The full `doctor` additionally checks database connectivity; network warnings
are separate from local dependency failures and do not prevent using local inputs.

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File "<skill-root>\scripts\run.ps1" doctor --local-only
powershell -NoProfile -ExecutionPolicy Bypass -File "<skill-root>\scripts\run.ps1" doctor
powershell -NoProfile -ExecutionPolicy Bypass -File "<skill-root>\scripts\run.ps1" dock --gene ZBTB10 --compound "槲皮素"
```

Default: human target, whole prepared receptor, 5 A padding per side,
exhaustiveness 64, 10 modes, energy range 5 kcal/mol, fixed recorded seed.
Rendering follows a successful docking automatically. Preserve the full box and
requested search strength on failure; report logs instead of silently reducing either.

Experimental target chains require >=95% sequence coverage and >=95% identity.
Otherwise use a full-coverage AlphaFold model. Record source, chain, residue mapping,
missing positions and prediction confidence. Structure selection never requires a pocket.

Use --protein to override automatic structure selection. With a gene, validate its
sequence against the selected local chain. Without a gene, a multi-chain local input
requires --chain. A user-selected fragment stays a fragment and is identified as such.

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File "<skill-root>\scripts\run.ps1" dock --gene ZBTB10 --protein "D:\input\target.pdb" --chain A --ligand "CC(=O)Oc1ccccc1C(=O)O"
powershell -NoProfile -ExecutionPolicy Bypass -File "<skill-root>\scripts\run.ps1" batch "D:\input\jobs.csv"
```

Use [assets/batch_template.csv](assets/batch_template.csv) as the CSV schema.
Its filled rows are real jobs, not an instruction to run a whole batch by default.

Use exact PubChem CID when names are ambiguous. Explicit --ligand overrides a named
compound structure; do not claim PubChem identity for an independently supplied ligand.
Batch CID takes precedence over the compound-name column.

## Parameters and Recovery

CLI > skill-root/config.yaml > built-in defaults. Supported controls:
--exhaustiveness, --num-modes, --energy-range, --seed, --cpu, --timeout,
--padding, --mode whole|local, --center X Y Z, --size SX SY SZ.
In whole mode, manual boxes must enclose every final receptor atom.
Use --mode local explicitly for legacy pocket/P2Rank/manual site docking.

--prepare-only validates structures and box without docking.
--no-render skips images; --no-resume requests a fresh computation.
--out-dir must resolve to D: and sets the parent for unique job directories.
Never reuse solely by job ID: compare prepared structures, actual box, parameters,
tool/code fingerprints and output checksum. Preserve failed attempts and their logs.

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File "<skill-root>\scripts\run.ps1" analyze
powershell -NoProfile -ExecutionPolicy Bypass -File "<skill-root>\scripts\run.ps1" analyze --export "D:\ChatGPT\geneprotein-dock-runs\ranked.csv"
powershell -NoProfile -ExecutionPolicy Bypass -File "<skill-root>\scripts\run.ps1" export-pose "D:\ChatGPT\geneprotein-dock-runs\work\<unique-job>"
powershell -NoProfile -ExecutionPolicy Bypass -File "<skill-root>\scripts\run.ps1" render "D:\ChatGPT\geneprotein-dock-runs\work\<unique-job>"
powershell -NoProfile -ExecutionPolicy Bypass -File "<skill-root>\scripts\run.ps1" render "D:\ChatGPT\geneprotein-dock-runs\batches\<batch>\results.csv"
```

## Figures and Interpretation

Before importing structures, verify that the rendering interpreter can load local
PyMOL and record its version/module path in figures/pymol_environment.json.
If unavailable, report the missing component; do not substitute synthetic graphics.
Import, label and ray-trace in PyMOL. Detail annotations contain only verified
hydrogen-bond lengths and names/numbers of protein residues shown as sticks.
Labels have fully transparent backgrounds and remain editable in scene.pse;
never paste label text or opaque text boxes over the exported structure image.
Use smooth opaque sticks with native lighting and keep the cartoon translucent.
Do not add outer borders. Use equal-width overview/detail layout regions, preserve
each image's native aspect ratio, and vertically centre the shorter detail image.
Mark the actual ligand position in the overview with a yellow dashed rectangle;
connect its upper/lower right corners to the detail image with yellow dashed lines.

Follow assets/reference_style.png: white background; overview on the left,
cartoon/sticks detail on the right; pale cyan protein, green ligand carbon,
red interacting residue carbon, element-coloured heteroatoms. Yellow dashed
connectors identify the zoom region. Yellow dashes show only PLIP-identified
hydrogen bonds with coordinate-verified donor-acceptor distances in angstroms.

Export per-job PNGs, a PyMOL PSE and replay script, all PLIP contacts as JSON/CSV/XML,
and ordered summary PNG/PDF pages (at most four rows/page). The displayed contact
subset is limited for legibility; keep all contacts in the data. Never invent
hydrogen bonds or hide successful jobs based on a score threshold.
Verify images visually as well as checking nonblank pixels and label layout.

Call the value a **Vina docking score (kcal/mol)**. It is not an experimental
binding free energy or proof of binding; do not use fixed -7/-9 cutoffs.
Group score comparisons by target. Interpret AlphaFold confidence and full-box
search uncertainty alongside poses. Do not substitute old scores for new runs.

Read [workflow guide](references/workflow_guide.md) for preparation and artifacts,
[grid strategy](references/pocket_selection.md) for geometry and local mode,
and [troubleshooting](references/troubleshooting.md) when execution fails.
