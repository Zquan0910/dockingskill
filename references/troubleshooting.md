# Troubleshooting

Run scripts/run.ps1 doctor. Default outputs and temp files are on
D:/ChatGPT/geneprotein-dock-runs; software may remain in D:/GeneProteinDocking/env.

Before starting a new job, run doctor --local-only and resolve all required
tool/library failures. dock/batch launchers enforce this check and stop before
retrieval/preparation if it fails. No automatic installation is performed.
The doctor itself uses only the Python standard library before checking imports,
so missing scientific packages are reported rather than imported prematurely.
Set DOCKINGSKILL_PYTHON to an existing python.exe to select another environment.

## Packaged environment

The skill itself is D:/ChatGPT/dockingskill. Scripts and config resolve relative
to that folder; the original skill folder is not required. The launcher prefers
D:/GeneProteinDocking/env/python.exe and redirects temporary/cache outputs to D:.
environment.yml records verified direct dependencies, including conda PyMOL and
OpenBabel; it is not a fully solved lockfile or a portable standalone installation.
Vina and optional P2Rank remain external, with their paths in config.yaml.
Reuse the existing environment. Do not run conda/pip installation during routine
docking. If rebuilding elsewhere is explicitly needed, first set TEMP/TMP/TMPDIR,
CONDA_PKGS_DIRS, CONDA_ENVS_PATH and PIP_CACHE_DIR to D: directories and specify a
D: environment prefix. Verify the solver result and run doctor before use.

This folder is not installed into a Codex discovery directory. To load it from
this location, ask Codex to read D:/ChatGPT/dockingskill/SKILL.md and use the skill.
Do not install or register a copy on C:.

## Vina memory or timeout failure

Old Vina 1.1.2 failed to allocate the full ZBTB10 grid in validation.
config.yaml now points to the official Vina 1.2.7 Windows binary under the new
run root's tools/ directory. Record its checksum in each run.
Inspect vina_stdout.txt, vina_stderr.txt and vina_run.json. Do not shrink the box,
lower exhaustiveness, or substitute an older score. Increase --timeout explicitly
when a longer run is appropriate.

## Incomplete receptor

Meeko runs without --allow_bad_res. Missing side chains or unsupported residues
produce a failure rather than silent deletion. Inspect receptor/meeko_receptor.log
and select a different validated structure. No missing-fragment repair is claimed.
8GN3 failed the residue-completeness check; explicit local 8GN4 prepared successfully
and preserved both target-associated Zn ions during validation.

## Rendering failure

Check result.json render_status, render_error.log and figures/pymol.log.
The D: Python environment needs pymol, PLIP, lxml, OpenBabel, RDKit and Pillow.
Use render <job-directory> to retry after fixing a dependency; do not redo docking.
No detected hydrogen bonds is a valid scientific outcome, not a render failure.
Look at interactions.json for every PLIP contact, not only the displayed subset.

## Name ambiguity, network and cached data

Specify an exact CID for ambiguous compounds. Cached UniProt/PubChem/PDB/AlphaFold
files reside under the new run root's cache/. Network errors are retained in logs.
Local paths can be supplied directly. A gene is needed to verify sequence coverage
for a local structure; otherwise choose an explicit chain in multi-chain structures.

## Resume

Changing a prepared structure, grid, parameter, tool or script invalidates resume.
A corrupted pose checksum also invalidates it. --no-resume always computes anew.
Batch CSV records include work_dir and fingerprint; never copy an old SUCCESS
status into a new task by job ID alone.
