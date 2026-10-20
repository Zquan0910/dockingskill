# Workflow and Artifacts

## Retrieval and Selection

Query UniProt with an exact human gene symbol (taxonomy 9606); prefer a unique
reviewed entry and reject ambiguity. Download linked experimental PDB files and
align observed protein-chain sequences to the target sequence using Biopython.
Require coverage >=0.95 and identity >=0.95 for automatic experimental selection.
Fallback to AlphaFold via its API; validate sequence coverage and record pLDDT.
Keep the first coordinate model and one target chain; local multi-chain inputs
without a gene must specify the chain. Preserve insertion codes in residue mapping.
Explicit local files override automatic selection but still undergo gene matching.

Retrieve compounds from PubChem by translated common Chinese name, English name
or CID; reject multi-CID names. Prefer 3D SDF and fallback to isomeric SMILES.
Direct ligand files/SMILES override compound lookup.

## Preparation

Protein: use Biopython for model/chain/altloc selection, preferring altloc A then
highest occupancy. Keep target-chain protein and metals within 4 A of it; remove
water, free ligands and unrelated polymers. Fail for unsupported polymer residues.
Meeko prepares the rigid PDBQT without permissive bad-residue deletion.
Compare original heavy atoms and residue identifiers against prepared coordinates;
verify retained metal positions. Fail and preserve logs on losses.

Ligand: preserve original file/identity; choose the largest fragment, preserve
stereochemistry, add H, use existing 3D coordinates or ETKDGv3 (20 conformers).
Optimise every candidate with MMFF94, fallback UFF, and keep the lowest energy.
Record convergence and charge. Meeko assigns docking atom types and torsions.
No pH microstate enumeration or missing-fragment modelling is performed.

## Docking and Resume

Calculate the whole box from the final receptor PDBQT (see grid strategy).
Use default exhaustiveness 64, num_modes 10, energy_range 5, seed 20260715,
cpu 0 (Vina automatic CPU selection), timeout 7200 seconds.
Capture complete stdout/stderr even on timeout. New attempts have unique directories.
SHA256 fingerprints include prepared files, actual grid, all docking parameters,
Vina binary, preparation package versions and script contents. Resume additionally
requires SUCCESS and the original output pose checksum.

Scores come from the actual output PDBQT; select the lowest scored model.
A render failure leaves docking SUCCESS with render_status FAILED and a nonzero
CLI exit, so an expensive successful calculation is still recoverable.

## Outputs

Each work/<unique-job>/ contains result.json, structure_selection.json, the original
receptor, receptor/ and ligand/ metadata/logs, grid.json, fingerprint.json,
vina_run.json, stdout/stderr and vina_out.pdbqt. Errors include error.log.

Pose export reconstructs chemistry using Meeko's SMILES/atom-index mapping and
checks the exported molecular identity against ligand_prepared.sdf. It assigns
a distinct ligand chain and noncolliding PDB serials, preserving receptor numbering.

figures/ holds overview.png, detail.png, detail_annotated.png, panel.png,
scene.pse, reproduce.py, PLIP report.xml and interactions.json/CSV, projection.json,
layout_qa.json, summary PNG/PDF. Hydrogen-bond lines show donor-acceptor distances;
PyMOL availability is checked before structure import and recorded in
pymol_environment.json. Residue names/numbers and hydrogen-bond lengths are native
PyMOL labels with fully transparent backgrounds, saved in scene.pse and ray-traced
into detail.png. detail_annotated.png is an identical compatibility copy, not a
post-processed label overlay. The composed overview and detail use a
borderless layout with equal-width regions and aspect-preserving image sizes.
The overview ligand is enclosed by a yellow dashed rectangle, connected to the
detail image by two yellow dashed lines. Label rectangles in layout_qa.json are conservative placement
footprints, not measured glyph masks; visual inspection is still required.
Other contact types remain available in the data. Up to four shortest hydrogen
bonds and their residues are drawn; when none exist, show detected contact residues
without inventing bonds. Each batch records results.csv in a unique batches/ folder.
