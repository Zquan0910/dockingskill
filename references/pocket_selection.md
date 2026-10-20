# Search Box Geometry

## Default: Whole Prepared Receptor

Use every ATOM/HETATM coordinate in receptor.pdbqt, including retained metals.
For each axis, centre=(min+max)/2, size=(max-min)+2*padding.
Default padding is 5 A per side, without a maximum edge cap.
Round centres to 0.001 A and round sizes outward so containment is preserved.
Verify all final receptor atoms are inside before invoking Vina.

Manual --center and --size must appear together, be finite, and have positive
dimensions. In whole mode they must enclose all atoms. A large box may require
substantial memory; fail with logs on resource limits instead of shrinking it.
Enclosing the downloaded coordinates does not imply coverage of the full sequence:
check structure_selection.json separately.

## Explicit Local Mode

--mode local enables site-restricted docking. A manual local box may exclude
protein atoms. Automatic local mode can use a bound small molecule contacting
the selected target chain, mapped UniProt functional residues, or P2Rank.
Experimental residue IDs must be mapped rather than equated to UniProt positions.
Local boxes are deliberately distinct from the default whole-protein contract.

The default full-protein workflow does not require P2Rank or Java.
