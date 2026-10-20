"""Behaviour tests use only isolated D: directories; no network or docking CPU required."""
import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import workflow as w
from preparation import whole_box, save_json


def atom_line(x, y, z):
    return f'ATOM      1  CA  ALA A   1    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C\n'


class WorkflowTests(unittest.TestCase):
    def setUp(self):
        (w.ROOT / 'tests').mkdir(exist_ok=True)
        self.temp = tempfile.TemporaryDirectory(dir=w.ROOT / 'tests')
        self.root = Path(self.temp.name)
        self.rec = self.root / 'receptor.pdbqt'
        self.rec.write_text(atom_line(-100, -40, 0) + atom_line(110, 70, 90))

    def tearDown(self):
        self.temp.cleanup()

    def test_box_has_no_30_angstrom_cap(self):
        center, size, bounds = whole_box(self.rec)
        self.assertEqual(center, [5, 15, 45])
        self.assertEqual(size, [220, 120, 100])
        self.assertTrue(bounds['all_atoms_inside'])

    def test_manual_box_must_enclose_receptor(self):
        with self.assertRaisesRegex(ValueError, 'enclose'):
            whole_box(self.rec, center=[0,0,0], size=[20,20,20])

    def test_invalid_grid(self):
        for kwargs in ({'padding': -1}, {'center': [0,0,0]}, {'size': [1,2,3]},
                       {'center': [float('nan'),0,0], 'size': [200,200,200]}):
            with self.subTest(kwargs=kwargs), self.assertRaises(ValueError):
                whole_box(self.rec, **kwargs)

    def test_configuration_precedence(self):
        with patch.object(w.sources, 'read_config', return_value={}):
            self.assertEqual(w.settings({})['exhaustiveness'], 64)
        with patch.object(w.sources, 'read_config', return_value={'docking': {'exhaustiveness': 80}}):
            self.assertEqual(w.settings({})['exhaustiveness'], 80)
            self.assertEqual(w.settings({'exhaustiveness': 96})['exhaustiveness'], 96)

    def test_outputs_reject_c_drive(self):
        with self.assertRaises(ValueError):
            w.d_path('C:/should-never-be-written')

    def test_local_doctor_reports_missing_library_without_network(self):
        import doctor
        def module(name, hint):
            return {'name': name, 'ok': name != 'pymol'}
        with patch.object(doctor, 'check_module', side_effect=module), \
             patch.object(doctor, 'check_python', return_value={'ok': True}), \
             patch.object(doctor, 'check_vina', return_value={'name': 'vina', 'ok': True}), \
             patch.object(doctor, 'check_p2rank', return_value={'name': 'p2rank', 'ok': False}), \
             patch.object(doctor, 'check_network') as network:
            report = doctor.run_diagnostics(local_only=True)
        self.assertFalse(report['healthy'])
        self.assertTrue(report['all_tools_ok'])
        self.assertTrue(report['network_skipped'])
        network.assert_not_called()

    def test_doctor_runs_without_scientific_packages(self):
        import subprocess
        script = Path(__file__).resolve().parents[1] / 'doctor.py'
        proc = subprocess.run([sys.executable, '-S', '-B', str(script), '--local-only'],
                              capture_output=True, text=True, timeout=60)
        self.assertEqual(proc.returncode, 1)
        self.assertIn('[FAIL] rdkit', proc.stdout)
        self.assertIn('[FAIL] pymol', proc.stdout)
        self.assertIn('DEGRADED', proc.stdout)
        self.assertNotIn('Traceback', proc.stderr)

    def test_launcher_blocks_docking_when_preflight_fails(self):
        import os
        import subprocess
        wrapper = self.root / 'python-without-packages.cmd'
        wrapper.write_text('@echo TOOL_CALL:%*\n@"' + sys.executable + '" -S %*\n')
        env = dict(os.environ, DOCKINGSKILL_PYTHON=str(wrapper))
        script = Path(__file__).resolve().parents[1] / 'run.ps1'
        proc = subprocess.run(['powershell', '-NoProfile', '-ExecutionPolicy', 'Bypass', '-File',
                               str(script), 'dock', '--gene', 'ZBTB10', '--compound', 'quercetin'],
                              env=env, capture_output=True, text=True, timeout=90)
        self.assertEqual(proc.returncode, 1)
        self.assertIn('doctor.py', proc.stdout)
        self.assertIn('Preflight failed', proc.stdout)
        self.assertNotIn('dock_pipeline.py', proc.stdout)

    def test_resume_checks_inputs_parameters_and_pose_integrity(self):
        vina = self.root / 'vina.exe'
        vina.write_bytes(b'test executable')
        def receptor(model, chain, job, meta):
            out = job / 'receptor.pdbqt'
            out.write_text(self.rec.read_text())
            return out
        def ligand(compound, inp, job, seed):
            out = job / 'ligand.pdbqt'
            out.write_text(str(inp))
            return out
        def dock(vina, rec, lig, center, size, job, params):
            (job / 'vina_out.pdbqt').write_text('REMARK VINA RESULT: -5.0 0 0\n')
            return [-5.0]
        meta = {'source': 'local_file', 'pdb_id': 'TEST'}
        options = dict(job_id='same', protein_path='unused', ligand_input='CCO', work_dir=self.root / 'work',
                       cache_dir=self.root / 'cache', render=False)
        with patch.object(w, 'select_structure', return_value=(None, 'A', meta)), \
             patch.object(w, 'prepare_receptor', side_effect=receptor), \
             patch.object(w, 'prepare_ligand', side_effect=ligand), \
             patch.object(w.sources, 'detect_vina', return_value=vina), \
             patch.object(w, 'run_vina', side_effect=dock) as execute:
            first = w.run_dock_job(**options)
            self.assertEqual(first['status'], 'SUCCESS')
            again = w.run_dock_job(**options)
            self.assertTrue(again['reused'])
            self.assertNotEqual(again['work_dir'], first['work_dir'])
            self.assertEqual(execute.call_count, 1)
            changed = w.run_dock_job(**options, exhaustiveness=65)
            self.assertEqual(changed['status'], 'SUCCESS')
            self.assertEqual(execute.call_count, 2)
            (Path(first['work_dir']) / 'vina_out.pdbqt').write_text('corrupt')
            (Path(again['work_dir']) / 'vina_out.pdbqt').write_text('corrupt')
            w.run_dock_job(**options)
            self.assertEqual(execute.call_count, 3)
            altered = dict(options, ligand_input='CCCC')
            w.run_dock_job(**altered)
            self.assertEqual(execute.call_count, 4)
            self.assertEqual(len(list((self.root / 'work').glob('*/result.json'))), 5)

    def test_failed_run_does_not_use_old_pose(self):
        with patch.object(w, 'select_structure', side_effect=ValueError('bad structure')):
            result = w.run_dock_job(job_id='same', protein_path='bad', ligand_input='CCO', work_dir=self.root / 'work')
        self.assertEqual(result['status'], 'FAILED')
        self.assertNotIn('best_binding_energy_kcal_mol', result)
        self.assertFalse((Path(result['work_dir']) / 'vina_out.pdbqt').exists())

    def test_analyze_empty_export_and_target_groups(self):
        from analyze_results import analyze
        self.assertEqual(analyze([], self.root / 'empty.csv'), [])
        rows = [{'status': 'SUCCESS', 'job_id': str(i), 'gene_name': g, 'best_binding_energy_kcal_mol': e}
                for i, (g,e) in enumerate([('A', -5), ('B', -10), ('A', -8)])]
        ranked = analyze(rows)
        self.assertEqual([r['gene_name'] for r in ranked], ['A','A','B'])
        self.assertEqual(float(ranked[0]['best_binding_energy_kcal_mol']), -8)

    def test_insertion_codes_survive_sequence_mapping(self):
        from Bio.PDB.Chain import Chain
        from Bio.PDB.Residue import Residue
        from preparation import chain_alignment
        chain = Chain('A')
        chain.add(Residue((' ', 5, ' '), 'ALA', ' '))
        chain.add(Residue((' ', 5, 'A'), 'CYS', ' '))
        mapping = chain_alignment(chain, 'AC')
        self.assertEqual(mapping['coverage'], 1)
        self.assertEqual(mapping['mapping'][1], {'uniprot': 2, 'pdb_number': 5, 'insertion': 'A'})

    def test_vina_timeout_preserves_logs(self):
        import subprocess
        from unittest.mock import Mock
        proc = Mock()
        proc.wait.side_effect = [subprocess.TimeoutExpired('vina', 1), 0]
        params = w.settings({'timeout': 1})
        with patch.object(w.subprocess, 'Popen', return_value=proc), self.assertRaisesRegex(RuntimeError, 'timed out'):
            w.run_vina(Path('fake.exe'), self.rec, self.rec, [0,0,0], [20,20,20], self.root, params)
        proc.kill.assert_called_once()
        self.assertTrue((self.root / 'vina_stdout.txt').is_file())
        args = json.loads((self.root / 'vina_run.json').read_text())['command']
        self.assertEqual(args[args.index('--exhaustiveness')+1], '64')

    def test_summary_paginates_after_four_rows(self):
        from PIL import Image
        from render_results import assemble_summary
        panel = self.root / 'panel.png'
        Image.new('RGB', (2800, 1080), 'white').save(panel)
        outputs = assemble_summary([panel] * 5, self.root / 'summary', 4)
        self.assertEqual(len(outputs), 3)
        with Image.open(outputs[0]) as first, Image.open(outputs[1]) as second:
            self.assertEqual(first.size, (2800, 4320))
            self.assertEqual(second.size, (2800, 1080))

    def test_render_checks_local_pymol_before_export(self):
        import render_results as render
        save_json(self.root / 'result.json', {'status': 'SUCCESS', 'best_model': 1})
        with patch.object(render, 'pymol_preflight', side_effect=RuntimeError('missing PyMOL')), \
             patch.object(render, 'export_pose') as export:
            with self.assertRaisesRegex(RuntimeError, 'missing PyMOL'):
                render.render_one(self.root)
            export.assert_not_called()

    def test_native_labels_are_created_in_pymol(self):
        from unittest.mock import Mock
        from render_results import native_labels
        cmd = Mock()
        cmd.get.return_value = '20'
        cmd.get_model.return_value.atom = []
        cmd.get_model.return_value.bond = []
        entries = [{'xyz': [0,0,0], 'text': 'ASP-461'}, {'xyz': [0,0,0], 'text': '2.80 A'}]
        view = [1,0,0,0,1,0,0,0,1,0,0,-50,0,0,0,0,0,0]
        boxes = native_labels(cmd, entries, view)
        self.assertEqual(cmd.label.call_count, 2)
        cmd.set.assert_any_call('label_bg_transparency', 1)
        self.assertEqual(len(boxes), 2)
        a, b = boxes
        self.assertTrue(a[2] < b[0] or a[0] > b[2] or a[3] < b[1] or a[1] > b[3])

    def test_compose_borderless_yellow_marker_and_preserved_proportions(self):
        from PIL import Image, ImageDraw
        from render_results import compose
        folder = self.root / 'figures'
        folder.mkdir()
        overview = Image.new('RGB', (1050,900), 'white')
        ImageDraw.Draw(overview).rectangle((500,400,510,410), fill=(0,200,0))
        overview.save(folder / 'overview.png')
        Image.new('RGB', (1600,900), 'white').save(folder / 'detail.png')
        save_json(folder / 'projection.json', {'labels': [], 'distances': [],
                  'overview_ligand_pixels': [[500,400],[510,410]]})
        compose(self.root, {'gene_name': 'TEST', 'compound_name': 'ligand',
                            'best_binding_energy_kcal_mol': -5})
        layout = json.loads((folder / 'layout_qa.json').read_text())
        self.assertEqual(layout['black_frames'], [])
        self.assertEqual(layout['connector_color'], '#e6b800')
        self.assertEqual(layout['marker_style'], 'dashed')
        left, right = layout['layout_regions']
        self.assertEqual(left[2]-left[0], right[2]-right[0])
        (origin, size, _), (_, detail_size, _) = layout['image_placements']
        self.assertAlmostEqual(size[0]/size[1], 1050/900, places=3)
        self.assertAlmostEqual(detail_size[0]/detail_size[1], 1600/900, places=3)
        box = layout['overview_marker']
        self.assertLess(box[0], origin[0]+500)
        self.assertGreater(box[2], origin[0]+510)
        self.assertEqual((folder / 'detail.png').read_bytes(), (folder / 'detail_annotated.png').read_bytes())


if __name__ == '__main__':
    unittest.main()
