#!/usr/bin/env python3
"""Unit test for the filter_feature_importance.py script.

This test simulates the inputs to filter_feature_importance.py to verify that it
runs and produces outputs in the expected format. It does not check the
statistical validity of the results, but rather the script's execution and
output structure.
"""
import unittest
import tempfile
import shutil
import os
import subprocess
import sys
import pandas as pd
import json

# Import the Filter module directly so the RF path can be unit-tested (the
# feature-importance test below only drives it through a subprocess and only
# checks output structure).
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import Filter  # noqa: E402


class TestFeatureImportance(unittest.TestCase):
    def setUp(self):
        """Set up a temporary directory and dummy input files."""
        self.temp_dir = tempfile.mkdtemp()
        self.data_path = os.path.join(self.temp_dir, "data.tsv")
        self.busco_path = os.path.join(self.temp_dir, "full_table.tsv")
        self.output_table_path = os.path.join(self.temp_dir, "feature_importance.tsv")
        self.output_json_path = os.path.join(self.temp_dir, "feature_importance.json")

        # Create dummy data.tsv.
        # semiSupRandomForest splits on the label column: anything != "None" is
        # training data, "None" is the unlabeled set it predicts on. Both halves
        # must be non-empty, and the training half must carry both classes, or
        # RandomForestClassifier has nothing to fit or nothing to predict.
        data = {
            'transcript_id': [f'tx{i}' for i in range(12)],
            'label': (['Keep'] * 4) + (['Discard'] * 4) + (['None'] * 4),
            'feature1': [0.90, 0.85, 0.92, 0.88, 0.10, 0.15, 0.11, 0.13, 0.80, 0.20, 0.75, 0.25],
            'feature2': [0.20, 0.25, 0.22, 0.18, 0.80, 0.75, 0.81, 0.83, 0.30, 0.70, 0.35, 0.65],
            'feature3': [0.70, 0.65, 0.72, 0.68, 0.30, 0.35, 0.31, 0.33, 0.60, 0.40, 0.55, 0.45]
        }
        pd.DataFrame(data).to_csv(self.data_path, sep='\t', index=False)

        # Create dummy full_table.tsv for BUSCO. Sequence must name transcripts
        # that appear in data.tsv, otherwise the kept/discarded BUSCO counts are
        # all zero and the test asserts nothing.
        busco_data = {
            '# Busco id': [f'busco{i}' for i in range(5)],
            'Status': ['Complete'] * 5,
            'Sequence': ['tx0', 'tx1', 'tx2', 'tx4', 'tx5'],
            'Score': [0.9] * 5,
            'Length': [100] * 5
        }
        with open(self.busco_path, "w") as f:
            f.write("# Some header lines\n")
            pd.DataFrame(busco_data).to_csv(f, sep='\t', index=False)

    def tearDown(self):
        """Clean up the temporary directory."""
        shutil.rmtree(self.temp_dir)

    def test_script_runs_and_creates_output(self):
        """Test if the script runs and creates the expected output files."""
        script_path = os.path.join(os.path.dirname(__file__), 'filter_feature_importance.py')
        
        # The script we are testing imports 'Filter' which is in the same directory.
        # We need to make sure python can find it.
        env = os.environ.copy()
        env['PYTHONPATH'] = os.path.dirname(__file__) + os.pathsep + env.get('PYTHONPATH', '')

        cmd = [
            sys.executable, script_path,
            self.data_path,
            self.busco_path,
            '--output-table', self.output_table_path,
            '--output-json', self.output_json_path,
            '--max-iter', '2', # Keep it fast
            # prepare_features emits one <col>_missing flag per feature, so the
            # 3-feature ablation runs see only 4 columns. The default
            # --predictors 6 would exceed that and RandomForestClassifier rejects
            # max_features > n_features.
            '--predictors', '3',
        ]

        result = subprocess.run(cmd, capture_output=True, text=True, env=env, check=False)
        
        self.assertEqual(result.returncode, 0, f"Script failed with exit code {result.returncode}\nstdout:\n{result.stdout}\nstderr:\n{result.stderr}")

        # Check if output files were created
        self.assertTrue(os.path.exists(self.output_table_path), "Output table file was not created.")
        self.assertTrue(os.path.exists(self.output_json_path), "Output json file was not created.")

        # Check the content of the TSV output
        output_df = pd.read_csv(self.output_table_path, sep='\t')
        self.assertEqual(len(output_df), 4) # baseline + 3 features
        expected_columns = [
            'feature_removed', 'num_features', 'oob_delta', 'iterations', 
            'final_kept', 'final_discarded', 'final_kept_buscos', 
            'final_discarded_buscos', 'final_oob_error'
        ]
        self.assertListEqual(list(output_df.columns), expected_columns)
        self.assertEqual(output_df.iloc[0]['feature_removed'], '(none)')

        # Check the content of the JSON output
        with open(self.output_json_path, 'r') as f:
            json_data = json.load(f)
        self.assertIn('runs', json_data)
        self.assertEqual(len(json_data['runs']), 4)
        self.assertEqual(json_data['runs'][0]['feature_removed'], '(none)')


class TestSemiSupRandomForest(unittest.TestCase):
    """Direct tests for Filter.semiSupRandomForest -- the core semi-supervised
    RF path. The existing feature-importance test only reaches it through a
    subprocess and asserts output shape, never that Keep/Discard decisions,
    prediction recycling, or the BUSCO safety net actually behave."""

    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.temp_dir)

    def _write_busco(self, complete_transcripts):
        """Minimal BUSCO full_table (busco_id, status, transcript_id columns)."""
        path = os.path.join(self.temp_dir, "busco.tsv")
        with open(path, "w") as f:
            f.write("# busco_id\tstatus\tsequence\n")
            for i, tx in enumerate(complete_transcripts):
                f.write(f"busco{i}\tComplete\t{tx}\n")
        return path

    @staticmethod
    def _separable_data():
        """A cleanly separable labelled + unlabelled set.

        Keep-like rows carry high feat1/feat2, Discard-like rows carry low
        values. Only numeric features are supplied, so the Pfam/TE rescue tiers
        stay inert and the assertions isolate the RF decision + recycling."""
        rows = []
        for i in range(8):
            rows.append((f"keep{i}", "Keep", 0.85 + 0.01 * (i % 5), 0.80 + 0.01 * (i % 5)))
        for i in range(8):
            rows.append((f"disc{i}", "Discard", 0.10 + 0.01 * (i % 5), 0.12 + 0.01 * (i % 5)))
        for i in range(4):  # unlabelled, looks like Keep
            rows.append((f"ukeep{i}", "None", 0.88 + 0.01 * (i % 3), 0.86 + 0.01 * (i % 3)))
        for i in range(4):  # unlabelled, looks like Discard
            rows.append((f"udisc{i}", "None", 0.09 + 0.01 * (i % 3), 0.11 + 0.01 * (i % 3)))
        return pd.DataFrame(rows, columns=["transcript_id", "label", "feat1", "feat2"])

    def test_rf_keeps_positive_discards_negative(self):
        """Keep seeds + Keep-like unlabelled genes survive; Discard seeds and
        Discard-like unlabelled genes are filtered; diagnostics are well-formed."""
        data = self._separable_data()
        busco = self._write_busco(["keep0", "keep1"])  # BUSCOs only among Keep seeds
        kept_df, process = Filter.semiSupRandomForest(
            data, predictors=2, busco_table=busco, num_trees=100, seed=1, maxiter=5,
        )
        self.assertIsInstance(kept_df, pd.DataFrame)
        self.assertIn("transcript_id", kept_df.columns)
        kept = set(kept_df["transcript_id"])

        for key in ("kept", "discarded", "kept_buscos", "discarded_buscos", "OOB"):
            self.assertIn(key, process)
            self.assertGreater(len(process[key]), 0)
        for oob in process["OOB"]:
            self.assertGreaterEqual(oob, 0.0)
            self.assertLessEqual(oob, 1.0)

        for i in range(8):
            self.assertIn(f"keep{i}", kept)
            self.assertNotIn(f"disc{i}", kept)
        for i in range(4):
            self.assertIn(f"ukeep{i}", kept)
            self.assertNotIn(f"udisc{i}", kept)

    def test_deterministic_with_seed(self):
        """Same seed -> identical kept set (guards accidental nondeterminism)."""
        data = self._separable_data()
        busco = self._write_busco(["keep0"])
        k1, _ = Filter.semiSupRandomForest(data, 2, busco, 100, seed=7, maxiter=5)
        k2, _ = Filter.semiSupRandomForest(data, 2, busco, 100, seed=7, maxiter=5)
        self.assertEqual(set(k1["transcript_id"]), set(k2["transcript_id"]))

    def test_busco_safety_net_rescues_discard_seed(self):
        """A Discard-seeded gene with a Complete BUSCO is force-kept; an
        ordinary Discard gene is not (the full-result BUSCO net, issue #20.2)."""
        rows = [(f"keep{i}", "Keep", 0.9, 0.9) for i in range(6)]
        rows += [(f"disc{i}", "Discard", 0.1, 0.1) for i in range(6)]
        data = pd.DataFrame(rows, columns=["transcript_id", "label", "feat1", "feat2"])
        busco = self._write_busco(["disc0"])  # Discard-seeded gene carrying a BUSCO
        kept_df, _ = Filter.semiSupRandomForest(data, 2, busco, 100, seed=3, maxiter=3)
        kept = set(kept_df["transcript_id"])
        self.assertIn("disc0", kept)
        self.assertNotIn("disc1", kept)

    def test_empty_seed_labels_keeps_everything(self):
        """No heuristic seeds (all 'None') -> keep all transcripts unfiltered
        rather than crash on a 0-row fit (issue #20.4 entry guard)."""
        rows = [(f"g{i}", "None", 0.5, 0.5) for i in range(6)]
        data = pd.DataFrame(rows, columns=["transcript_id", "label", "feat1", "feat2"])
        busco = self._write_busco(["g0"])  # any non-empty table; the guard returns early
        kept_df, _ = Filter.semiSupRandomForest(data, 2, busco, 50, seed=0, maxiter=3)
        self.assertEqual(set(kept_df["transcript_id"]), set(data["transcript_id"]))


if __name__ == '__main__':
    unittest.main()
