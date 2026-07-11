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


if __name__ == '__main__':
    unittest.main()
