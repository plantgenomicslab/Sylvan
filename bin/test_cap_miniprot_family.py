#!/usr/bin/env python3
"""Opt-in family-diversity tests for cap_miniprot_depth.py."""
import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from cap_miniprot_depth import cap_miniprot_lines, family_key  # noqa: E402
from test_cap_miniprot_depth import _ids, _seg  # noqa: E402


class TestFamilyCap(unittest.TestCase):
    def test_ortho_isoforms_share_family(self):
        self.assertEqual(family_key("AL1006U10010.t1", "MP1"),
                         family_key("AL1006U10010.t2", "MP2"))

    def test_unknown_targets_are_conservative_singletons(self):
        self.assertNotEqual(family_key("opaque", "MP1"),
                            family_key("opaque", "MP2"))

    def test_opt_in_limits_family_but_preserves_diversity(self):
        lines = [
            _seg("MP1", 100, 200, .99, target="AL1006U10010.t1"),
            _seg("MP2", 110, 210, .98, target="AL1006U10010.t2"),
            _seg("MP3", 120, 220, .80, target="AL1006U99999.t1"),
        ]
        kept = cap_miniprot_lines(lines, 10, 10, family_cap=1)
        self.assertEqual(_ids(kept), {"MP1", "MP3"})

    def test_default_disabled_retains_previous_admission(self):
        lines = [
            _seg("MP1", 100, 200, .99, target="AL1006U10010.t1"),
            _seg("MP2", 110, 210, .98, target="AL1006U10010.t2"),
        ]
        self.assertEqual(cap_miniprot_lines(lines, 10, 10), lines)

    def test_default_bytes_match_prechange_fixture_expectation(self):
        """Freeze the existing suite's file-interface fixture byte for byte."""
        lines = [_seg("MP1", 100, 200, .9),
                 _seg("MP2", 150, 250, .8),
                 _seg("MP3", 180, 280, .7),
                 _seg("MP4", 900, 950, .95)]
        expected = "\n".join((lines[0], lines[1], lines[3])) + "\n"
        actual = "\n".join(cap_miniprot_lines(lines, 2, 9)) + "\n"
        self.assertEqual(actual.encode(), expected.encode())

    def test_invalid_cap_rejected(self):
        with self.assertRaises(ValueError):
            cap_miniprot_lines([_seg("MP1", 1, 10, .9)], 2, 2,
                               family_cap=0)


if __name__ == "__main__":
    unittest.main()
