"""Regression test for issue #40: GFF3Clear resolves overlapping loci in favour
of whichever input file is listed FIRST, so argument order is load-bearing.

Getting pickBetterModels' order backwards (picked_evidence first) dropped 3,038
of geneModels.a's correct intron chains on Ath_ODB (IC F1 65.25 -> 53.81);
listing a first recovers 3,134 chains (F1 64.57). These tests parse
bin/Snakefile_annotate and fail if either guarded ordering regresses.
"""
import os
import re
import unittest

SNAKEFILE = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                         "Snakefile_annotate")


def _rule_block(text, rule_name):
    """Return the source block of one rule, from its header to the next rule."""
    match = re.search(rf"^rule {rule_name}:$.*?(?=^rule |^def |\Z)",
                      text, re.MULTILINE | re.DOTALL)
    if match is None:
        raise AssertionError(f"rule {rule_name} not found in Snakefile_annotate")
    return match.group(0)


class TestGff3ClearOrder(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        with open(SNAKEFILE, encoding="utf-8") as handle:
            cls.text = handle.read()

    def test_pick_better_models_lists_genemodels_a_before_picked_evidence(self):
        """#40: in the GFF3Clear call that writes geneModels.d, {input.a} must
        come before picked_evidence_geneModels.gff3."""
        block = _rule_block(self.text, "pickBetterModels")
        # Drop comment lines first: the rationale comment above the call also
        # mentions GFF3Clear and picked_evidence_geneModels.gff3.
        code_only = "\n".join(line for line in block.splitlines()
                              if not line.lstrip().startswith("#"))
        call = re.search(r"GFF3Clear.*?>\s*\{output\.d\}", code_only, re.DOTALL)
        self.assertIsNotNone(
            call, "GFF3Clear call producing {output.d} not found in pickBetterModels")
        call_text = call.group(0)
        pos_a = call_text.find("{input.a}")
        pos_picked = call_text.find("picked_evidence_geneModels.gff3")
        self.assertGreater(pos_a, -1, "{input.a} missing from the GFF3Clear call")
        self.assertGreater(pos_picked, -1,
                           "picked_evidence_geneModels.gff3 missing from the GFF3Clear call")
        self.assertLess(
            pos_a, pos_picked,
            "issue #40 regression: picked_evidence_geneModels.gff3 is listed before "
            "{input.a}; GFF3Clear keeps the FIRST file's models on overlap, so this "
            "order costs ~3,100 correct intron chains (Ath_ODB F1 64.57 -> 53.81)")

    def test_transfrag_orf_inputs_list_stranded_first(self):
        """Companion guard: _transfrag_orf_inputs must keep the stranded ORF file
        first so spliced models win over single-exon ones at shared loci."""
        match = re.search(r"^def _transfrag_orf_inputs\(.*?(?=^rule |^def (?!_transfrag))",
                          self.text, re.MULTILINE | re.DOTALL)
        self.assertIsNotNone(match, "_transfrag_orf_inputs not found")
        body = match.group(0)
        pos_strand = body.find("transdecoder2ORF.gff3")
        pos_nostrand = body.find("transdecoder2ORF.nostrand")
        self.assertGreater(pos_strand, -1, "stranded ORF gff3 missing")
        if pos_nostrand != -1:
            self.assertLess(pos_strand, pos_nostrand,
                            "stranded ORF file must be listed before the unstranded one")


if __name__ == "__main__":
    unittest.main(verbosity=2)
