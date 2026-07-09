#!/usr/bin/env bash
# rerun_genewise.sh <run-dir> [--apply]
#
# Re-run the GeneWise homolog branch of an existing Sylvan annotate run so it
# picks up the corrected miniprot region selection.
#
# WHY THIS EXISTS
#   Before commit 5bc3325, bin/miniprot2Genewise.py keyed regions by the
#   `Target=` attribute and let the LAST mRNA line win. miniprot emits one mRNA
#   per alignment ranked best-first (Rank=1 is the best), so every protein with
#   more than one hit was handed to GeneWise at its WORST locus -- sometimes on
#   a different chromosome. Measured on a finished run: 9% of proteins moved.
#   Runs started before that commit therefore carry degraded protein evidence.
#
# WHAT IT DOES
#   1. Refuses to touch a run whose bin/ predates the fix.
#   2. Backs up every artifact the re-run will overwrite (per-group GeneWise
#      GFFs, the merged genewise.gff/.gff3, and the EVM-converted evidence).
#   3. Touches the `miniprot2genewise` checkpoint's input so Snakemake re-runs
#      the checkpoint, then GeneWise, then the merge, then the EVM conversion.
#      Nothing is deleted: the controller's `--rerun-triggers mtime` does the
#      rest on its next invocation.
#
# COST: one GeneWise job per ~100 proteins (1,000-1,500 per run, 2 CPUs each);
# 1-4 h wall per run. Only worth doing BEFORE the EVM consensus step has run --
# check that results/EVM/ holds no consensus output first.
#
# Dry-run by default. Pass --apply to make changes.
set -euo pipefail

RUN="${1:?usage: rerun_genewise.sh <run-dir> [--apply]}"
APPLY="${2:-}"
[ -d "$RUN" ] || { echo "not a directory: $RUN" >&2; exit 1; }
RUN="$(cd "$RUN" && pwd)"

say() { printf '%s\n' "$*"; }
die() { printf 'ERROR: %s\n' "$*" >&2; exit 1; }

# ---- 1. refuse to run against pre-fix code -------------------------------
SCRIPT="$RUN/bin/miniprot2Genewise.py"
[ -f "$SCRIPT" ] || die "missing $SCRIPT"
grep -q 'Rank=' "$SCRIPT" || die \
"$SCRIPT predates the best-hit fix (commit 5bc3325).
       Copy bin/miniprot2Genewise.py and bin/fasta_utils.py from a checkout
       that contains it, or re-seed the run, then try again."

# ---- 2. locate the checkpoint's input ------------------------------------
shopt -s nullglob
CLEAN=( "$RUN"/results/PROTEIN/*.miniprot.clean.gff3 )
shopt -u nullglob
[ "${#CLEAN[@]}" -eq 1 ] || die "expected exactly one results/PROTEIN/*.miniprot.clean.gff3, found ${#CLEAN[@]}"
CLEAN="${CLEAN[0]}"

TMPDIR_GW="$RUN/results/GETA/homolog/geneRegion_genewise.tmp"
[ -d "$TMPDIR_GW" ] || die "no $TMPDIR_GW -- this run has not reached miniprot2genewise yet, nothing to redo"

N_GROUPS=$(find "$TMPDIR_GW" -maxdepth 1 -name '*.faa' | wc -l)
DONE=$(find "$TMPDIR_GW" -maxdepth 1 -name '*.gff' | wc -l)

# ---- 3. warn if the EVM consensus already consumed this evidence ----------
CONSENSUS=$(find "$RUN/results/EVM" -maxdepth 1 \( -name 'evm.out*' -o -name '*.EVM.gff3' \) 2>/dev/null | wc -l)

say "run              : $RUN"
say "checkpoint input : ${CLEAN#$RUN/}"
say "groups           : $N_GROUPS  (GeneWise jobs the re-run will submit)"
say "per-group GFFs   : $DONE  (will be regenerated)"
say "EVM consensus    : $([ "$CONSENSUS" -eq 0 ] && echo 'not run -- safe' || echo "$CONSENSUS file(s) -- downstream WILL be invalidated")"
say ""

if [ "$APPLY" != "--apply" ]; then
	say "DRY RUN. Nothing changed. Re-run with --apply to:"
	say "  1. back up per-group *.gff, genewise.gff, genewise.gff3, EVM/genewise* "
	say "  2. touch $(basename "$CLEAN") so the checkpoint re-runs"
	exit 0
fi

# ---- 4. back up everything the re-run overwrites --------------------------
STAMP=$(date +%Y%m%d-%H%M%S)
BK="$RUN/results/GETA/homolog/genewise_backup_$STAMP"
mkdir -p "$BK"
say "backing up to $BK"
tar -cf "$BK/group_gff.tar" -C "$TMPDIR_GW" --files-from <(cd "$TMPDIR_GW" && find . -maxdepth 1 -name '*.gff' -printf '%P\n')
BK_N=$(tar -tf "$BK/group_gff.tar" | wc -l)
[ "$BK_N" -eq "$DONE" ] || die "backup incomplete: tar has $BK_N entries, expected $DONE"
for f in "$RUN/results/GETA/homolog/genewise.gff" "$RUN/results/GETA/homolog/genewise.gff3"; do
	[ -f "$f" ] && cp -p "$f" "$BK/"
done
for f in "$RUN"/results/EVM/genewise*; do
	[ -f "$f" ] && cp -p "$f" "$BK/"
done
cp -p "$SCRIPT" "$BK/miniprot2Genewise.py"
say "backed up $BK_N per-group GFFs + merged products"

# ---- 5. trigger ----------------------------------------------------------
touch "$CLEAN"
say ""
say "touched $CLEAN"
say "The next Snakemake invocation will re-run: miniprot2genewise (checkpoint)"
say "  -> geneRegion2Genewise x $N_GROUPS -> merge_geneRegion2Genewise -> genewiseGFF2GFF3"
say ""
say "To roll back: restore $BK/* and \`touch\` the restored files so they are"
say "newer than $(basename "$CLEAN")."
