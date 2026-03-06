#!/usr/bin/env bash
#
# Combined local pipeline runner — runs annotate + filter + benchmark sequentially
#
# Usage:
#   ./bin/run_local.sh                    # full pipeline
#   ./bin/run_local.sh -n                 # dry-run both phases
#   ./bin/run_local.sh --forceall         # force rerun everything
#
# Environment variables:
#   SYLVAN_CONFIG           Override annotate config (default: toydata/config/config_annotate_local.yml)
#   SYLVAN_FILTER_CONFIG    Override filter config (default: toydata/config/config_filter_local.yml)
#   SYLVAN_SKIP_ANNOTATE    Set to 1 to skip Phase 1 (annotate)
#   SYLVAN_SKIP_FILTER      Set to 1 to skip Phase 2 (filter)
#   SYLVAN_SKIP_BENCHMARK   Set to 1 to skip Phase 3 (benchmark)

set -e

echo "============================================="
echo " Sylvan Local Pipeline"
echo " $(date)"
echo "============================================="

EXTRA_ARGS="$@"

# Phase 1: Annotate
if [ "${SYLVAN_SKIP_ANNOTATE:-0}" != "1" ]; then
	echo ""
	echo "=== Phase 1: Annotate ==="
	./bin/annotate_local.sh $EXTRA_ARGS
	echo "=== Phase 1 complete ==="
else
	echo "=== Phase 1: SKIPPED (SYLVAN_SKIP_ANNOTATE=1) ==="
fi

# Phase 2: Filter
if [ "${SYLVAN_SKIP_FILTER:-0}" != "1" ]; then
	echo ""
	echo "=== Phase 2: Filter ==="
	./bin/filter_local.sh $EXTRA_ARGS
	echo "=== Phase 2 complete ==="
else
	echo "=== Phase 2: SKIPPED (SYLVAN_SKIP_FILTER=1) ==="
fi

# Phase 3: Benchmark
if [ "${SYLVAN_SKIP_BENCHMARK:-0}" != "1" ]; then
	echo ""
	echo "=== Phase 3: Benchmark ==="
	./bin/benchmark_local.sh $EXTRA_ARGS
	echo "=== Phase 3 complete ==="
else
	echo "=== Phase 3: SKIPPED (SYLVAN_SKIP_BENCHMARK=1) ==="
fi

echo ""
echo "============================================="
echo " Sylvan Pipeline Complete"
echo " $(date)"
echo "============================================="
echo ""
echo "Key outputs:"
echo "  Annotation: results/PREFILTER/Sylvan.gff3"
echo "  Filtered:   results/FILTER/filtered.gff3"
echo "  Benchmark:  results/BENCHMARK/benchmark_summary.tsv"
