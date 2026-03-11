#!/usr/bin/env bash
#
# Helixer-only gene prediction (no RNA-seq, no homology evidence needed)
# Runs Helixer ab initio prediction inside the Sylvan Singularity container.
#
# Usage:
#   ./bin/helixer.sh <genome.fasta[.gz]> [output.gff3] [lineage]
#
# Arguments:
#   genome       : Input genome FASTA file (required, supports .gz)
#   output.gff3  : Output GFF3 path (default: results/AB_INITIO/Helixer/helixer.gff3)
#   lineage      : Helixer model lineage (default: land_plant)
#                  Options: land_plant, vertebrate, fungi
#
# Environment variables:
#   SYLVAN_SINGULARITY_ARGS  Override singularity bind/GPU args
#   SYLVAN_HELIXER_SUBSEQ    Subsequence length (default: auto from lineage)
#   SYLVAN_HELIXER_MODEL_DIR Model directory inside container
#                            (default: /usr/local/src/Helixer/models)
#
# Examples:
#   ./bin/helixer.sh data/genome.fna.gz
#   ./bin/helixer.sh data/genome.fna.gz data/helixer_out.gff3 land_plant
#   SYLVAN_SINGULARITY_ARGS="--nv -B /scratch" ./bin/helixer.sh data/genome.fna.gz
#

set -euo pipefail

# --- Arguments ---
GENOME="${1:?Usage: $0 <genome.fasta[.gz]> [output.gff3] [lineage]}"
OUTPUT="${2:-results/AB_INITIO/Helixer/helixer.gff3}"
LINEAGE="${3:-land_plant}"

# --- Validate input ---
if [ ! -f "$GENOME" ]; then
	echo "ERROR: Genome file not found: $GENOME" >&2
	exit 1
fi

# --- Singularity image ---
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_DIR="$(dirname "$SCRIPT_DIR")"
SIF="${REPO_DIR}/singularity/sylvan.sif"

if [ ! -f "$SIF" ]; then
	echo "ERROR: Singularity image not found: $SIF" >&2
	echo "Build it with: singularity build --fakeroot singularity/Sylvan.sif singularity/Sylvan.def" >&2
	exit 1
fi

# --- Subsequence length (auto from lineage if not set) ---
if [ -n "${SYLVAN_HELIXER_SUBSEQ:-}" ]; then
	SUBSEQ="$SYLVAN_HELIXER_SUBSEQ"
else
	case "$LINEAGE" in
		land_plant)  SUBSEQ=64152  ;;
		vertebrate)  SUBSEQ=213840 ;;
		fungi)       SUBSEQ=21384  ;;
		*)
			echo "ERROR: Unknown lineage '$LINEAGE'. Use: land_plant, vertebrate, fungi" >&2
			exit 1
			;;
	esac
fi

# --- Model path ---
MODEL_DIR="${SYLVAN_HELIXER_MODEL_DIR:-/usr/local/src/Helixer/models}"
MODEL_FILE="${MODEL_DIR}/${LINEAGE}.h5"

# --- Singularity args ---
REAL_PWD="$(pwd -P)"
SINGULARITY_ARGS="${SYLVAN_SINGULARITY_ARGS:---nv -B $(pwd) -B ${REAL_PWD} -B /tmp}"

# --- Create output directory ---
mkdir -p "$(dirname "$OUTPUT")"

# --- Run ---
echo "=== Helixer ab initio gene prediction ==="
echo "Genome:     $GENOME"
echo "Output:     $OUTPUT"
echo "Lineage:    $LINEAGE"
echo "Subseq:     $SUBSEQ"
echo "Model:      $MODEL_FILE"
echo "Container:  $SIF"
echo "==========================================="

singularity exec $SINGULARITY_ARGS "$SIF" bash -c "
eval \"\$(/usr/local/bin/micromamba shell hook -s bash)\"
micromamba activate helixer
# TF XLA JIT needs libdevice.10.bc — point to pip-installed nvidia cuda_nvcc
export XLA_FLAGS=\"--xla_gpu_cuda_data_dir=\$(python -c 'import nvidia.cuda_nvcc; import os; print(os.path.dirname(nvidia.cuda_nvcc.__file__))' 2>/dev/null || echo '')\"
Helixer.py \
  --fasta-path $GENOME \
  --model-filepath $MODEL_FILE \
  --subsequence-length $SUBSEQ \
  --gff-output-path $OUTPUT
"

# --- Verify output ---
if [ ! -s "$OUTPUT" ]; then
	echo "WARNING: Output GFF3 is empty. Helixer may have failed (GPU/CUDA issue?)." >&2
	exit 1
fi

GENE_COUNT=$(grep -c $'\tgene\t' "$OUTPUT" 2>/dev/null || echo 0)
echo ""
echo "=== Done: $GENE_COUNT genes predicted ==="
echo "Output: $OUTPUT"
