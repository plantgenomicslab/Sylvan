#!/usr/bin/env bash
#
# Setup non-conda tools for the Sylvan pipeline
#
# Usage:
#   ./bin/setup_tools.sh [--geta-home /path/to/geta] [--pfam-dir /path/to/pfamdb]
#
set -e

GETA_HOME="${GETA_HOME:-tools/geta}"
PFAM_DIR="${PFAM_DIR:-databases/Pfam}"

# Parse arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        --geta-home) GETA_HOME="$2"; shift 2;;
        --pfam-dir) PFAM_DIR="$2"; shift 2;;
        *) echo "Unknown option: $1"; exit 1;;
    esac
done

echo "=== Sylvan Non-Conda Tool Setup ==="

# 1. Clone GETA (Perl/Python scripts, no compilation needed)
if [ ! -d "$GETA_HOME" ]; then
    echo "Cloning GETA to $GETA_HOME ..."
    git clone https://github.com/plantgenomicslab/geta "$GETA_HOME"
else
    echo "GETA already exists at $GETA_HOME"
fi

# 2. Download Pfam database (required for GETA's HMM scanning and filter phase)
if [ ! -f "$PFAM_DIR/Pfam-A.hmm" ]; then
    echo "Downloading Pfam-A.hmm to $PFAM_DIR ..."
    mkdir -p "$PFAM_DIR"
    wget -P "$PFAM_DIR" https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
    gunzip "$PFAM_DIR/Pfam-A.hmm.gz"
    echo "Running hmmpress ..."
    hmmpress "$PFAM_DIR/Pfam-A.hmm"
else
    echo "Pfam-A.hmm already exists at $PFAM_DIR"
fi

# 3. Pre-create all conda environments (optional, speeds up first pipeline run)
echo ""
echo "To pre-create conda environments, run:"
echo "  snakemake --use-conda --conda-frontend micromamba --conda-create-envs-only -s bin/Snakefile_annotate"
echo "  snakemake --use-conda --conda-frontend micromamba --conda-create-envs-only -s bin/Snakefile_filter"

echo ""
echo "=== Setup complete ==="
echo "GETA_HOME: $GETA_HOME"
echo "PFAM_DIR:  $PFAM_DIR"
