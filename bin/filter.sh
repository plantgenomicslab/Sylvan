#!/usr/bin/env bash

set -x
set -e

snakemake -p --rerun-incomplete --cluster-config config_filter.yml \
		--rerun-triggers mtime \
		--snakefile bin/Snakefile_filter \
		--use-singularity \
		--use-conda \
		--max-jobs-per-second 50 \
		--max-status-checks-per-second 50 \
		--jobs 150 \
		--latency-wait 30 \
		--cluster "sbatch -N {cluster.nodes} --mem={cluster.memory} --cpus-per-task={cluster.ncpus} \
				-J {cluster.name} \
				--parsable -A {cluster.account} -p {cluster.partition} \
				-t {cluster.time} -o {cluster.output} -e {cluster.error}" \
		--singularity-args "--cleanenv --env PYTHONNOUSERSITE=1" \
		"$@"
