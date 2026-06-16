# Running Helixer on GPU

Helixer (the `helixer` rule in `bin/Snakefile_annotate`, and the standalone
`bin/helixer.sh`) is a deep-learning ab initio gene predictor. It runs on CPU
by default, but **CPU inference on large (multi-Gb) genomes is impractical** —
it is extremely slow and routinely exhausts host RAM (OOM). For anything larger
than a small genome, run Helixer on a GPU.

This page explains how to configure the pipeline config YAML to request a GPU,
and what to do when your SLURM cluster has no GPU partition.

---

## 1. Enable GPU for the `helixer` rule in the config YAML

The config file (`config/config_annotate.yml`, or your project copy such as
`config/config_deschampsia.yml`) doubles as Snakemake's `--cluster-config`.
Each rule block under it maps to `sbatch` options. The `helixer` rule block
exposes three GPU-related keys (commented out by default):

```yaml
helixer:
  account: gpu-account          # GPU-specific SLURM account (if different from __default__)
  partition: gpu-partition      # GPU-specific SLURM partition (e.g. gpu-s1-pgl-0)
  extra_args: "--gres=gpu:1"    # Request 1 GPU from SLURM. Remove/empty for CPU-only.
  ncpus: 12
  threads: 12
  memory: 48g
```

To enable GPU:

1. Set `partition` (and `account`, if your GPU partition uses a different one)
   to your cluster's GPU partition/account.
2. Set `extra_args: "--gres=gpu:1"` to request one GPU. Adjust the GRES string
   to match your scheduler (e.g. `--gres=gpu:a100:1`).
3. Leave `memory` generous (host RAM is still used for genome numerification and
   `HelixerPost`); 48 GB is a reasonable floor, more for very large genomes.

Everything else (Singularity `--nv` GPU passthrough) is already handled by the
entry scripts.

> **CPU-only:** leave `extra_args` empty (or omit the block) and Helixer runs on
> CPU. Only viable for small genomes. For large genomes expect OOM — see §4.

---

## 2. XLA / libdevice requirement (handled automatically)

Helixer's TensorFlow backend JIT-compiles some ops with XLA, which needs
`libdevice.10.bc` and `ptxas` from the CUDA toolkit. These ship in the
pip-installed `nvidia-cuda-nvcc` package inside the `helixer` conda env.

The pipeline locates them automatically and sets:

```bash
CUDA_NVCC_DIR=$(python -c 'import nvidia.cuda_nvcc as m; print(list(m.__path__)[0])')
export XLA_FLAGS="--xla_gpu_cuda_data_dir=$CUDA_NVCC_DIR"
export PATH="$CUDA_NVCC_DIR/bin:$PATH"
```

> `nvidia.cuda_nvcc` is a *namespace* package — its `__file__` is `None`, so
> `__path__[0]` (not `__file__`) must be used to find its directory. Using
> `__file__` yields an empty `XLA_FLAGS` and the run fails with
> `libdevice not found at ./libdevice.10.bc` / `JIT compilation failed`.

If you build your own container, make sure `nvidia-cuda-nvcc` is importable in
the `helixer` env (i.e. `python -c 'import nvidia.cuda_nvcc'` succeeds).

---

## 3. Standalone run on a GPU node

To run Helixer outside Snakemake (e.g. on a dedicated GPU host):

```bash
SYLVAN_SINGULARITY_ARGS="--nv -B /path/to/data -B /tmp" \
  ./bin/helixer.sh <genome.fasta> [output.gff3] [lineage]
# lineage: land_plant (default) | vertebrate | fungi
```

`bin/helixer.sh` applies the same XLA/libdevice detection as the rule.

---

## 4. When the SLURM cluster has no GPU partition

If your annotation cluster is CPU-only, do **not** run the `helixer` rule there
on a large genome — it will OOM and (without the silent fallback, which has been
removed) the rule will now fail loudly instead of producing an empty
`helixer.gff3`.

Recommended pattern:

1. Run Helixer on a separate GPU host (e.g. `./bin/helixer.sh` as in §3), with
   the **same genome FASTA** used by the pipeline (identical sequence IDs).
2. Copy the result into the pipeline tree:
   `results/AB_INITIO/Helixer/helixer.gff3`.
3. Run/resume the pipeline. Because the supplied `helixer.gff3` is newer than
   its input genome, Snakemake (`--rerun-triggers mtime`) will **not** re-run the
   `helixer` rule; it uses the provided file. Downstream rules that consume it
   (`pasaPrep`, `combineEVMInputs`) re-run as needed.

> Verify sequence-ID compatibility: the seqids in `helixer.gff3` must match the
> pipeline genome (`results/GETA/genome.fasta`). A GFF3 produced against a
> different assembly (e.g. NCBI accession names vs. `Chromosome1`/`ptg…`) is not
> usable.

---

## 5. Verifying a successful run

```bash
# Non-empty output with gene features
[ -s results/AB_INITIO/Helixer/helixer.gff3 ] && \
  awk -F'\t' '$3=="gene"' results/AB_INITIO/Helixer/helixer.gff3 | wc -l

# In the run log, look for completion and the absence of XLA errors
grep -iE "libdevice|JIT compilation failed" <run-log>   # should be empty
```
