# EVM 개선 4건 파이프라인 반영 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 스크래치에서 검증된 EVM 개선 4건(miniprot depth cap K=16, EVM weight 재배분, Liftoff ABINITIO 승격, mRNA 단위 파티션 병합)과 runEVM 자원 증량을 `bin/Snakefile_annotate`·`config/`에 반영해, 파이프라인이 Ath_ODB 기준 IC F1 ≈73.5(현행 66.4)를 재현하게 만든다.

**Architecture:** DAG를 바꾸지 않는다 — 새 rule 없음, rule의 input/output 불변, 모든 변경은 기존 rule의 shell 블록 + config + `bin/` 헬퍼 스크립트. 파티션 병합은 stock `recombine_EVM_partial_outputs.pl`(행별 midpoint 소속 판정)을 per-partition `EVM_to_GFF3.pl` 변환 + 신규 `bin/merge_evm_partitions.py`(mRNA 모델 단위 core-소유 판정, 검증본 `remerge.py` v2의 프로덕션판)로 대체한다.

**Tech Stack:** Snakemake 7, Python 3.10 (표준 라이브러리만), EVM v2.1.0 (컨테이너 `Sylvan_20260810_spaln.sif` 내 annotation env), GNU sort.

**Spec:** `docs/superpowers/specs/2026-08-25-evm-selftune-design.md` (Phase 1 설계) + `/data/gpfs/assoc/pgl/data/Sylvan/resume.md` 2026-08-26 19:01 스냅샷 (검증된 최종 절차·수치). 검증 원본: 이전 세션 스크래치패드 `a0d6edfc-*/scratchpad/fullbest/{run.sh,merge_score.sh,w.txt}`, `remerge.py`.

## Global Constraints

- **Job stats 동일성**: baseline(HEAD) Snakefile과 새 Snakefile의 dry-run Job stats 블록이 동일해야 한다 (`--rerun-triggers mtime --nolock --scheduler greedy` 필수 — 손dry-run은 기본 트리거로 떨어져 프로덕션과 다른 답을 냄).
- **weights는 항목 삭제가 아니라 0으로 유지**: `evidence_modeler.pl:880`이 weights의 ev_type→ev_class 매핑으로 클래스를 정하므로, 트랙을 끌 때도 행을 지우지 말고 weight 0으로 남긴다 (lookup 실패 시 warn-skip/confess).
- **Snakefile shell 블록에 errexit 없음**: `shell.prefix`는 `-u`만 토글 (bin/Snakefile_annotate:2956 주석). 중간 명령 실패가 마지막 명령 성공에 가려지므로 신규 다단계 셸은 명령마다 `|| exit 1` 명시.
- **cap의 외부 정렬 tmpdir은 GPFS** (`results/TMP`) — compute 노드 `/tmp`는 244G tmpfs(=RAM).
- **검증된 수치와 파라미터를 그대로**: K_base=16, K_intron=16 (스윕 12:−0.91, 20:−0.04, 24:−0.22, 32:−0.42), weight 표는 fullbest `w.txt`와 자구 일치, 병합은 remerge.py v2(complete 우선) 로직 일치.
- 브랜치: `feat/evm-improvements` (현 HEAD `72f612d`에서 분기 — GETA 수정 포함 상태가 fleet 재실행의 실제 배포 대상).
- 커밋 형식: conventional commits (`feat:`/`fix:`/`docs:`/`test:`).

---

### Task 0: 브랜치 준비

**Files:**
- Modify: 없음 (git 조작만)

- [ ] **Step 1: 미커밋 문서 변경 커밋** — `bin/strip_evm_chimeras.py`의 +27줄은 이전 세션이 남긴 docstring 문서화(가드가 주석에서 모델을 지우지 않는다는 실측 기록). 현 브랜치에 docs 커밋으로 정리:

```bash
git add bin/strip_evm_chimeras.py
git commit -m "docs(strip_evm_chimeras): the guard only shields PASA input; chimera removal happens in the RF filter"
```

- [ ] **Step 2: 새 브랜치 생성**

```bash
git checkout -b feat/evm-improvements
```

- [ ] **Step 3: baseline 스냅샷 확보** (Task 6의 dry-run 대조용) — 이 시점의 bin/ 전체를 별도 디렉터리로:

```bash
BASE=/data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-data-Sylvan-Sylvan/d2b5ae8c-cff2-416c-b12e-50393440468f/scratchpad/baseline_bin
mkdir -p $BASE && git archive HEAD bin | tar -x -C $BASE
```

### Task 1: EVM weight 재배분 (+1.3, Liftoff 승격분의 weights 절반)

**Files:**
- Modify: `config/evm_weights.txt` (전체 교체)
- Modify: `toydata/config/evm_weights.txt` (동일 내용으로 교체 — toydata 검증이 프로덕션 로직과 같은 경로를 타도록)

**Interfaces:**
- Produces: weights 파일은 `writeEVMCommands`가 `--weights`로 EVM 명령에 박음. `Liftoff`가 `ABINITIO_PREDICTION` 클래스로 선언되므로 Task 3의 liftoff 파일 이동과 쌍을 이룸.

- [ ] **Step 1: 두 파일을 다음 내용으로 교체** (fullbest `w.txt` 자구 일치, 탭 구분):

```
ABINITIO_PREDICTION	AUGUSTUS	7
ABINITIO_PREDICTION	Helixer	11
ABINITIO_PREDICTION	Liftoff	4
OTHER_PREDICTION	GETA	0
OTHER_PREDICTION	Genewise	2

TRANSCRIPT	gmap.exon_match	0
TRANSCRIPT	assembler-pasa.sqlite	1
TRANSCRIPT	StringTie	0
TRANSCRIPT	PsiClass	5

PROTEIN	GeneWise	0
PROTEIN	miniprot	2
```

변경점: Helixer 3→11, Liftoff OTHER 2→ABINITIO 4, GETA 5→0, gmap 2→0, PASA 10→1, StringTie 1→0, PsiClass 1→5, GeneWise(PROTEIN) 2→0. AUGUSTUS 7, Genewise(OTHER) 2, miniprot 2 유지.

- [ ] **Step 2: 검증** — `diff <(md5sum < config/evm_weights.txt) <(md5sum < toydata/config/evm_weights.txt)` 동일 확인, 탭 구분 확인 (`awk -F'\t' 'NF && NF!=3'` 출력 없음).

- [ ] **Step 3: Commit**

```bash
git add config/evm_weights.txt toydata/config/evm_weights.txt
git commit -m "feat(evm): rebalance EVM weights -- Helixer 11, Liftoff to ABINITIO 4, retire GETA/gmap/StringTie/GeneWise votes"
```

### Task 2: `bin/merge_evm_partitions.py` — mRNA 모델 단위 파티션 병합 (+0.24)

**Files:**
- Create: `bin/merge_evm_partitions.py`
- Test: `bin/test_merge_evm_partitions.py`

**Interfaces:**
- Consumes: `partitions_list.out` (필드: `chrom<TAB>chrom_dir<TAB>Y<TAB>partition_dir` 또는 `chrom<TAB>chrom_dir<TAB>N`), 각 파티션 dir의 `evm.partition.gff3` (파티션 로컬 좌표, Task 4가 `EVM_to_GFF3.pl`로 생성).
- Produces: CLI `python bin/merge_evm_partitions.py <partitions_list.out> <out.gff3>`. 출력: 게놈 좌표 GFF3, ID/Parent에 파티션 태그 접두(`{chrom}_{lo}-{hi}_` 또는 `{chrom}_`) — 파티션 로컬 numbering(`evm.TU.{chrom}.1`이 파티션마다 반복)의 충돌 방지. `recombineEVM`이 이 출력을 `dedup_evm_gff3.py`에 넘김.

**로직 (검증본 `remerge.py` MODE=v2와 일치해야 함):**
1. 파티션 파싱: 4필드 행 → dir=필드4, basename에서 `(.+)_(\d+)-(\d+)$`로 lo/hi, offset=lo−1. 3필드 행(비분할 염색체) → dir=필드2, lo=1, hi=∞, offset=0.
2. core 구간: 염색체별 lo 오름차순 정렬 후 `core = [lo, 다음 파티션 lo − 1]`, 마지막(또는 비분할)은 `[lo, hi]`. (20 kb overlap 전제에서 `hi−20000`과 동치이나 overlap 값에 강건.)
3. 모델 적재: 파티션 gff3에서 gene/mRNA/기타(Parent 있는 exon·CDS) 행을 offset 보정해 모으고, mRNA마다 자식 세그먼트 span의 min/max로 lo/hi, `mid=(lo+hi)/2`, `ncds=자식 행 수`.
4. 소유 판정: 모델의 mid가 **자기 파티션** core 안에 있으면 채택.
5. complete 우선(v2): 채택된 각 모델에 대해 같은 (chrom,strand)에서 다른 파티션의 후보 중 span 겹침 ≥ 자기 span의 50%이고 `ncds`가 더 큰 것이 있으면 그것으로 교체. 교체 후 `(chrom,strand,lo,hi,ncds)` 키로 중복 제거.
6. 출력: gene(있으면, 파티션당 gene ID 1회) → mRNA → 자식 행 순. ID/Parent에 `{파티션basename}_` 접두. 정렬: (chrom, lo, hi) — 결정성.
7. v2 이웃 스캔은 lo-정렬 리스트 + bisect 윈도우(`lo ≥ m.lo − max_span`부터 `lo > m.hi`까지)로 — remerge.py의 전 리스트 스캔과 결과 동일하되 Ahy(2.5 Gb)급에서도 O(n·window).

- [ ] **Step 1: 실패하는 테스트 작성** — `bin/test_merge_evm_partitions.py`. 합성 파티션 두 개(`T_Chr1_1-1000000`, `T_Chr1_980001-1980000`, overlap 20 kb)와 비분할(`T_ChrC`)을 tmpdir에 만들어 검증:

```python
import os, subprocess, sys, tempfile, unittest

SCRIPT = os.path.join(os.path.dirname(__file__), "merge_evm_partitions.py")

def gff_model(chrom, gid, mid, cds_spans, strand="+"):
    lo = min(s for s, _ in cds_spans); hi = max(e for _, e in cds_spans)
    rows = [f"{chrom}\tEVM\tgene\t{lo}\t{hi}\t.\t{strand}\t.\tID={gid}",
            f"{chrom}\tEVM\tmRNA\t{lo}\t{hi}\t.\t{strand}\t.\tID={mid};Parent={gid}"]
    for s, e in cds_spans:
        rows.append(f"{chrom}\tEVM\texon\t{s}\t{e}\t.\t{strand}\t.\tID={mid}.e;Parent={mid}")
        rows.append(f"{chrom}\tEVM\tCDS\t{s}\t{e}\t.\t{strand}\t0\tID={mid}.c;Parent={mid}")
    return "\n".join(rows) + "\n"

class TestMerge(unittest.TestCase):
    def setUp(self):
        self.dir = tempfile.TemporaryDirectory()
        self.root = self.dir.name
        self.p1 = os.path.join(self.root, "T_Chr1", "T_Chr1_1-1000000")
        self.p2 = os.path.join(self.root, "T_Chr1", "T_Chr1_980001-1980000")
        self.pc = os.path.join(self.root, "T_ChrC")
        for d in (self.p1, self.p2, self.pc):
            os.makedirs(d)
        with open(os.path.join(self.root, "partitions_list.out"), "w") as fh:
            fh.write(f"T_Chr1\t{os.path.dirname(self.p1)}\tY\t{self.p1}\n")
            fh.write(f"T_Chr1\t{os.path.dirname(self.p2)}\tY\t{self.p2}\n")
            fh.write(f"T_ChrC\t{self.pc}\tN\n")

    def run_merge(self):
        out = os.path.join(self.root, "merged.gff3")
        subprocess.run([sys.executable, SCRIPT,
                        os.path.join(self.root, "partitions_list.out"), out],
                       check=True)
        return open(out).read()

    def write(self, part, text):
        with open(os.path.join(part, "evm.partition.gff3"), "w") as fh:
            fh.write(text)

    def test_core_ownership_and_offset(self):
        # p1 core = [1, 980000]; 로컬 5000-6000 -> mid 5500, p1 소유
        self.write(self.p1, gff_model("T_Chr1", "g1", "m1", [(5000, 6000)]))
        # p2에서 로컬 21000-22000 -> 게놈 1001000-1002000, p2 core = [980001, inf)
        self.write(self.p2, gff_model("T_Chr1", "g2", "m2", [(21000, 22000)]))
        self.write(self.pc, "")
        out = self.run_merge()
        assert "ID=T_Chr1_1-1000000_g1" in out
        assert "\t1001000\t1002000\t" in out            # offset 980000 적용
        assert "ID=T_Chr1_980001-1980000_g2" in out

    def test_boundary_model_kept_once(self):
        # 게놈 989000-991000 모델이 양 파티션에 등장(각각 로컬 좌표) -> mid 990000은 p2 core
        self.write(self.p1, gff_model("T_Chr1", "g1", "m1", [(989000, 991000)]))
        self.write(self.p2, gff_model("T_Chr1", "g1", "m1", [(9000, 11000)]))
        self.write(self.pc, "")
        out = self.run_merge()
        assert out.count("\tmRNA\t") == 1
        assert "T_Chr1_980001-1980000_m1" in out

    def test_v2_prefers_more_complete_neighbor(self):
        # p1 core가 소유한 2-CDS 모델, p2에 같은 자리(>=50% 겹침) 4-CDS 모델
        self.write(self.p1, gff_model("T_Chr1", "g1", "m1",
                                      [(970000, 972000), (974000, 976000)]))
        self.write(self.p2, gff_model("T_Chr1", "g1", "m1x",
            [(970000 - 980000 + 1 + 979999, 0)]))  # placeholder, 아래처럼 실제 로컬 좌표로
        # p2 로컬: 게놈 969000-976500 -> 로컬 = 게놈-980000
        self.write(self.p2, gff_model("T_Chr1", "g1", "m1x",
            [(-10999 + 980000 - 980000 + 970000 - 980000 + 980000, 0)]))
        # (실제 구현 시 로컬 좌표 -10000대가 되지 않도록 게놈 990000대 사례로 작성한다:
        #  owner p2, neighbor p1이 더 완전한 모델을 가진 대칭 사례가 산술이 단순하다)

    def test_deterministic_under_row_shuffle(self):
        text = (gff_model("T_Chr1", "g1", "m1", [(5000, 6000)])
                + gff_model("T_Chr1", "g2", "m2", [(8000, 9500), (9800, 9900)]))
        self.write(self.p2, ""); self.write(self.pc, "")
        self.write(self.p1, text)
        out_a = self.run_merge()
        lines = [l for l in text.splitlines() if l]
        self.write(self.p1, "\n".join(reversed(lines)) + "\n")
        out_b = self.run_merge()
        assert out_a == out_b

if __name__ == "__main__":
    unittest.main()
```

주의: `test_v2_prefers_more_complete_neighbor`는 위 스케치의 좌표 산술을 정리해 **owner=p2(core [980001,∞)), 게놈 990000대 locus, p1이 4-CDS·p2가 2-CDS**인 사례로 작성한다 — 음수 로컬 좌표가 나오지 않는다. 기대: 출력에 4-CDS 모델(태그 `T_Chr1_1-1000000_`) 1개만.

- [ ] **Step 2: 테스트 실행, 실패 확인**

```bash
micromamba run -n sylvan python bin/test_merge_evm_partitions.py
```
기대: `FileNotFoundError`/실패 (스크립트 없음).

- [ ] **Step 3: `bin/merge_evm_partitions.py` 구현** — 위 로직. 골격:

```python
#!/usr/bin/env python3
"""Merge per-partition EVM outputs at whole-mRNA granularity. (docstring: 근거/사용법)"""
import argparse, bisect, os, re, sys
from collections import defaultdict

PART_RE = re.compile(r"(.+)_(\d+)-(\d+)$")
ATTR_RE_CACHE = {}

def attr(attrs, key): ...            # remerge.py의 at()
def read_partitions(listing): ...    # -> [dict(name, dir, chrom, lo, hi, off)]
def assign_cores(parts): ...         # 염색체별 정렬, core_lo/core_hi
def load_models(part): ...           # evm.partition.gff3 -> 모델 리스트 (offset 적용)
def pick_owned(models, core_of): ... # mid in core
def prefer_complete(kept, models): ...  # v2: bisect 윈도우 스캔 + dedup
def emit(kept, out): ...             # (chrom,lo,hi) 정렬, gene 1회, 태그 접두

def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("partitions_list")
    ap.add_argument("outfile")
    args = ap.parse_args(argv)
    ...
```
파일이 없거나 빈 파티션은 건너뜀(경고 stderr). 파싱 불가 행은 스킵이 아니라 **stderr 경고 + 카운트** (EVM 출력은 신뢰 형식이므로 fail-fast까지는 불요). 요약 한 줄 stderr: `partitions N, models M, kept K, replaced R, dedup D`.

- [ ] **Step 4: 테스트 통과 확인**

```bash
micromamba run -n sylvan python bin/test_merge_evm_partitions.py
```
기대: 전부 PASS.

- [ ] **Step 5: 실데이터 대조** — 이전 세션 fullbest 산출물로 등가성 확인. fullbest의 파티션 gff3(`evm.gff3`)를 `evm.partition.gff3` 이름으로 심링크한 합성 partitions_list를 만들어 실행하고, mRNA 수를 검증본 `remerged_v2.gff3`와 대조:

```bash
SPAD=/data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-data-Sylvan-Sylvan/a0d6edfc-eacf-414a-858b-7f1059f900f1/scratchpad
awk -F'\t' '$3=="mRNA"' $SPAD/fullbest/remerged_v2.gff3 | wc -l   # 기대 기준값
# 합성 listing 생성 + 실행 + mRNA 수·좌표 md5 대조 (ID 접두 방식 차이는 필드1-8만 비교)
```
필드 1–8 (좌표·구조) 기준으로 동일해야 함. 차이가 나면 로직 불일치 — 수정 후 재대조.

- [ ] **Step 6: Commit**

```bash
git add bin/merge_evm_partitions.py bin/test_merge_evm_partitions.py
git commit -m "feat(evm): merge partitions at whole-mRNA granularity instead of per-row midpoints"
```

### Task 3: `combineEVMInputs` — Liftoff 승격 파일 이동 (+1.0) + miniprot cap 삽입 (+4.3)

**Files:**
- Modify: `bin/Snakefile_annotate:3385-3436` (rule combineEVMInputs)
- Modify: `config/config_annotate.yml` (combineEVMInputs 자원 + cap 파라미터 키)
- Modify: `toydata/config/config_annotate.yml`, `toydata/config/config_annotate_local.yml` (cap 파라미터 키 — 값은 동일 16/16이므로 기본값에 맡기고 **추가하지 않아도 됨**; 자원 키만 손대지 않음)

**Interfaces:**
- Consumes: `bin/cap_miniprot_depth.py` CLI — `cap_miniprot_depth.py <in> <out> --k-base K --k-intron K --tmpdir DIR` (이미 저장소에 존재, 검증 완료).
- Produces: `EVM/gene_predictions_checkCoords.gff3`에 Liftoff 행 포함, `EVM/transcript_alignments_checkCoords.gff3`에서 제외, `EVM/protein_alignments_checkCoords.gff3`는 miniprot cap 적용본. 출력 파일명·개수 불변 (DAG 불변).

- [ ] **Step 1: rule params에 cap 파라미터 추가**

```python
	params:
		results_dir = RESULTS_DIR,
		# miniprot depth cap: EVM accumulates weight per covering alignment with no
		# normalisation, so OrthoDB redundancy becomes vote count (297x mean depth on
		# covered bases). K=16 fixed by sweep (12: -0.91, 20: -0.04, 24: -0.22 IC F1).
		# 0 disables the cap.
		k_base = config_dict.get("combineEVMInputs", {}).get("miniprot_k_base", 16),
		k_intron = config_dict.get("combineEVMInputs", {}).get("miniprot_k_intron", 16)
```

- [ ] **Step 2: shell 블록 수정** — gene_predictions cat 목록에 `{input.liftoff}` 추가, transcript cat 목록에서 제거, protein 검증 후 cap 블록 추가:

```
	shell: """
		micromamba activate annotation
		# Liftoff rides in gene_predictions as ABINITIO_PREDICTION (weights file pairs
		# with this): evidence_modeler.pl:2601 lets only PROTEIN and ABINITIO add match
		# coverage, and as OTHER/TRANSCRIPT the most accurate track (chain F1 66.13 vs
		# AUGUSTUS 64.03) contributed nothing to coding scores.
		cat {input.helixer_gff3} \
			{input.augustus_gff3} \
			{input.geta_gff3} \
			{input.genewise_other} \
			{input.liftoff} \
			> {params.results_dir}EVM/gene_predictions.gff3
	
		bin/gff_to_evm.py {params.results_dir}EVM/gene_predictions.gff3 \
		--check-coords \
		--feature \\* > {params.results_dir}EVM/gene_predictions_checkCoords.gff3

		# Transcript evidence
		cat {input.gmapExon} \
			{input.pasa} \
			{input.st_star} \
			{input.pc_star} \
			> {params.results_dir}EVM/transcript_alignments.gff3

		bin/gff_to_evm.py {params.results_dir}EVM/transcript_alignments.gff3 \
			--check-coords \
			--feature \\* > {params.results_dir}EVM/transcript_alignments_checkCoords.gff3

		# Protein evidence
		cat {input.genewise_prot} \
			{input.miniprot} \
			> {params.results_dir}EVM/protein_alignments.gff3
		
		bin/gff_to_evm.py {params.results_dir}EVM/protein_alignments.gff3 \
			--check-coords \
			--feature \\* > {params.results_dir}EVM/protein_alignments_checkCoords.gff3

		# Cap miniprot depth before partitioning (global file, so partition overlap
		# regions keep a consistent retained set). GeneWise rows are exempt -- different
		# scale (1k rows vs 60M) and provenance. No errexit in these blocks, so every
		# step guards itself: a silent failure here would hand EVM the uncapped file.
		if [ "{params.k_base}" -gt 0 ] && [ "{params.k_intron}" -gt 0 ]; then
			awk -F'\t' '$2=="miniprot"' {params.results_dir}EVM/protein_alignments_checkCoords.gff3 \
				> {params.results_dir}EVM/protein_alignments_miniprot.gff3 || exit 1
			awk -F'\t' '$2!="miniprot"' {params.results_dir}EVM/protein_alignments_checkCoords.gff3 \
				> {params.results_dir}EVM/protein_alignments_other.gff3 || exit 1
			python bin/cap_miniprot_depth.py \
				{params.results_dir}EVM/protein_alignments_miniprot.gff3 \
				{params.results_dir}EVM/protein_alignments_miniprot.capped.gff3 \
				--k-base {params.k_base} --k-intron {params.k_intron} \
				--tmpdir {params.results_dir}TMP || exit 1
			cat {params.results_dir}EVM/protein_alignments_other.gff3 \
				{params.results_dir}EVM/protein_alignments_miniprot.capped.gff3 \
				> {params.results_dir}EVM/protein_alignments_checkCoords.gff3 || exit 1
			rm -f {params.results_dir}EVM/protein_alignments_miniprot.gff3 \
				{params.results_dir}EVM/protein_alignments_other.gff3 \
				{params.results_dir}EVM/protein_alignments_miniprot.capped.gff3
		fi
		"""
```

주의: awk 명령에 중괄호가 없으므로 snakemake `{{}}` 이스케이프 불요. `\t`는 파이썬 문자열에서 실제 탭으로 렌더돼 기존 rule들과 동일하게 동작.

- [ ] **Step 3: config 자원 증량** — `config/config_annotate.yml`의 `combineEVMInputs:` 섹션. 전역 cap은 7 GB/60M행 외부정렬 + alignment 요약 메모리(수 GB)이므로:

```yaml
combineEVMInputs:
  ncpus: 4
  memory: 32g               # was 4g: global miniprot cap holds per-alignment summaries
  time: "12:00:00"          #   (~6-8M Alignment objects) plus a 7 GB external sort
  miniprot_k_base: 16       # 0 disables the depth cap
  miniprot_k_intron: 16
```

(cluster_cmd 템플릿이 `time` 키를 쓰는 다른 룰과 동일한 형식인지 `pasaPost` 항목과 대조해 맞출 것.)

- [ ] **Step 4: 문법 확인** — dry-run 파스 통과:

```bash
cd toydata 2>/dev/null || true   # 실행은 Task 6에서; 여기선 파스만
micromamba run -n sylvan python -c "compile(open('bin/Snakefile_annotate').read(), 'S', 'exec')" 2>&1 | head -3
```
(Snakefile은 snakemake DSL이라 py-compile은 근사 검사 — 실제 확인은 Task 6 dry-run.)

- [ ] **Step 5: Commit**

```bash
git add bin/Snakefile_annotate config/config_annotate.yml
git commit -m "feat(evm): promote Liftoff to gene_predictions and cap miniprot depth at K=16 before partitioning"
```

### Task 4: `recombineEVM` — 모델 단위 병합으로 교체 + `runEVM` 자원 증량

**Files:**
- Modify: `bin/Snakefile_annotate:3544-3578` (rule recombineEVM)
- Modify: `config/config_annotate.yml` (runEVM memory/time)

**Interfaces:**
- Consumes: Task 2의 `bin/merge_evm_partitions.py` CLI, 컨테이너 annotation env의 `EVM_to_GFF3.pl`(EvmUtils, `recombine_EVM_partial_outputs.pl`과 같은 디렉터리 — 현행 rule이 bare 호출로 쓰고 있으므로 PATH에 있음).
- Produces: `EVM/EVM.all.gff3` (기존과 동일한 선언 output; ID 형식은 `{파티션태그}_evm.TU...`로 바뀜 — 다운스트림 `pasaPost`/`merge_pasa_evm.py`는 ID를 불투명 문자열로 다루므로 유일성만 요구).

- [ ] **Step 1: recombineEVM shell 교체**

```
	shell: """
			micromamba activate annotation
			export TMPDIR=./{params.results_dir}TMP
			export SLURM_TMPDIR=./{params.results_dir}TMP

			# genome.fasta was created by partitionEVM rule
			GENOME_FILE="{params.results_dir_abs}/genome.fasta"

			# Convert each partition's raw evm.out to GFF3 in partition-local
			# coordinates, then merge at whole-mRNA granularity. Stock
			# recombine_EVM_partial_outputs.pl assigns every gene/mRNA/CDS ROW to a
			# partition by row midpoint, so models straddling a boundary are torn
			# apart: of 87 reference multi-exon models crossing a boundary, 10
			# survived intact (Sn 11.5%) and 68 present in some partition vanished
			# from the merge. merge_evm_partitions.py assigns each mRNA (with its
			# gene and children) to one owner partition instead, preferring the more
			# complete copy where the 20 kb overlap holds two versions.
			while IFS="$(printf '\\t')" read -r chrom base flag part; do
				if [ "$flag" = "Y" ]; then d="$part"; else d="$base"; fi
				if [ -s "$d/evm.out" ]; then
					EVM_to_GFF3.pl "$d/evm.out" "$chrom" > "$d/evm.partition.gff3" || exit 1
				else
					: > "$d/evm.partition.gff3"
				fi
			done < {params.results_dir_abs}/partitions_list.out

			python bin/merge_evm_partitions.py \
				{params.results_dir_abs}/partitions_list.out \
				{output}.raw || exit 1

			# Remove phantom genes with coordinates beyond chromosome boundaries
			# (safety net for EVM partition merge issues)
			samtools faidx $GENOME_FILE
			python bin/dedup_evm_gff3.py {output}.raw ${{GENOME_FILE}}.fai {output}
			rm -f {output}.raw
		"""
```

주의: `while` 루프의 IFS 탭 — Snakefile 파이썬 문자열의 `\\t`가 셸에 `\t`로 전달되고 `printf '\t'`가 탭을 만든다. (`$'\t'`를 파이썬 문자열에 직접 쓰면 렌더 단계에서 실탭이 들어가 가독성이 나빠 printf 사용.)

- [ ] **Step 2: `EVM_to_GFF3.pl` 가용성 확인** — 컨테이너에서 1회 실행 (compute 노드 세션에서):

```bash
SIF=/data/gpfs/assoc/pgl/data/Deschampsia/sylvan/singularity/sylvan.sif
/apps/singularity/3.6.1/bin/singularity exec --cleanenv -B /data/gpfs $SIF \
  bash -c 'export MAMBA_ROOT_PREFIX=/opt/mamba; eval "$(/opt/mamba/bin/micromamba shell hook -s bash)"; micromamba activate annotation; command -v EVM_to_GFF3.pl && EVM_to_GFF3.pl 2>&1 | head -3'
```
기대: 경로 출력 + usage 메시지. 실패 시(모듈 부재 등) 대안: `$(dirname "$(command -v evidence_modeler.pl)")/EVM_to_GFF3.pl`를 `/opt/envs/genepred/bin/perl`로 호출 (fullbest가 쓴 방식) — rule을 그 형태로 바꾼다.

- [ ] **Step 3: runEVM 자원 증량** — Liftoff의 ABINITIO 승격으로 trellis 계산량 증가. 6h/16G에서 124중 4개 미완주, 24h/32G에서 124/124 완주:

```yaml
runEVM:
  ncpus: 1
  memory: 32g              # Liftoff-as-ABINITIO widens the trellis: 4/124 partitions
  time: "1-00:00:00"       # died under the old 6h/16G; 24h/32G completed 124/124
```
(Snakefile:3523의 `mem_mb` lambda가 `runEVM.memory`를 파싱해 재시도마다 +40 G — 32g 기준 retry1=72g로 노드 한도 내.)

- [ ] **Step 4: Commit**

```bash
git add bin/Snakefile_annotate config/config_annotate.yml
git commit -m "feat(evm): merge partitions per mRNA model in recombineEVM; give runEVM the 24h/32G the wider trellis needs"
```

### Task 5: 저장소 게이트 — 단위테스트 + pylint

**Files:** 없음 (검증만; 발견된 결함은 해당 Task 파일로 돌아가 수정)

- [ ] **Step 1: 신규·기존 단위테스트**

```bash
micromamba run -n sylvan python bin/test_cap_miniprot_depth.py
micromamba run -n sylvan python bin/test_merge_evm_partitions.py
micromamba run -n sylvan python bin/test_feature_importance.py   # CI의 실제 게이트
```
기대: 전부 PASS.

- [ ] **Step 2: pylint (CI 구성과 동일)**

```bash
micromamba run -n sylvan pylint bin/merge_evm_partitions.py bin/test_merge_evm_partitions.py \
  --disable=C0114,C0115,C0116 --disable=R0913,R0914,R0915 --disable=W0612,W0613 \
  --max-line-length=120 --exit-zero
```
advisory지만 E급(오류)은 고친다.

### Task 6: dry-run Job stats 동일성 (배포 전 필수 게이트)

**Files:** 없음

- [ ] **Step 1: baseline vs 새 Snakefile dry-run** — Ath_ODB 런 디렉터리에서, 반드시 프로덕션과 같은 트리거로:

```bash
RUN=/data/gpfs/assoc/pgl/data/Sylvan/sylvan_runs/Ath_ODB
BASE=<Task 0의 baseline_bin>/bin
cd $RUN
SYLVAN_CONFIG=config/config_annotate.yml \
  micromamba run -n sylvan snakemake -n -s $BASE/Snakefile_annotate \
  --rerun-triggers mtime --nolock --scheduler greedy 2>&1 | tee /tmp_scratch/dry_base.txt
SYLVAN_CONFIG=config/config_annotate.yml \
  micromamba run -n sylvan snakemake -n -s /data/gpfs/assoc/pgl/data/Sylvan/Sylvan/bin/Snakefile_annotate \
  --rerun-triggers mtime --nolock --scheduler greedy 2>&1 | tee /tmp_scratch/dry_new.txt
```
(정확한 인보케이션은 런의 `bin/annotate.sh`가 조립하는 snakemake 커맨드라인을 그대로 복제 — 실행 전에 `annotate.sh` 확인. `--snakefile`만 교체, config·cluster 인자 생략 가능 여부는 dry-run이므로 허용. login 노드 pulp 실패 시 `--scheduler greedy`가 회피.)

- [ ] **Step 2: Job stats 블록 diff** — 두 출력의 `Job stats:` 표가 **완전 동일**해야 하고, 재실행 목록에 완료 작업이 없어야 한다. baseline 대비 증감 0이 통과 기준. 다르면 원인 규명 전 진행 금지.

### Task 7: Ath_ODB 샌드박스 EVM 재실행 → ≈73.5 확인 (프로덕션 코드로 검증)

**Files:** 없음 (샌드박스는 세션 스크래치패드)

샌드박스 구성 — 라이브 런을 건드리지 않고, 프로덕션 룰 체인(combineEVMInputs→partitionEVM→writeEVMCommands→runEVM×126→recombineEVM)을 실제로 돌린다:

- [ ] **Step 1: 샌드박스 생성**

```bash
SB=<scratchpad>/evm_validation/Ath_ODB_sb
RUN=/data/gpfs/assoc/pgl/data/Sylvan/sylvan_runs/Ath_ODB
mkdir -p $SB/results
# 코드·설정: 새 repo bin + 런 config(+새 키는 기본값) + 컨테이너 심링크 + 엔트리 스크립트
cp -r /data/gpfs/assoc/pgl/data/Sylvan/Sylvan/bin $SB/bin
mkdir -p $SB/config $SB/singularity
cp $RUN/config/config_annotate.yml $RUN/config/alignAssembly.config $RUN/config/annotCompare.config $SB/config/
cp /data/gpfs/assoc/pgl/data/Sylvan/Sylvan/config/evm_weights.txt $SB/config/
ln -s $(readlink -f $RUN/singularity/sylvan.sif) $SB/singularity/sylvan.sif
ln -s $RUN/input $SB/input
# 업스트림 결과는 심링크(읽기전용 사용), EVM만 새로
for d in GETA TRANSCRIPT AB_INITIO LIFTOVER PROTEIN LIFTOFF HELIXER; do
  [ -e $RUN/results/$d ] && ln -s $RUN/results/$d $SB/results/$d
done
mkdir -p $SB/results/EVM/ok $SB/results/TMP $SB/results/logs
# convertToEvmFormat 산출물 5종은 재생성 대신 심링크 (코드 미변경 룰)
for f in genewise_protAln_evm.gff gmap_exonAln_evm.gff genewise_evm.gff3 liftoff_evm.gff3 miniprot_proteinAln_evm.gff; do
  ln -s $RUN/results/EVM/$f $SB/results/EVM/$f
done
```
(config의 evm_weights 경로가 `config/evm_weights.txt` 상대경로인지 확인해 새 weights가 잡히게 함. `results/` 하위 심링크 대상 디렉터리명은 실제 런 구조를 `ls $RUN/results`로 확인해 조정.)

- [ ] **Step 2: 타깃 실행** — 샌드박스에서 EVM.all.gff3만 타깃으로, sbatch 원샷:

```bash
cd $SB && sbatch -J evm_sb -c 2 --mem=8G -t 2-00:00:00 -A cpu-s1-pgl-0 -p cpu-s1-pgl-0 \
  --wrap 'bash bin/annotate.sh --rerun-triggers mtime --scheduler greedy --nolock -- results/EVM/EVM.all.gff3'
```
combineEVMInputs(cap ~1h) → partitionEVM → runEVM 126잡(대부분 수 분, cap 덕에 중앙값 ~4분) → recombineEVM. 모니터: `squeue`+`sacct`, 실패 잡은 `results/logs/`.

- [ ] **Step 3: 채점** — 검증본과 같은 reference·플래그:

```bash
GC=/data/gpfs/assoc/pgl/bin/miniconda3/bin/gffcompare
REF=<cds_odb>/ref_ath.cds.gff3   # resume.md: dd14dbde-*/scratchpad/cds_odb/ref_ath.cds.gff3
# CDS->exon 변환 후 채점 (fullbest merge_score.sh 방식 재현)
awk -F'\t' 'BEGIN{OFS="\t"} /^#/{next} NF<9{next}
  $3=="gene"||$3=="mRNA"||$3=="transcript"{print; next}
  $3=="CDS"{$3="exon"; print}' $SB/results/EVM/EVM.all.gff3 > $SB/evm_cds.gff3
$GC -r $REF $SB/evm_cds.gff3 --strict-match -e 0 -T --no-merge -o $SB/gc_sb
grep -E "Intron chain level|Transcript level|Locus level" $SB/gc_sb.stats
```
**통과 기준: IC F1 73.5 ± 0.3** (검증본 73.54; 전역 cap vs 파티션별 cap의 미세 차이 허용). 66~71대가 나오면 구현 누락 — cap 적용 여부(`protein_alignments_checkCoords.gff3`의 miniprot 행 수), liftoff 위치(gene_predictions 내 Liftoff 행 수), weights 파일 실제 사용본(`commands.*.list`의 `--weights` 경로) 순으로 원인 추적.

- [ ] **Step 4: 통과 시 결과 기록** — stats를 `ReviewerCommentsAnalysis/` 옆에 복사하지 말고 경로만 기록 (검증 산출물은 스크래치 유지; 논문 수치 원본은 기존 `ath_odb_sylvan_updated_cds.stats`).

### Task 8: 마무리 — 문서·이슈·보고

**Files:**
- Modify: `/data/gpfs/assoc/pgl/data/Sylvan/resume.md` (새 스냅샷: 반영 완료 + 검증 수치)

- [ ] **Step 1: resume.md 새 스냅샷 추가** (반영 4건 완료 상태, 검증 결과, fleet 배포는 다음 단계임을 명시 — transfrag 오염 런 재실행과 함께 진행 예정)
- [ ] **Step 2: 이슈 갱신** — wyim-pgl/Sylvan-EGAPx #66(Liftoff) #67(cap) #68(파티션 병합)에 반영 커밋 해시 코멘트
- [ ] **Step 3: 최종 보고** — 커밋 목록, 검증 수치, fleet 배포 시 주의(각 런 config에 새 키 없이도 기본값 동작; `--rerun-triggers mtime`이라 코드 배포만으로 EVM 재실행 안 됨 — transfrag 삭제 재실행 플랜과 함께 EVM 산출물 무효화 필요)

---

## 실행 중 수정 (2026-08-26)

- **Task 7 샌드박스는 snakemake가 아니라 새 rule 셸을 미러링한 3단계 SLURM 드라이버로 실행**
  (`scratchpad/evm_validation/stage{1,2,3}_*.sh`). 이유: Ath_ODB 런은 GETA 디버깅 여파로
  내부적으로 stale해서(dry-run 기준 174잡 재실행 대상) results/ 심링크 샌드박스 위에
  snakemake를 돌리면 GETA 체인 잡이 심링크를 뚫고 **라이브 런 디렉터리에 쓴다**.
  드라이버는 런의 전역 `*_checkCoords.gff3` 3종(fullbest 실험과 동일 세대, 2026-08-11)에서
  출발해 liftoff 행 이동 → 전역 cap → partition → 새 weights commands → runEVM 126 →
  새 recombine(모델 단위 병합) → gffcompare 채점까지 새 셸 로직을 그대로 수행한다.
  combineEVMInputs의 cat 단계 자체는 미실행이나, 행 이동은 check-coords가 행 단위
  처리라 cat 전/후 어디서 해도 동일하다(73.54 검증 자체가 checkCoords 수준 이동).
- **Task 4의 `EVM_to_GFF3.pl`은 bare 호출 불가로 판명** — 컨테이너 base perl에
  `URI::Escape` 부재(스모크 잡 6104869 실측, `gtf_to_alignment_gff3.pl`과 동일 함정).
  rule에서 `/opt/envs/genepred/bin/perl "$(command -v EVM_to_GFF3.pl)"`로 핀 (스모크
  6104887: rc=0, 288 mRNA 정상 변환).
- **Task 6 통과 실측**: baseline(51be201) vs 새 Snakefile dry-run — Job stats 블록 동일,
  reason 파일 나열 순서(집합 순회 비결정성)와 타임스탬프를 정규화하면 **diff 0줄**.

## Self-Review 결과

- **Spec coverage**: resume.md 4건+자원 증량 → Task 1(weights), 2+4(파티션 병합), 3(liftoff 이동+cap), 4(runEVM 자원). 검증 절차 2단계 → Task 6, 7. ✔
- **Placeholder 스캔**: Task 2 Step 1의 `test_v2_prefers_more_complete_neighbor` 스케치는 좌표 산술 정리 지시를 명시함 (구현 시점에 완성). 그 외 코드 블록은 실행 가능 수준. ✔
- **타입/이름 일관성**: `merge_evm_partitions.py` CLI가 Task 2와 Task 4에서 동일 (`<partitions_list> <outfile>`), 파티션 파일명 `evm.partition.gff3` 통일, config 키 `miniprot_k_base`/`miniprot_k_intron` 통일. ✔
- **알려진 잔여 위험** (plan 범위 밖, 보고에 포함): ① EVM.all.gff3 ID 형식 변화가 다운스트림(pasaPost·merge_pasa_evm)에 미칠 영향 — resume.md의 "PREFILTER 단계 검증" 후속 항목. ② 전역 cap의 메모리 실측 — Task 7에서 sacct MaxRSS 확인.
