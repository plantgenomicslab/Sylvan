# psiClass per-chromosome 병렬화 — 설계

날짜: 2026-07-22 · 대상: `bin/Snakefile_annotate` (`psiClass_STAR`, `prepPsiClass_STAR`), `bin/splitBam.py`
관련 이슈: wyim-pgl/Sylvan-EGAPx#25

## 문제
`psiClass_STAR`가 whole-bam을 `-p 48`(config threads=48)로 실행 → 고복잡도 영역
(Osa Chr8 2.8–4.3Mb: subexon 2,509 / atCnt 2.3M)에서 psiClass `classes` 단계의
**멀티스레딩 hang** → 3일 walltime TIMEOUT → 마커 미생성 → 컨트롤러 데드락.
psiClass는 EVM 필수(`combineEVMInputs.pc_star`)라 해당 런의 EVM을 완전 차단.
`prepPsiClass_STAR` 체크포인트가 크로모좀 분할을 하지만 psiClass_STAR가 **미사용**(설계 갭).

## 목표
psiClass를 **크로모좀별 단일스레드(threads=1) 잡의 병렬**로 재구성:
race hang 제거 + 크로모좀 간 병렬성 회복 + hotspot(Chr8) 격리(느려도 나머지 독립 완주).

## 설계

### 1. `bin/splitBam.py` — bin 구성 (개선)
- bam 헤더 @SQ에서 (seq, length) 수집.
- `THRESHOLD = 5_000_000` (5Mb, 인자로 조정 가능).
- `big` = length > 5Mb, `small` = length ≤ 5Mb.
- **big 없으면**: 전체 1 bin (`bin_all.bam`) — 병렬화 안 함(현행 유지).
- **big 있으면**:
  - `T = min(length of big)` — bin 목표 크기.
  - big 각각 → 개별 bin: `bin_{sanitized_name}.bam`.
  - small → greedy bin-packing(내림차순, 각 bin 합 ≤ T) → `bin_small_{K}.bam`.
- **순수 함수 `compute_bins(seqs) -> list[(binid, [seqnames])]`** 로 분리(단위 테스트 대상).
  I/O(samtools view -b)는 별도.
- 출력: `splitBam/bin_*.bam` (+ 각 bin의 seq 목록은 파일명/매핑으로).

### 2. `rule psiClass_bin` (신규, per-bin)
- input `splitBam/bin_{binid}.bam` → `samtools view -bq 10`(MAPQ≥10) → `psiclass -p 1 -c 0.05`
  → `splitBam/{binid}/psi_vote.gtf`. threads=1(멀티스레딩 hang 회피).

### 3. `rule psiClass_STAR` (aggregation으로 재구성)
- 체크포인트 output에서 bin 동적 수집(`glob_wildcards`).
- 각 bin vote.gtf의 gene_id/transcript_id에 **binid 접두어** 부여(ID 충돌 방지) 후 concat.
- `gtf_to_alignment_gff3.pl`(genepred perl) → `psiclass.STAR_vote.gff` (소비처 combineEVMInputs 불변).

### 4. config (`config_annotate.yml`)
- `psiClass_bin`: ncpus 2, threads 1, memory 64g (per-chromosome, 병렬 다수 → 노드당 여러개).
- `psiClass_STAR`(aggregation): ncpus 1, memory 8g (cat+convert만).
- `prepPsiClass_STAR`: 기존 유지(split I/O).

## 효과
- 총 시간 ≈ max(크로모좀별 시간) (기존 whole-bam hang 대비 대폭 단축).
- Chr8 hotspot은 자기 bin에 격리 — single-thread라 hang 안 함(느릴 뿐), 나머지 병렬 완주.
- memory를 bin당 64g로 낮춰 노드당 다수 bin 동시 실행(기존 192g는 노드 독점).

## 트레이드오프 / 후속
- Chr8 bin의 single-thread walltime 내 완주는 미검증 → 미완주 시 #25 옵션2(초고-subexon
  영역 복잡도 가드)를 그 bin에 추가(별도 작업).
- ID 접두어 방식은 psiclass 출력 네이밍 확인 후 확정.

## 검증
- `compute_bins` 단위 테스트(각 게놈 크기 분포로 bin 구성·균형 확인).
- Osa로 dry-run(체크포인트→per-bin→aggregation DAG 형성) 후 실 적용.
