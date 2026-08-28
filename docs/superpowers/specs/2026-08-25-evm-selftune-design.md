# EVM reference-free self-tuning — 설계

작성 2026-08-25 · 상태: **Codex 검토 반영 완료, Phase 1 구현 대기** · 대상 repo `plantgenomicslab/Sylvan`

## 1. 문제

CDS 기준 gffcompare(리뷰어 R1-B 플래그 `--strict-match -e 0 -T --no-merge`)에서
Sylvan-ODB가 BRAKER3-ODB에 intron-chain F1 8/8 전패한다.

| 종 | Sylvan-ODB IC F1 | BRAKER-ODB IC F1 |
|---|---:|---:|
| Ath | 66.7 | 73.0 |
| Osa | 43.0 | 50.7 |
| Mtr | 45.1 | 52.0 |
| Ptr | 63.1 | 72.2 |
| Sly | 33.0 | 36.2 |
| Dca | 46.5 | 52.9 |
| Svi | 53.2 | 64.6 |
| Spe | 43.5 | 51.3 |

sensitivity가 8/8 낮고 precision은 3종(Osa·Sly·Svi)에서만 우위다.

## 2. 원인 귀속 — 실측으로 확정

동일 reference·동일 플래그로 세 단계를 각각 채점했다
(`scratchpad/stage_attribution.sh`, 결과 `scratchpad/stage_attr/`).

intron-chain Sn:

| 종 | EVM | PREFILTER | FILTER | BRAKER-ODB | EVM 시점 격차 |
|---|---:|---:|---:|---:|---:|
| Ath | 60.3 | 61.3 | 61.3 | 67.2 | −6.9 |
| Osa | 37.5 | 39.8 | 39.7 | 44.1 | −6.6 |
| Mtr | 37.3 | 38.5 | 38.3 | 42.1 | −4.8 |
| Ptr | 57.4 | 59.1 | 58.9 | 67.5 | −10.1 |
| Sly | 30.0 | 30.2 | 30.1 | 33.0 | −3.0 |
| Dca | 53.3 | 53.8 | 53.8 | 58.1 | −4.8 |
| Svi | 49.1 | 51.8 | 51.7 | 56.5 | −7.4 |
| Spe | 38.0 | 39.0 | 38.7 | 44.1 | −6.1 |

**결론: BRAKER와의 격차는 EVM 단계에서 이미 전부 발생한다.** 후단은 8/8 전부
sensitivity를 +0.1~+2.6 올린다.

FILTER는 병목이 아니다. gene을 Svi 25.9%·Osa 20.5% 제거하면서도 intron-chain Sn
손실은 −0.1이고, locus precision은 크게 오른다(Svi 45.0→57.5, Osa 47.1→53.4).
버려지는 것은 애초에 reference와 매칭되지 않던 모델이다.

> gene count 감소를 근거로 FILTER attrition을 주원인으로 본 초기 판단은 틀렸다.
> count는 20% 줄어도 chain accuracy는 유지된다.

## 3. EVM 내부의 지배 구조 — 소스로 확정

컨테이너의 EVM은 공식 릴리스 그대로다(`singularity/Sylvan.def:492-515`가
tarball을 받아 `tar zxf` + `make`만 수행, 로컬 패치 없음). 따라서 공식
v2.1.0 소스 분석이 컨테이너에 그대로 적용된다.

### 3.1 점수는 염기 단위로 누적된다

`EvmUtils/evidence_modeler.pl`:

```perl
# add_match_coverage
for (my $i = $end5; $i <= $end3; $i++) {
    $CODING_SCORES[$i] = $current_coding_score + $weight;
}
# add_introns
$INTRONS_TO_SCORE{$intron_key} += $intron_score;   # intron_score = weight × intron 길이
```

source별 cap·평균·saturation·redundancy 정규화가 **없다**. 같은 염기를 N개
alignment가 덮으면 `N × weight`가 쌓인다.

### 3.2 실측 depth (Ath_ODB, `Ath_Chr1:1-1,000,000`)

| 지표 | 값 |
|---|---:|
| miniprot segment | 654,453 |
| 고유 Target(query protein) | 124,450 |
| GeneWise row | 1,075 |
| 평균 depth (1Mb 전체) | 166.5× |
| 평균 depth (covered base) | 297.4× |
| 최대 depth | 2,316× |

miniprot weight 2의 per-base 기여는 covered base 평균 약 595, 최대 4,632다.
AUGUSTUS 단일 모델의 weight 7과 같은 척도가 아니다.

**가중치 스칼라로는 보정 불가능하다.** AUGUSTUS를 400 이상으로 올려야 균형이
맞는데 정상적 튜닝 범위가 아니고, 고품질 단일 정렬 locus와 OrthoDB homolog가
수천 개 겹친 locus를 하나의 스칼라로 동시에 보정할 수 없다.

### 3.3 클래스별 기여 경로 (설계에 영향)

```perl
unless ($ev_class =~ /^(PROTEIN|ABINITIO_PREDICTION)$/) { return; }   # add_match_coverage
```

- per-base coding score: **PROTEIN(miniprot, GeneWise) + ABINITIO(AUGUSTUS, Helixer)만**
- OTHER_PREDICTION(Liftoff, GETA, Genewise) / TRANSCRIPT(PASA, StringTie, PsiClass, gmap):
  exon·intron 점수에만 기여
- intergenic 점수: ABINITIO만

즉 weights에서 가장 높은 PASA(TRANSCRIPT, 10)와 GETA(OTHER_PREDICTION, 5)는
base coding score에 한 점도 기여하지 않는다.

### 3.4 Liftoff 카테고리는 버그가 아니다

`evidence_modeler.pl:880`이 ev_class를 입력 파일이 아니라 weights의
`ev_type → ev_class` 매핑에서 결정한다. `OTHER_PREDICTION Liftoff 2` 선언대로
정상 처리되며(`evm.out` 137,330건과 일치), lookup 실패 시엔 경고 후 skip(892행)
또는 `confess`(900행)다. 조용한 실패 경로가 없다. **튜닝 대상에서 제외.**

## 4. 목표와 비목표

**목표**

- reference를 목적함수에서 완전히 배제한 self-tuning을 파이프라인 기능으로 만든다.
- 새 게놈에서 그 종의 RNA-seq만으로 EVM 입력·가중치를 스스로 조정한다.
- 조정 결과를 사후에 reference로 독립 검증한다.

**비목표**

- 논문 수치 최적화(사용자 확정: 목적은 파이프라인 실질 품질).
- FILTER RF 재설계 — 현재 잘 작동하므로 이번 범위 밖.
- Liftoff 카테고리 변경 — 버그가 아님.

## 5. Phase 1 — miniprot 입력 정규화 (최우선)

### 5.1 대상

`results/EVM/miniprot_proteinAln_evm.gff` — Ath_ODB 기준 **7.08 GB / 60,038,290 행**.
`rule convertToEvmFormat`(Snakefile_annotate:3163)의 출력이며
`rule combineEVMInputs`(:3216)가 `protein_alignments.gff3`로 cat한다.

행 형식:

```
Ath_Chr5  miniprot  miniprot  13203494  13203829  530  +  0  ID=MP000001;Rank=1;Identity=0.8661;Target=AL1006U10010.t1 1 112
```

필터에 쓸 수 있는 필드가 모두 있다: score(6열), `Rank`, `Identity`, `Target`.

### 5.2 Rank 필터는 무효 — 측정으로 기각

파티션 실측 Rank 분포:

| Rank | 행 수 | 누적 |
|---:|---:|---:|
| 1 | 527,763 | 80.6% |
| 2 | 92,009 | 94.6% |
| 3 | 23,908 | 98.3% |
| ≥4 | 10,773 | 100% |

`Rank≤1`로 잘라도 20%만 준다. **중복은 순위가 아니라 서로 다른 query protein
124,450개가 같은 locus에 몰리는 데서 온다.** OrthoDB의 데이터베이스 중복도와
계통별 단백질 수가 그대로 투표수가 되는 구조다.

### 5.3 입력 실태 — 설계 전제를 바꾸는 실측

구현 전에 세 가지를 실측했고, 초안의 전제 두 개가 틀렸음이 드러났다.

| 확인 항목 | 실측 결과 | 설계에 미치는 영향 |
|---|---|---|
| 좌표 정렬 | **정렬 안 됨.** `Chr5:13203494 → 13201546 → 17233752 → … → Chr2`로 seqid까지 섞임 | "seqid별 청크 스윕" 불가 → external sort 필요 |
| Identity 단위 | **segment별.** multi-segment ID 29,791개 중 **29,498개(99%)에서 변동** | 행별 identity 필터는 alignment를 조각냄 → alignment-level 집계 필수 |
| ID 범위 | 단일 seqid·단일 strand에 고정(500k행 검사, 위반 0) | alignment ID를 원자 단위로 써도 안전 |

용어를 분리한다. 초안이 `Target 원자성`과 `동일 ID`를 혼용했으나 둘은 다르다.

```
segment    GFF 한 행
alignment  동일 ID를 가진 모든 segment — 선택의 원자 단위
query      Target의 첫 토큰 — 한 query가 Rank 1,2,… 여러 alignment를 가짐
```

원자 단위는 **alignment ID**다. query 전체를 원자적으로 다루면 한 단백질의 여러
genomic hit이 함께 유지/제거되어 cap의 locus 선택 기능을 잃는다.

### 5.4 채택 설계 — per-base + exact-intron 이중 cap

EVM의 두 누적 경로에 각각 대응하는 제약을 건다.

```
add_match_coverage  → 염기별 weight 누적   → per-base coding depth cap
add_introns         → intron_key별 누적    → exact intron multiplicity cap
```

두 제약은 동치가 아니다. 서로 다른 splice boundary를 가진 alignment들이 같은 exon
body를 덮을 수 있으므로 별도 파라미터로 둔다.

```
각 genomic base b:                    retained coding segment 수 ≤ K_base
각 exact intron (seqid,donor,acceptor,strand): 지지 alignment 수 ≤ K_intron
```

윈도우 평균 cap은 쓰지 않는다. 평균이 K 이하여도 좁은 exon이나 splice boundary에서
국소 폭증을 허용해 EVM이 실제로 보는 값을 제한하지 못한다.

#### alignment-level 품질 지표

segment별 값을 그대로 쓸 수 없으므로 ID 단위로 집계한다.

```
aligned_aa  = Target 좌표 union 길이
identity_w  = Σ(seg_identity × seg_aa) / Σ(seg_aa)     # 길이 가중 평균
score_total = Σ(seg_score)
```

우선순위(사전식):

```
1. identity_w   — 길이 정규화 품질
2. aligned_aa   — 짧은 보존 domain hit이 capacity를 독점하지 않도록
3. Rank         — 같은 query 안에서만 유효한 순위 (query 간 비교에 쓰지 않음)
4. score_total  — raw score는 정렬 길이 비례라 후순위
5. ID 문자열    — 결정적 tie-break
```

identity 하한(`min_identity`)도 **alignment-level `identity_w`에 적용**한다.
segment별로 적용하면 낮은-identity exon만 빠진 chimeric evidence가 만들어진다.

#### 선택 알고리즘 — component 내 우선순위 admission

"과밀 구간에서 낮은 순위부터 제거"는 순회 순서에 의존하고 긴 alignment를 과도하게
제거할 수 있다. 대신 admission 방식을 쓴다.

```
overlap component(연결된 겹침 덩어리) 안에서:
  품질 내림차순으로 alignment를 검토
  그 alignment의 모든 coding segment와 exact intron을 추가해도
  K_base·K_intron을 넘지 않으면 통째로 accept, 넘으면 reject
```

고품질이 먼저 capacity를 예약하고, ID 원자성이 자연히 보장되며, 제약 만족을
불변식으로 테스트할 수 있다. 전역 최적은 아니므로, **작은 component(수십 alignment
이하)에서 exact solver와 대조해 greedy 손실을 측정**한다 — reference 없이 검증
가능한 항목이다.

`min_identity` 참고용 실측 분포(파티션, segment 기준):

| Identity | 행 수 | 비율 |
|---|---:|---:|
| ≥0.9 | 89,584 | 13.7% |
| 0.8–0.9 | 107,205 | 16.4% |
| 0.7–0.8 | 111,077 | 17.0% |
| 0.5–0.7 | 211,655 | 32.3% |
| <0.5 | 134,932 | 20.6% |

identity 단독 필터는 원격 상동 단백질을 통째로 날려 evidence 다양성을 해치므로
**주 수단은 depth cap**, identity는 우선순위와 하한에만 쓴다. `min_identity`는
기본값으로 굳히지 않고 K와 독립적인 sweep 축으로 유지한다.

### 5.5 구현

3-pass + external sort 구조로 7 GB를 상수 메모리에 가깝게 처리한다.

```
Pass 1  sort -k1,1 -k<ID> 로 동일 ID를 인접시킴
        → 스트리밍으로 alignment 요약 생성 (ID, seqid, strand, span, identity_w,
          aligned_aa, score_total, Rank, intron key 목록)
Pass 2  요약을 (seqid, start)로 sort → overlap component 단위 admission
        → accept된 ID 집합 산출
Pass 3  원본을 다시 읽어 accept ID의 segment만 출력
        → EVM parser가 요구하는 형식·순서로 재작성
```

- 새 모듈 `bin/cap_miniprot_depth.py`. 외부 정렬은 GNU `sort`에 위임(`-S`, `-T`는
  GPFS scratch 지정 — compute 노드 `/tmp`는 tmpfs라 RAM을 먹는다).
- **GeneWise는 cap 대상에서 제외**한다. 파티션당 1,075행으로 규모가 다르고
  provenance가 달라, 같은 capacity pool에 넣으면 miniprot을 줄이는 과정에서
  함께 제거될 수 있다. GeneWise의 상대적 강화는 Phase 3의 weight 축에서 다룬다.
- cap 단위는 EVM 내부 스코어 컨텍스트와 맞춰 `(seqid, strand)`로 계산한다.
- 파이프라인 편입: `convertToEvmFormat`과 `combineEVMInputs` 사이에 새 rule
  `capMiniprotDepth`. 기본 비활성(`Internal.miniprot_depth_cap: null`)으로 두어
  기존 15런의 DAG를 바꾸지 않는다.
- 단위 테스트(불변식):
  - 모든 base에서 depth ≤ K_base
  - 모든 exact intron에서 multiplicity ≤ K_intron
  - accept된 ID의 segment가 하나도 빠지지 않음(원자성)
  - 입력 순서를 섞어도 출력이 동일(결정성)
  - malformed·ID 누락·Identity 결측 행 처리 정책

### 5.6 Phase 1 효과 확인 (reference 미사용)

depth 통계로 확인 **가능한 것**과 **불가능한 것**을 구분한다.

가능:

1. 구현이 의도한 capacity 제약을 실제로 만족하는가 (불변식 검사)
2. miniprot의 수치적 지배가 얼마나 줄었는가 — depth 분포 전후
3. 제거가 특정 염색체·identity/coverage 계층·repeat 구간에 편향되는가
4. EVM 실행시간·메모리 개선
5. 산출 gene/intron/coding-base 수가 붕괴하지 않는가

불가능: **depth 정상화만으로 annotation 품질 개선을 결론지을 수 없다.**
높은 depth는 잘 보존된 진짜 유전자에서도 발생한다.

기록할 통계는 segment depth뿐 아니라 alignment 수, 고유 query 수,
identity/coverage 분포, intron multiplicity를 함께 남긴다.

#### K 후보 범위의 물리적 근거

miniprot weight가 2이고 주요 ab initio weight가 7이므로, `K_base=4`면 miniprot의
coding 기여가 8이 되어 단일 AUGUSTUS 모델과 같은 규모가 된다. 이를 기준으로
로그 스케일 후보를 둔다.

```
K_base   : 2, 4, 8, 16, 32
K_intron : 2, 4, 8, 16
```

`K≥64`는 기여가 128 이상이라 정규화 효과가 약하고, `K=1`은 ortholog 다양성을
지나치게 잃을 위험이 있다. 이는 최적값의 근거가 아니라 **탐색 범위의 근거**다.

Phase 1에서는 K를 확정하지 않는다. 구현 검증 + 무의미한 K 제거 + Phase 3에 넘길
후보 생성까지가 범위이며, 확정은 Phase 3에서 evaluator로 한다.

## 6. Phase 2 — reference-free evaluator

### 6.1 기존 코드를 재사용하지 않는 이유

`bin/refine_boundaries.py`가 이미 splice junction을 다루지만 목적함수로는 부적합하다:

- `load_splice_junctions(splice_path, tolerance=5)` — **tolerance 인자가 본문에서 미사용**.
- `score_splice_support()` — donor/acceptor 각각 ±5를 허용해 **intron당 121개 좌표 조합**을
  정답 처리. alternative splice site를 같은 것으로 합쳐 평가를 부풀린다.
- 로그 `n_junctions // 11` — junction이 `(chrom,strand)`와 `(chrom,".")` 두 키에
  중복 저장되므로 총합은 실제의 2배인데 11로 나눈다. 실제 228,314행이 41,511로 출력된다.
- 첫 isoform만 평가, 대규모 `set(range(...))` 기반 coverage.

따라서 별도 경량 evaluator를 만든다. 재사용은 `get_introns()` 골격 정도로 한정한다.

### 6.2 설계

```
입력: 예측 GFF(EVM.all.gff3), STAR SJ.out.tab 집합
고신뢰 junction: uniq reads ≥ t, 최소 m개 샘플에서 재현, canonical motif
observability mask: 양쪽 exon flank의 최소 coverage + unique split-read 하한
                    (exon body coverage만으로는 그 intron의 검출 가능성을 보장 못 함)
지표: tolerance 0 이 primary
      junction recall    = |I ∩ J_eval| / |J_eval ∩ observable|
      junction precision = |I ∩ J_eval| / |I ∩ observable|
      junction-level micro F1 + gene-level macro F1 병기
      ±1/±2 는 diagnostic 으로만 기록
```

### 6.3 이 지표는 "held-out test"가 아니라 RNA-consistency objective다

초안은 RNA-seq 샘플을 `J_train`/`J_test`로 나눠 evidence overfitting을 차단한다고
했으나, §10이 evidence 생성에는 전체 샘플을 유지한다고 명시하므로 모순이었다.
전체 샘플이 PASA/StringTie/PsiClass 입력에 쓰이면 어느 junction도 이미 후보 모델에
영향을 준 상태다. junction 테이블만 나눠서는 독립성이 생기지 않는다.

진짜 held-out을 만들려면 평가용 샘플을 **그 tuning run의 모든 RNA-derived evidence
생성에서 제외**해야 하고, 이는 PASA·StringTie·PsiClass 재실행을 뜻해 Phase 1–3의
반복 루프에서는 비현실적이다.

따라서 이 지표를 **RNA-consistency objective**로 정직하게 규정한다. 후보 가중치·K를
비교하는 상대 신호로는 유효하지만, 독립적 일반화 성능으로 해석하지 않는다. 독립
평가는 Phase 5의 reference 채점이 담당한다.

샘플 분할은 다음 목적으로만 쓴다 — 지표의 분산 추정과 안정성 확인:

- 11개 샘플에 고정 8/3 분할은 쓰지 않는다(평가 3개는 조직 특이 junction의
  재현성 기준이 불안정).
- **stratified 7/4 split을 여러 번 반복**하고, 최종 한 번 6/5 swap으로 확인한다.
- 동일 생물학적 조건의 replicate는 같은 fold에 묶는다(흩뜨리면 junction 중복으로
  leakage가 커진다).
- "최소 m개 샘플 재현" 기준은 평가 표본 수에 맞춰 조정한다(평가 4개면 m=2 + read
  support 임계 병용).

### 6.4 보조 제약 (단일 지표 과최적화 방지)

splice F1만 최적화하면 다중 exon·고발현 유전자만 편애한다. 다음을 penalty로 둔다:

- complete ORF 비율, 내부 stop/frameshift 비율
- 인접 유전자 fusion 비율, exon 수 폭증, 과도하게 긴 유전자
- single-exon gene 수 급변, coding base 총량 급변

BUSCO/OMArk는 목적함수에서 제외하고 사후 QC로만 쓴다(reference-free 일관성 유지,
OMArk는 이슈 #23의 빈 결과 문제도 있음).

## 7. Phase 3 — stratified partition sweep

파티션은 이질성이 크다(Ath_ODB 126개 실측):

| 지표 | min | p10 | median | p90 | max |
|---|---:|---:|---:|---:|---:|
| gene-prediction CDS | 0 | 1,748 | 4,703 | 5,638 | 6,079 |
| transcript rows | 0 | 1,939 | 10,491 | 12,758 | 14,455 |
| miniprot rows | 0 | 85,283 | 547,181 | 717,759 | 915,733 |

따라서 무작위 10~20개는 위험하다.

```
1. 파티션별 feature 계산: 염색체, gene/CDS density, trusted junction density,
   RNA coverage·샘플 수, miniprot segment 수·depth, repeat fraction, GC/N,
   single/multi-exon 비율
2. 표준화 후 stratified sampling — 염색체별 최소 1개, gene density·repeat·RNA·
   miniprot depth의 low/mid/high strata, 극단값 의도적 포함
3. tune 16~20개 / validation 16~20개로 분리
4. 탐색 축: K(depth cap), min_identity, EVM 가중치(AUGUSTUS·Helixer·GETA·PASA·miniprot)
   successive halving으로 후보 축소
5. 20 kb overlap 영역은 한 번만 평가(경계 유전자·junction 이중집계 방지)
6. bootstrap으로 순위 안정성 확인 — 평균 점수뿐 아니라 승률, 염색체 holdout,
   worst-stratum 성능을 함께 본다
```

최종 후보 2~3개만 더 넓은 40~60개 또는 전체로 확대한다.

## 8. Phase 4 — 전체 적용

확정된 (K, min_identity, weights)로 전체 EVM 재실행 후
PASA → `merge_pasa_evm` → `refine_boundaries` → AGAT → PREFILTER → FILTER를 돌린다.

**raw EVM에서 개선됐으나 최종에서 사라지는 경우를 분리 측정한다.** 실제 순서상
PASA와 refine이 구조를 바꾸며(Ath_ODB에서 truncated candidate 9,907 중 531개 실제 교체),
EVM 가중치의 raw 효과가 가려지거나 역전될 수 있다. 두 지점 모두에서 Phase 2 지표를
측정해 EVM tuner 문제와 downstream interaction 문제를 구분한다.

## 9. Phase 5 — reference 사후 검증

CDS + R1-B 플래그로 gffcompare, BRAKER-ODB와 대조.
목적함수가 reference를 한 번도 보지 않았으므로 이 수치는 독립 평가로 유효하다.

기준선(현재): §2 표의 FILTER 열.

## 10. 위험

| 위험 | 완화 |
|---|---|
| depth cap이 진짜 보존 영역의 evidence까지 깎음 | K를 Phase 3에서 데이터로 결정, identity 하한은 낮게(0.5)·sweep 축으로 유지 |
| 반복서열의 고identity 짧은 domain hit이 capacity 독점 | coverage·repeat stratum별 제거율을 통계로 감시 |
| 거의 동일한 ortholog가 상위 K 독점(계통 다양성 상실) | Phase 1은 모니터링만, family cap은 메타데이터 확보 시 도입 |
| K 경계에서 결과가 불연속 | K 주변 후보를 촘촘히, 결정적 tie-break로 순서 독립성 확보 |
| 서브셋 최적이 전체와 불일치 | stratified + validation 분리 + bootstrap 승률 |
| 7 GB 처리 성능·메모리 | external sort 3-pass, `sort -T`는 GPFS scratch(compute `/tmp`는 tmpfs=RAM) |
| cap 후 EVM parser 형식 위반 | Pass 3에서 형식·순서 재작성 + 실제 파티션 EVM 실행으로 검증 |
| 새 rule이 기존 15런 DAG를 흔듦 | 기본 비활성 + `--rerun-triggers mtime` dry-run으로 baseline 대조 |

## 11. 검토로 해소된 결정사항

Codex 검토(2026-08-25, gpt-5.6-sol) + 실측으로 확정한 항목.

| 질문 | 결론 | 근거 |
|---|---|---|
| depth cap 단위 | per-base cap + exact-intron cap **이중**. 윈도우 평균은 국소 폭증을 허용해 부적합 | EVM의 두 누적 경로(`add_match_coverage`/`add_introns`)에 1:1 대응 |
| 우선순위 기준 | raw score 아님. **length-weighted identity → aligned_aa → Rank → score_total** | score·Identity가 segment별 값이고, 99% ID에서 Identity가 변동 |
| GeneWise | Phase 1 cap **제외** | 1,075행 규모·다른 provenance. 상대적 강화는 Phase 3 weight 축에서 처리 |
| 샘플 분할 | 고정 8/3 기각. **반복 7/4 stratified + 6/5 swap**, replicate는 같은 fold | 평가 3개는 재현성 기준이 불안정 |
| Phase 1 신뢰 범위 | 구현 검증·지배 완화·편향 점검까지. **품질 결론은 불가** | 높은 depth는 진짜 보존 유전자에서도 발생 |
| K 탐색 범위 | K_base 2·4·8·16·32 / K_intron 2·4·8·16 | K=4에서 miniprot 기여 8 ≈ AUGUSTUS 단일 모델(7) |

초안에서 **수정된 결함**:

1. `Target 원자성` → `alignment ID 원자성` (§5.3)
2. 행별 identity 필터 → alignment-level `identity_w` (§5.4)
3. "좌표 정렬 입력" 전제 → external sort 3-pass (§5.5)
4. per-base 단일 cap → per-base + exact-intron 이중 cap (§5.4)
5. "과밀 구간에서 낮은 순위 제거" → component 내 우선순위 admission (§5.4)
6. `J_test` held-out 주장 → RNA-consistency objective로 재규정 (§6.3)

## 12. 남은 열린 항목

- **query coverage**: `Target` 필드에 정렬 구간은 있으나 query 전체 길이가 없다.
  `input/proteins/all_proteins.fa`(12.4M seq)에서 길이 테이블을 만들지, `aligned_aa`
  대용으로 갈지 Phase 1 구현 중 비용을 보고 결정한다.
- **partition 경계 일관성**: 같은 genomic region이 서로 다른 파티션에서 독립적으로
  cap되면 경계에서 retained evidence가 달라진다. cap을 파티션 이전(전역 GFF)에
  적용하므로 현재 설계에서는 발생하지 않으나, Phase 3의 파티션 sweep에서
  재확인한다.
- **phylogenetic diversity**: per-base cap은 투표 수만 제한하고 거의 동일한 ortholog가
  상위 K를 독점하는 것을 막지 못한다. `Target` 문자열로 family를 안전하게 복원할 수
  없어 Phase 1에서는 다루지 않고, 통계로 모니터링만 한다.
- **greedy 손실**: 작은 overlap component에서 exact solver와 대조해 정량화한다.
