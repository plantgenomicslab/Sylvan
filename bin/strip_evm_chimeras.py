#!/usr/bin/env python3
"""EVM 출력에서 키메라 모델을 걸러 PASA 로 넘길 입력을 만든다.

## 왜 필요한가

EVM 은 증거가 밀집한 구간에서 인접 유전자의 증거를 연쇄적으로 이어붙여 하나의
거대 모델을 만든다. 개별 증거는 멀쩡하다 — `Sly` 문제 구간의 증거 최대 span 은
AUGUSTUS 49.7 kb / Liftoff 45.4 kb / miniprot 4.8 kb 인데, EVM 이 낸 모델은
**188 kb 에 575 exon, 38,009 aa** 였다.

이런 모델은 PASA 의 alignment assembler 를 죽인다 (2026-07-30, `pasaPost.Sly`):

    Fatal Error reading input file.
    PASA died on input file.  See pasa_killer.input
    Thread 4 terminated abnormally: Died at PASA_alignment_assembler.pm line 262

`pasa_killer.input` 은 `evm.model.Sly_11.1097`(575 exon, `-` strand)과 그 안의
반대 strand transcript `asmbl_14097` 를 지목했다. 입력이 변하지 않으므로 재시도는
**같은 지점에서 7 시간을 다시 태우고 죽는다** — psiclass(#41)와 같은 무한 루프다.

## 기준과 그 근거

    exon > 110  또는  단백질 > 20,000 aa

두 축을 모두 쓰는 이유는 PASA 가 정확히 무엇에 걸려 죽는지 모르기 때문이다.
2026-07-30 fleet 실측(EVM 출력 `EVM.all.gff3` 기준, 8 런 전수):

    런       최대exon  최대aa   기준초과
    Ath          92    5,394       0
    Osa          79    5,384       0
    Mtr          90   18,012       0
    Ptr         107   18,700       0
    Dca          94    7,748       0
    Svi_A        79    5,382       0
    Spe_A       117   18,242       1   <- exon 117 > 110
    Sly         669   42,592      57

**Sly 57 개 외에 Spe_A 에도 1 개가 걸린다**(`evm.model.Spe_ch10.1257`,
117 exon / 18,242 aa). 초판 주석은 "Sly 뿐"이라고 적었으나 오측이었다 —
Spe_A 의 최대 exon 을 107(Ptr 값)로 잘못 옮겼다.

다만 Spe_A 의 PREFILTER 는 이 가드가 생기기 전에 만들어졌고, 그 모델은
**random-forest filter 가 이미 걸러냈다**(`filtered.gff3` 최대 exon 103 /
9,511 aa, 초과 0). 발표되는 산출물에는 남아 있지 않다. 그렇더라도 PASA 에
넘기는 입력이 달라지면 해당 좌위 주변 assembly 가 달라질 수 있으므로
"결과가 동일하다"고까지는 말할 수 없다 — 재실행하지 않았다는 사실과 함께
방법론에 적는다.

### aa 10,000~20,000 잔존분 — 단계마다 수가 다르다

기준 미만이라 제거되지 않은 것들. **발표 수치는 `filtered.gff3` 기준**이므로
그 열을 봐야 한다:

    런       EVM  PREFILTER  filtered(발표)
    Mtr       3       3          2
    Ptr       2       2          0
    Spe_A     2       2          0

Mtr 에 남는 2 개(`MtrG00110180.t1` 18,012 aa, `MtrG00270160.t1` 10,918 aa)는
**인트론 없는 단일 exon**(CDS 54.0 kb / 32.8 kb)이다. Sly 의 다중 exon 연쇄
병합과 기전이 다르므로 "잔존 키메라"가 아니라 **단일 exon 거대 ORF**로
따로 부른다. Ptr 의 2 개는 다중 exon(94·58 exon)이었으나 filter 가 제거했다.

## 사용법

    strip_evm_chimeras.py <in.gff3> <out.gff3> [--max-exon 110] [--max-aa 20000]

제거 대상이 없으면 입력을 그대로 복사한다(하류가 항상 같은 경로를 읽게).
요약은 stderr 로 나가 규칙 로그에 남는다.
"""
import argparse
import re
import shutil
import sys
from collections import defaultdict

ID_RE = re.compile(r"(?:^|;)ID=([^;]+)")
PARENT_RE = re.compile(r"(?:^|;)Parent=([^;]+)")


def _attr(pat, s):
    m = pat.search(s)
    return m.group(1) if m else None


def find_chimeras(path, max_exon, max_aa):
    """(제거할 mRNA 집합, 진단용 dict) 반환."""
    exon_n = defaultdict(int)
    cds_bp = defaultdict(int)
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            if f[2] == "exon":
                p = _attr(PARENT_RE, f[8])
                if p:
                    exon_n[p] += 1
            elif f[2] == "CDS":
                p = _attr(PARENT_RE, f[8])
                if p:
                    cds_bp[p] += int(f[4]) - int(f[3]) + 1
    bad = {
        m for m in set(exon_n) | set(cds_bp)
        if exon_n[m] > max_exon or cds_bp[m] / 3 > max_aa
    }
    return bad, exon_n, cds_bp


def strip(src, dst, bad):
    """bad 에 든 mRNA 와 그 자식을 제거. 자식이 전부 사라지는 gene 도 함께 제거."""
    gene_children = defaultdict(set)
    with open(src) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "mRNA":
                continue
            mid, gid = _attr(ID_RE, f[8]), _attr(PARENT_RE, f[8])
            if mid and gid:
                gene_children[gid].add(mid)
    # 정상 isoform 이 하나라도 남는 gene 은 보존한다
    bad_gene = {g for g, kids in gene_children.items() if kids and kids <= bad}

    kept = dropped = 0
    with open(src) as fin, open(dst, "w") as fout:
        for line in fin:
            if line.startswith("#"):
                fout.write(line)
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                fout.write(line)
                continue
            fid, par = _attr(ID_RE, f[8]), _attr(PARENT_RE, f[8])
            if ((f[2] == "gene" and fid in bad_gene)
                    or (f[2] == "mRNA" and fid in bad)
                    or (f[2] not in ("gene", "mRNA") and par in bad)):
                dropped += 1
            else:
                fout.write(line)
                kept += 1
    return kept, dropped, len(bad_gene)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("src")
    ap.add_argument("dst")
    ap.add_argument("--max-exon", type=int, default=110)
    ap.add_argument("--max-aa", type=int, default=20000)
    a = ap.parse_args()

    bad, exon_n, cds_bp = find_chimeras(a.src, a.max_exon, a.max_aa)
    if not bad:
        shutil.copyfile(a.src, a.dst)
        print(f"[strip_evm_chimeras] 키메라 없음 (기준 exon>{a.max_exon} 또는 "
              f"aa>{a.max_aa}) — 입력 그대로 사용", file=sys.stderr)
        return 0

    kept, dropped, n_gene = strip(a.src, a.dst, bad)
    aas = sorted(cds_bp[m] // 3 for m in bad)
    exs = sorted(exon_n[m] for m in bad)
    print(f"[strip_evm_chimeras] 키메라 {len(bad)} mRNA / {n_gene} gene 제거 "
          f"(기준 exon>{a.max_exon} 또는 aa>{a.max_aa})", file=sys.stderr)
    print(f"[strip_evm_chimeras]   exon {exs[0]}~{exs[-1]}, "
          f"단백질 {aas[0]:,}~{aas[-1]:,} aa", file=sys.stderr)
    print(f"[strip_evm_chimeras]   라인 보존 {kept:,} / 제거 {dropped:,} -> {a.dst}",
          file=sys.stderr)
    # 제거된 ID 는 감사 흔적으로 남긴다
    with open(a.dst + ".removed.txt", "w") as fh:
        for m in sorted(bad):
            fh.write(f"{m}\texon={exon_n[m]}\taa={cds_bp[m] // 3}\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
