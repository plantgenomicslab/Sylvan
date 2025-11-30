
#!/usr/bin/env python
import sys
import argparse
import re

parser = argparse.ArgumentParser(description='Usage:')

parser.add_argument('augustus_gff3', help='augustus.gff3')
parser.add_argument('transfrag_gff3', help='transfrag.gff3')
parser.add_argument('genewise_gff3', help='genewise.gff3')
parser.add_argument('intron_gff', help='intron.gff')
parser.add_argument('-o', '--overlap', type=int, default=30, help='overlap')
parser.add_argument('-m', '--min_augustus_transcriptSupport_percentage', type=float, default=10.0, help='min_augustus_transcriptSupport_percentage')
parser.add_argument('-n', '--min_augustus_intronSupport_number', type=int, default=1, help='min_augustus_intronSupport_number')
parser.add_argument('-r', '--min_augustus_intronSupport_ratio', type=float, default=0.01, help='min_augustus_intronSupport_ratio')
parser.add_argument('--more_strict', action='store_true', help='more_strict')

args = parser.parse_args()

augustus = {}
transfrag = {}
genewise = {}
intron = {}
geneRegion = {}
aug_gene = {}
aug_gene_list = []


def combine_transfrag_and_genewise(gene_id, ID):
    """Combine transfrag and genewise gene models."""
    gff3_out = {}
    gff3_out[gene_id] = {"CDS": {}, "exon": {}, "strand": None, "chr": None}

    if gene_id in transfrag:
        transfrag_info = transfrag[gene_id]
        for line in transfrag_info.split("\n"):
            if not line.strip():
                continue
            fields = line.split("\t")
            if len(fields) < 9:
                continue
            if fields[2] == "CDS":
                gff3_out[gene_id]["CDS"][f"{fields[3]}\t{fields[4]}\t{fields[5]}\t{fields[6]}\t{fields[7]}"] = 1
                gff3_out[gene_id]["exon"][f"{fields[3]}\t{fields[4]}"] = 1
            elif fields[2] == "gene" or fields[2] == "mRNA":
                gff3_out[gene_id]["chr"] = fields[0]
                gff3_out[gene_id]["strand"] = fields[6]

    return gff3_out


def get_intron_support(gene_id):
    augustus_info = augustus[gene_id]
    augustus_info = augustus_info.split("\n")
    total_num, support_num = 0, 0
    for line in augustus_info:
        if "\tintron\t" in line:
            total_num += 1
            data = line.split("\t")
            if f"{data[0]}\t{data[3]}\t{data[4]}" in intron:
                support_num += 1
    return f"{support_num}/{total_num}"


def add_UTR(augustus_id, *id):
    augustus_info = augustus[augustus_id]
    gff3_out, gene_id, strand, augustus_CDS = {}, None, None, []

    match = re.search(r"(\S+?)\t\S+?\tmRNA\t\d+?\t\d+?\t\S+?\t(\S+?)\t\S+?\tID=([^;\s]+)", augustus_info)
    if match:
        gene_id = match.group(3)
        strand = match.group(2)
        gff3_out[gene_id] = {"chr": match.group(1), "strand": strand}

    for line in augustus_info.split("\n"):
        data = line.split("\t")
        if data[2] == "start_codon":
            gff3_out[gene_id]["start_codon"] = {f"{data[3]}\t{data[4]}\t{data[5]}\t{data[6]}\t{data[7]}": 1}
        elif data[2] == "stop_codon":
            gff3_out[gene_id]["stop_codon"] = {f"{data[3]}\t{data[4]}\t{data[5]}\t{data[6]}\t{data[7]}": 1}
        elif data[2] == "CDS":
            gff3_out[gene_id]["CDS"] = {f"{data[3]}\t{data[4]}\t{data[5]}\t{data[6]}\t{data[7]}": 1}
            augustus_CDS.append(f"{data[3]}\t{data[4]}")
    augustus_CDS.sort()

    utr5_cds = augustus_CDS[0].split("\t")
    utr3_cds = augustus_CDS[-1].split("\t")

    augustus_CDS.pop(0)
    augustus_CDS.pop(-1)
    for exon in augustus_CDS:
        data = exon.split("\t")
        gff3_out[gene_id]["exon"] = {f"{data[0]}\t{data[1]}": 1}

    utr5_status, utr3_status = 0, 0
    for transfrag_id in id:
        if transfrag_id in transfrag:
            transfrag_info = transfrag[transfrag_id]
            utr5_ok, utr3_ok = 0, 0
            transfrag_exon = []
            for line in transfrag_info.split("\n"):
                data = line.split("\t")
                if data[2] == "CDS":
                    if int(data[3]) <= int(utr5_cds[0]) and int(data[4]) >= int(utr5_cds[1]):
                        utr5_ok = 1
                    if int(data[3]) <= int(utr3_cds[0]) and int(data[4]) >= int(utr3_cds[1]):
                        utr3_ok = 1
                elif data[2] == "exon":
                    transfrag_exon.append(f"{data[3]}\t{data[4]}")

            transfrag_exon.sort()

            if utr5_ok == 1:
                for exon in transfrag_exon:
                    data = exon.split("\t")
                    if int(data[0]) <= int(utr5_cds[1]):
                        if int(data[1]) > int(utr5_cds[1]):
                            gff3_out[gene_id]["exon"] = {f"{data[0]}\t{utr5_cds[1]}": 1}
                        else:
                            gff3_out[gene_id]["exon"] = {f"{data[0]}\t{data[1]}": 1}
                utr5_status = 1

            if utr3_ok == 1:
                for exon in transfrag_exon:
                    data = exon.split("\t")
                    if int(data[1]) >= int(utr3_cds[0]):
                        if int(data[0]) < int(utr5_cds[0]):
                            gff3_out[gene_id]["exon"] = {f"{utr5_cds[0]}\t{data[1]}": 1}
                        else:
                            gff3_out[gene_id]["exon"] = {f"{data[0]}\t{data[1]}": 1}
                utr3_status = 1

    if utr5_status == 0:
        gff3_out[gene_id]["exon"] = {f"{utr5_cds[0]}\t{utr5_cds[1]}": 1}
    if utr3_status== 0:
        gff3_out[gene_id]["exon"] = {f"{utr3_cds[0]}\t{utr3_cds[1]}": 1}

    print(f"{gff3_out[gene_id]['chr']}\t.\tstart_codon\t{gff3_out[gene_id]['start_codon']}\tID={gene_id}.t1.start_codon;Parent={gene_id}.t1")
    print(f"{gff3_out[gene_id]['chr']}\t.\tstop_codon\t{gff3_out[gene_id]['stop_codon']}\tID={gene_id}.t1.stop_codon;Parent={gene_id}.t1")
    for cds in sorted(gff3_out[gene_id]["CDS"]):
        print(f"{gff3_out[gene_id]['chr']}\t.\tCDS\t{cds}\tID={gene_id}.t1.cds;Parent={gene_id}.t1")
    for exon in sorted(gff3_out[gene_id]["exon"]):
        print(f"{gff3_out[gene_id]['chr']}\t.\texon\t{exon}\t.\t{gff3_out[gene_id]['strand']}\t.\tID={gene_id}.t1.exon;Parent={gene_id}.t1")

    return gff3_out

file_list = [sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]]

for i,file in enumerate(file_list):
    with open(file) as f:
        for line in f:
            if line.startswith("#") or line.isspace():
                continue
            if "gene" in line and "ID=" in line:
                id = re.search("ID=([^;\s]+)", line).group(1)
                _ = line.strip().split("\t")
                key = f"{_[3]}\t{_[4]}"
                geneRegion.setdefault(key, {})[id] = 1
            if i==0:
                if "mRNA" in line and "ID=" in line:
                    id = re.search("ID=([^;\s]+)", line).group(1)
                    aug_gene[id] = {}
                    _ = line.strip().split("\t")
                    key = f"{_[3]}\t{_[4]}"
                    geneRegion.setdefault(key, {})[id] = 1
                augustus[id] = augustus.get(id, "") + line
            elif i==1:
                transfrag[id] = transfrag.get(id, "") + line
            elif i==2:
                genewise[id] = genewise.get(id, "") + line
            elif i==3:
                if "intron" in line:
                    _ = line.strip().split("\t")
                    intron[f"{_[0]}\t{_[3]}\t{_[4]}"] = 1

gene_out = {}
for region in sorted(geneRegion.keys()):
    ID = list(geneRegion[region].keys())
    ID_str = ','.join(ID)

    status_augustus, status_transfrag, status_genewise = (0, 0, 0)
    for id in ID:
        if id in augustus:
            status_augustus = 1
        if id in transfrag:
            status_transfrag = 1
        if id in genewise:
            status_genewise = 1

    if status_augustus == 1:
        for gene_id in ID:
            if gene_id in augustus:
                intron_support_result = get_intron_support(gene_id)
                gff3_out = add_UTR(gene_id, ID)
                gene_out.setdefault(gene_id, {"CDS": {}, "exon": {}})

                for cds in gff3_out[gene_id]["CDS"]:
                    gene_out[gene_id]["CDS"][cds] = 1
                for exon in gff3_out[gene_id]["exon"]:
                    gene_out[gene_id]["exon"][exon] = 1

                gene_out[gene_id]["strand"] = gff3_out[gene_id]["strand"]
                gene_out[gene_id]["chr"] = gff3_out[gene_id]["chr"]
                gene_out[gene_id]["intron_support"] = intron_support_result

                if intron_support_result.startswith("0"):
                    gene_out[gene_id]["pfam"] = 1
    elif status_transfrag == 1:
        for gene_id in ID:
            if gene_id in transfrag:
                gff3_out = combine_transfrag_and_genewise(gene_id, ID)
                gene_out.setdefault(gene_id, {"CDS": {}, "exon": {}})
                for cds in gff3_out[gene_id]["CDS"]:
                    gene_out[gene_id]["CDS"][cds] = 1
                for exon in gff3_out[gene_id]["exon"]:
                    gene_out[gene_id]["exon"][exon] = 1

                gene_out[gene_id]["strand"] = gff3_out[gene_id]["strand"]
                gene_out[gene_id]["chr"] = gff3_out[gene_id]["chr"]
    elif status_genewise == 1:
        for gene_id in ID:
            if gene_id in genewise:
                genewise_info = genewise[gene_id]
                gene_out.setdefault(gene_id, {"CDS": {}, "exon": {}})
                for line in genewise_info.split("\n"):
                    fields = line.split("\t")
                    if fields[2] == "CDS":
                        gene_out[gene_id]["CDS"][f"{fields[3]}\t{fields[4]}\t{fields[5]}\t{fields[6]}\t{fields[7]}"] = 1
                        gene_out[gene_id]["exon"][f"{fields[3]}\t{fields[4]}"] = 1
                    elif fields[2] == "gene":
                        gene_out[gene_id]["chr"] = fields[0]
                        gene_out[gene_id]["strand"] = fields[6]


for gene_id in aug_gene:
    out = {}
    gene_pos = {}
    if_intron_support = None

    first_line = augustus[gene_id].split("\n")[0]
    fields = first_line.split("\t")
    chr, strand, gene_score = fields[0], fields[6], fields[5]
    augustus_transcriptSupport_percentage, augustus_intronSupport_number, augustus_intronSupport_ratio, augustus_intronSupport = None, None, None, None
    if re.search(r"exonHintRatio=([^;\s]+)", fields[8]):
        augustus_transcriptSupport_percentage = re.search(r"exonHintRatio=([^;\s]+)", fields[8]).group(1)
    if re.search(r"intronSupport=(\d+)\/(\d+)", fields[8]):
        match = re.search(r"intronSupport=(\d+)\/(\d+)", fields[8])
        augustus_intronSupport_number = match.group(1)
        augustus_intronSupport = match.group(0)
        if match.group(2) == 0:
            augustus_intronSupport_ratio = 0
        else:
            augustus_intronSupport_ratio = match.group(1) / match.group(2)

    mRNA_number = 0
    for mRNA_id in sorted(aug_gene[gene_id].keys()):
        if_intron_support = 1 if mRNA_id not in gene_out or "pfam" not in gene_out[mRNA_id] else 0
        out = {}
        position = {}
        mRNA_number += 1

        if strand == "+":
            number = 0
            for key in sorted(gene_out[mRNA_id]["CDS"].keys()):
                number += 1
                fields = key.split("\t")
                position[fields[0]] = 1
                gene_pos[fields[0]] = 1
                position[fields[1]] = 1
                gene_pos[fields[1]] = 1
                out[f"{chr}\t.\tCDS\t{key}\tID={mRNA_id}.CDS{number};Parent={mRNA_id};"] = fields[0]

            number = 0
            for key in sorted(gene_out[mRNA_id]["exon"].keys(), reverse=True):
                number += 1
                fields = key.split("\t")
                position[fields[0]] = 1
                gene_pos[fields[0]] = 1
                position[fields[1]] = 1
                gene_pos[fields[1]] = 1
                out[f"{chr}\t.\texon\t{key}\t.\t{strand}\t.\tID={mRNA_id}.exon{number};Parent={mRNA_id};"] = fields[0]

