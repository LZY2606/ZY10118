#!/usr/bin/env python3
"""Generate deterministic offline test fixtures.

Outputs:
  resources/fixtures/refs/hg19_chr1.fasta(.fai)
  resources/fixtures/refs/hg38_chr1.fasta(.fai)
  resources/fixtures/transcripts/transcripts.tsv
  resources/fixtures/chains/hg19_to_hg38.chain.tsv
  resources/fixtures/sample_batch.vcf

All coordinates are 1-based closed. Transcript exon field is comma-separated
"start-end" pairs; CDS is a single "start-end" range.
"""
import os, random

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.join(HERE, "..", "src", "main", "resources", "fixtures")
REF = os.path.join(BASE, "refs")
TX = os.path.join(BASE, "transcripts")
CH = os.path.join(BASE, "chains")
for d in (BASE, REF, TX, CH):
    os.makedirs(d, exist_ok=True)

def seeded_seq(n, seed):
    rnd = random.Random(seed)
    b = "ACGT"
    return "".join(b[rnd.randrange(4)] for _ in range(n))

def write_fasta(path, name, seq, line=70):
    with open(path, "w") as f:
        f.write(">" + name + "\n")
        for i in range(0, len(seq), line):
            f.write(seq[i:i+line] + "\n")
    with open(path + ".fai", "w") as f:
        f.write("%s\t%d\t%d\t%d\t%d\n" % (name, len(seq), len(name) + 2, line, line + 1))

# ---------------------------------------------------------------- hg19
N19 = 3000
s = list(seeded_seq(N19, 3701))
def ovw(pos1, motif):
    for i, ch in enumerate(motif):
        s[pos1 - 1 + i] = ch

# NM_DEMO_P (+): exons 201-289 (89), 370-459 (90); intron 290-369 (80)
# CDS 210-456 = 80 + 87 = 167 -> not divisible by 3; choose 211..456 = 168
# exon1 CDS 211-289 (79)... use 211-288? simpler: define CDS 210-456 length 168?
# 456-210+1 = 247. Let's set CDS 213..456 => 289-213+1=77 + 456-370+1=87 =>164 no.
# Solve: a..289 length L1, 370..b length L2, L1+L2=168, stop inside exon2.
# pick L1=81 => a=209; L2=87 => b=456. stop codon = 454..456.
P_EX1, P_EX2 = (201, 289), (370, 459)
P_CDS_A, P_CDS_B = 209, 456
# Build a CDS of length 168: start ATG, 54 interior codons, stop TAA.
codons = ["ATG"] + (["GCT"] * 40 + ["TTC"] * 14) + ["TAA"]
assert len(codons) == 56
cds = "".join(codons)
assert len(cds) == 168
ovw(P_CDS_A, cds)  # intron splits it: 209..289 = 81 bp then 370..456 = 87 bp

# Variant landmarks (positions chosen inside/around NM_DEMO_P)
# SNV inside exon1 (missense): base at 220 -> change
# insertion in exon1 homopolymer: position 230 area
# deletion across exon boundary: 285..374 removes 285-289 + 370-374 (10bp)
ovw(2001, "CAAAAAAAC")            # poly-A inside seg2-only 2000-2009 (FLIPPED)
# SNV at 2060 is written directly in sample VCF (overlap zone -> MANY_TO_ONE)
ovw(900, "CATCATCATC")            # tandem deletion demo (seg1 UNIQUE)
ovw(950, "ATG" + "GGGG")          # composite multiallelic locus at 950 (seg1)
ovw(2500, "ACGTAC")               # unmapped (no chain segment covers 2500)

# NM_DEMO_N (-): exons 1501-1580, 1661-1740; CDS 1510-1731
# length = 1580-1510+1=71 + 1731-1661+1=71 =>142 not div3. use 1512..1731:
# 1580-1512+1=69, 1731-1661+1=71 => 140. pick CDS 1510..1732 =>71+72=143.
# Solve 69+72=141 => a=1512, b=1732.
N_EX1, N_EX2 = (1501, 1580), (1661, 1740)
N_CDS_A, N_CDS_B = 1512, 1732
n_codons = ["ATG"] + ["CTA"] * 45 + ["TAG"]
assert len("".join(n_codons)) == 141
ovw(N_CDS_A, "".join(n_codons))

s19 = "".join(s)
write_fasta(os.path.join(REF, "hg19_chr1.fasta"), "chr1", s19)

# ---------------------------------------------------------------- hg38
N38 = 9000
t = list(seeded_seq(N38, 3802))
# chain segments hg19 -> hg38:
#  seg1: 100-1999  -> 1000-2899  (+)
#  seg2: 2000-2099 -> 5000-5099  (-) (flipped)
#  seg3: 2100-2199 -> 1000-1099  (+) (many-to-one with seg1 target 1000-2899)
#  seg4: 2200-2299 -> 6000-6079  (+) (short, target gap 6080+ unmapped => INTERRUPTED edge)
#  seg5: 2300-2499 -> 7000-7199  (+)
def copyseg(src_a, dst_a, length, flip=False):
    block = s19[src_a-1:src_a-1+length]
    if flip:
        comp = str.maketrans("ACGT", "TGCA")
        block = block.translate(comp)[::-1]
    for i, ch in enumerate(block):
        t[dst_a - 1 + i] = ch
copyseg(100, 1000, 1900, False)
copyseg(2000, 4000, 100, True)
copyseg(2010, 5050, 90, True)
copyseg(2200, 6000, 41, False)
copyseg(2250, 7000, 250, False)
s38 = "".join(t)
write_fasta(os.path.join(REF, "hg38_chr1.fasta"), "chr1", s38)

# chain TSV: build  fromBuild toBuild chrom srcStart srcEnd strand dstChrom dstStart dstEnd
with open(os.path.join(CH, "hg19_to_hg38.chain.tsv"), "w") as f:
    f.write("# fromBuild=hg19 toBuild=hg38 segments=5\n")
    for row in [
        ("hg19","hg38","chr1",100,1999,"+","chr1",1000,2899),
        ("hg19","hg38","chr1",2000,2099,"-","chr1",4000,4099),  # seg2 exclusive dst -> FLIPPED
        ("hg19","hg38","chr1",2010,2099,"-","chr1",5050,5139),  # seg3 overlaps seg2 source -> MANY_TO_ONE
        ("hg19","hg38","chr1",2200,2240,"+","chr1",6000,6040),  # seg4, gap 2241-2249 -> INTERRUPTED edge
        ("hg19","hg38","chr1",2250,2499,"+","chr1",7000,7249),
    ]:
        f.write("\t".join(str(x) for x in row) + "\n")

# transcripts TSV
with open(os.path.join(TX, "transcripts.tsv"), "w") as f:
    f.write("# transcriptSetVersion=tx-fixture-1\n")
    f.write("txId\tgene\tchrom\tstrand\texons\tcdsStart\tcdsEnd\tbuild\n")
    f.write("NM_DEMO_P\tDEMO1\tchr1\t+\t%d-%d,%d-%d\t%d\t%d\thg19\n"
            % (P_EX1[0], P_EX1[1], P_EX2[0], P_EX2[1], P_CDS_A, P_CDS_B))
    f.write("NM_DEMO_N\tDEMO2\tchr1\t-\t%d-%d,%d-%d\t%d\t%d\thg19\n"
            % (N_EX1[0], N_EX1[1], N_EX2[0], N_EX2[1], N_CDS_A, N_CDS_B))

# sample VCF (hg19). REF bases are read from the generated reference.
def ref(pos, length=1):
    return s19[pos-1:pos-1+length]

p_snv = 220
alt_snv = {"A":"T","T":"A","C":"G","G":"A"}[ref(p_snv)]
p_ins = 2002        # left-aligns to 2000; seg2-only -> FLIPPED; CAAAAAAAC
p_rep = 2060        # SNV inside MANY_TO_ONE overlap zone 2010-2099
p_del = 900          # CATCATCATC : delete CAT -> 900 CA? choose 900 "CA" ALT "C" (trims later)
p_mix = 950          # ATGGGG: biallelic composite-style line
p_bad = 2500
p_exn = 1550
p_cross = 285  # exactly spans exon1 end (289) into the intron
p_int = 2238  # delete 2239-2241: REF 2238-2241 crosses seg4 end 2240 into gap -> INTERRUPTED
lines = [
  "##fileformat=VCFv4.2",
  "##source=fixture",
  "##reference=hg19",
  "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
  "chr1\t%d\t.\t%s\t%s\t.\tPASS\tSRC=labA;NOTE=exonic_snv" % (p_snv, ref(p_snv), alt_snv),
  "chr1\t%d\t.\t%s\t%s\t.\tPASS\tSRC=labB;NOTE=homopolymer_insert_left_align" % (p_ins, ref(p_ins,2), ref(p_ins) + "A" + ref(p_ins+1)),
  "chr1\t%d\t.\t%s\t%s\t.\tPASS\tSRC=labA;NOTE=many_to_one_overlap_snv" % (p_rep, ref(p_rep,1), {"A":"T","T":"A","C":"G","G":"A"}[ref(p_rep,1)]),
  "chr1\t%d\t.\t%s\t%s\t.\tPASS\tSRC=labC;NOTE=tandem_deletion" % (p_del, ref(p_del,4), ref(p_del,1)),
  "chr1\t%d\t.\t%s\t%s\t.\tPASS\tSRC=labA;NOTE=composite_multiallelic" % (p_mix, ref(p_mix,3), ref(p_mix,1) + "," + ref(p_mix,2)),
  "chr1\t%d\t.\t%s\t%s\t.\tPASS\tSRC=labB;NOTE=ref_mismatch" % (p_bad, "T" if ref(p_bad) != "T" else "G", "A"),
  "chr1\t%d\t.\t%s\t%s\t.\tPASS\tSRC=labC;NOTE=negative_strand_exon" % (p_exn, ref(p_exn), {"A":"T","T":"A","C":"G","G":"A"}[ref(p_exn)]),
  "chr1\t%d\t.\t%s\t%s\t.\tPASS\tSRC=labA;NOTE=cross_exon_boundary_deletion" % (p_cross, ref(p_cross,10), ref(p_cross,5)),
  "chr1\t%d\t.\t%s\t%s\t.\tPASS\tSRC=labC;NOTE=liftover_interrupted_at_segment_end" % (p_int, ref(p_int,4), ref(p_int,1)),
  "chr1\t2500\t.\t%s\t%s\t.\tPASS\tSRC=labA;NOTE=unmapped_no_segment" % (ref(2500,1), {"A":"T","T":"A","C":"G","G":"A"}[ref(2500,1)]),
]
with open(os.path.join(BASE, "sample_batch.vcf"), "w") as f:
    f.write("\n".join(lines) + "\n")

print("fixtures generated; hg19 len=%d hg38 len=%d" % (len(s19), len(s38)))
