#!/usr/bin/env python3
"""Offline fixture generator for the variant normalizer demo.

The application itself has no Python dependency: this script only renders the
fixed FASTA/TSV/VCF fixtures shipped under src/main/resources/fixtures and
prints expected normalization values used to cross-check the Kotlin tests.
"""
import json, os, sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIX = os.path.join(ROOT, "src", "main", "resources", "fixtures")
os.makedirs(FIX, exist_ok=True)

COMP = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}


def revcomp(s):
    return "".join(COMP[b] for b in reversed(s))


def filler(length, seed=0):
    """Deterministic 3-base cycle with no long homopolymers."""
    bases = "ACG"
    return "".join(bases[(i + seed) % 3] for i in range(length))


def write_fasta(path, contigs):
    lines = []
    for name, seq in contigs:
        lines.append(">%s" % name)
        for i in range(0, len(seq), 60):
            lines.append(seq[i:i + 60])
    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")


# ---- b37 chr11: plus-strand gene TX_BETA ------------------------------------
# Layout (1-based, inclusive): exons 21-80 and 121-180; CDS 40-170.
CDS = (
    "ATG"
    "GCT" "TTC" "GTT" "CAA"
    "TCT" "GGG" "ACT" "ACC"
    "ATC" "CTG" "GCA" "AGT"
    "TTA"                                  # codon 14 spans g79,g80,g121
    "GGC" "ATC" "TAT" "CGT"
    "GCC" "GAT" "CGT" "CAG"
    "TGG" "AAG" "TGC"
    "TTA" "TTG" "ACT" "TTG"                # unique 3bp deletion GACT at g160
    "TAA"
)
assert len(CDS) == 90, len(CDS)
UTR5 = "ACGTACGTACGGTACCAATT"        # 20 bases, exon positions 21-40
INTRON = filler(40, seed=2)          # 40 bases, positions 81-120
UTR3 = "CATTGGCATTA"                 # 11 bases, exon positions 170-180
GENE11 = UTR5[:19] + CDS[:41] + INTRON + CDS[41:] + UTR3[:11]
assert len(GENE11) == 160, len(GENE11)
LEFT11 = filler(20, seed=1)
RIGHT11 = filler(120, seed=4)
B37_CHR11 = LEFT11 + GENE11 + RIGHT11
assert len(B37_CHR11) == 300

# ---- b37 chr17: minus-strand gene TX_DELTA ----------------------------------
# Mirror layout on the plus strand: exons 121-180 and 221-280; CDS 131-260.
LEFT17 = filler(120, seed=5)
INTRON17 = filler(40, seed=6)
RIGHT17 = filler(20, seed=7)
STREAM17 = CDS  # minus-strand coding stream mirrors the plus-strand CDS
assert len(STREAM17) == 90
U5_17 = "ACGTACGTACGGTACCAATT"
U3_17 = "CATTGGCATTA"
# Plus-strand layout mirrored from TX_BETA: coding stream is STREAM17.
# exon1 plus bases 121-180 hold revcomp(STREAM tail 49 + 11bp 3' UTR context);
# exon2 plus bases 221-280 hold revcomp(19bp 5' UTR context + STREAM head 41).
EXON1_PLUS = revcomp(U3_17 + STREAM17[41:])            # coding g132-180
EXON2_PLUS = revcomp(STREAM17[:41] + U5_17[:19])       # coding g240-280
assert len(EXON1_PLUS) == 60 and len(EXON2_PLUS) == 60
MINUS_GENE = EXON1_PLUS + INTRON17 + EXON2_PLUS
assert len(MINUS_GENE) == 160
B37_CHR17 = LEFT17 + MINUS_GENE + RIGHT17
assert len(B37_CHR17) == 300

# ---- b37 chr20: tandem repeats and chromosome-start anchor ------------------
START20 = "ACGTACGTAC"                                   # 1-10
POLYA = "AAAAAAAA"                                       # 11-18
MID20 = filler(40, seed=8)                               # 19-58
POLYT = "TTTTTTTT"                                       # 59-66
TAIL20 = filler(64, seed=9)                              # 67-130
B37_CHR20 = START20 + POLYA + MID20 + POLYT + TAIL20
assert len(B37_CHR20) == 130

write_fasta(os.path.join(FIX, "ref_b37.fa"), [
    ("chr11", B37_CHR11), ("chr17", B37_CHR17), ("chr20", B37_CHR20),
])

# ---- b38 references: offset, strand flip, gap + many-to-one ----------------
B38_CHR11 = ("G" * 10) + B37_CHR11                       # +10 prefix
B38_CHR17 = revcomp(B37_CHR17)                           # strand flip
# chr20: first 50 bases align one-to-one at +5, then a 10-base insertion gap;

B38_CHR20 = (
    "C" * 5 + B37_CHR20[:50]
    + "T" * 10
    + B37_CHR20[50:120]
    + "A" * 10
)
assert len(B38_CHR20) == 145
write_fasta(os.path.join(FIX, "ref_b38.fa"), [
    ("chr11", B38_CHR11), ("chr17", B38_CHR17), ("chr20", B38_CHR20),
])
print("FASTA fixtures written to", FIX)


# ---- normalization prototype (mirrors app/norm/Normalizer.kt, vt-style) ----
def _normalize_lists(seq, pos, ref, alt):
    p, r, a = pos, list(ref), list(alt)
    while len(r) > 1 and len(a) > 1 and r[-1] == a[-1]:
        r.pop(); a.pop()
    while len(r) > 1 and len(a) > 1 and r[0] == a[0]:
        r.pop(0); a.pop(0); p += 1
    while True:
        if len(r) > len(a):
            if not (r[0] == r[-1] and p > 1 and seq[p-2] == r[0]):
                break
            r.pop(); a.pop(); r.insert(0, seq[p-2]); a.insert(0, seq[p-2]); p -= 1
        elif len(a) > len(r):
            if not (a[0] == a[-1] and p > 1 and seq[p-2] == a[0]):
                break
            r.pop(); a.pop(); r.insert(0, seq[p-2]); a.insert(0, seq[p-2]); p -= 1
        else:
            break
    while len(r) > 1 and len(a) > 1 and r[-1] == a[-1]:
        r.pop(); a.pop()
    while len(r) > 1 and len(a) > 1 and r[0] == a[0]:
        r.pop(0); a.pop(0); p += 1
    return p, "".join(r), "".join(a)


def normalize(ref_seq, pos, ref, alt):
    return _normalize_lists(ref_seq, pos, ref, alt)


def placements(ref_seq, pos, ref, alt):
    cp, cr, ca = normalize(ref_seq, pos, ref, alt)
    if len(cr) == len(ca):
        return [cp], (cp, cr, ca)
    out = [cp]
    p, r, a = cp, list(cr), list(ca)
    while p - 1 + len(r) < len(ref_seq):
        tail = ref_seq[p - 1 + len(r)]
        if len(r) > len(a):
            r2, a2 = r[1:] + [tail], a[1:] + [tail]
        else:
            r2, a2 = [tail] + r[:-1], [tail] + a[:-1]
        if normalize(ref_seq, p + 1, "".join(r2), "".join(a2)) != (cp, cr, ca):
            break
        r, a, p = r2, a2, p + 1
        out.append(p)
    return sorted(out), (cp, cr, ca)


expected = {}
# poly-A run occupies 1-based positions 11-18; B10=C, B19=G.
cases = {
    "chr20_polya_del_anchor10": (10, "CA", "C"),   # delete first A (pos 11)
    "chr20_polya_del_anchor11": (11, "AA", "A"),   # delete interior A
    "chr20_polya_del_anchor18": (18, "AG", "G"),   # delete last A (pos 18)
    "chr20_polya_ins_anchor14": (14, "A", "AA"),   # insert A inside run
    "chr20_polya_ins_anchor19": (19, "G", "AG"),   # insert A before G19
    "chr20_start_sub": (1, "A", "T"),
}
for name, (pp, rr, aa) in cases.items():
    pls, canon = placements(B37_CHR20, pp, rr, aa)
    expected[name] = {"input": [pp, rr, aa], "canonical": list(canon),
                      "equivalent_anchors": pls}
# poly-T run occupies positions 59-66; B58=G, B67=A.
pls, canon = placements(B37_CHR20, 66, "TA", "A")  # delete T at 66
expected["chr20_polyt_del"] = {"input": [66, "TA", "A"],
                               "canonical": list(canon),
                               "equivalent_anchors": pls}
with open(os.path.join(FIX, "expected_norm.json"), "w") as fh:
    json.dump(expected, fh, indent=2, sort_keys=True)
print(json.dumps(expected, indent=2, sort_keys=True))


# ---- transcript fixtures ----------------------------------------------------
def write_tsv(path, header, rows):
    with open(path, "w") as fh:
        fh.write("\t".join(header) + "\n")
        for row in rows:
            fh.write("\t".join(str(x) for x in row) + "\n")


# b37 plus strand: exons 21-80 / 121-180; CDS 40-170.
write_tsv(os.path.join(FIX, "transcripts_b37.tsv"),
          ["transcript_id", "contig", "strand", "exons_start1", "exons_end1",
           "cds_start1", "cds_end1", "cds_exon_starts", "cds_exon_ends",
           "build"],
          [
              ["TX_BETA", "chr11", "+", "21,121", "80,180", 40, 169, "40,121", "80,169", "b37"],
              ["TX_DELTA", "chr17", "-", "121,221", "180,280", 121, 280, "121,240", "169,280", "b37"],
          ])
# b38 plus strand: +10 offset on chr11.
write_tsv(os.path.join(FIX, "transcripts_b38.tsv"),
          ["transcript_id", "contig", "strand", "exons_start1", "exons_end1",
           "cds_start1", "cds_end1", "cds_exon_starts", "cds_exon_ends",
           "build"],
          [["TX_BETA38", "chr11", "+", "31,131", "90,190", 50, 179, "50,131", "90,179", "b38"]])

# ---- liftover map b37 -> b38 ------------------------------------------------
# columns: sourceContig srcStart1 srcEnd1 destContig destStart1 destEnd1 strand
write_tsv(os.path.join(FIX, "map_b37_to_b38.tsv"),
          ["src_contig", "src_start1", "src_end1", "dst_contig",
           "dst_start1", "dst_end1", "dst_strand"],
          [
              ["chr11", 1, 300, "chr11", 11, 310, "+"],
              ["chr17", 1, 300, "chr17", 1, 300, "-"],
              # primary chr20 chain: 1-50 -> 6-55, then 10-base gap, 51-120 -> 66-135
              ["chr20", 1, 50, "chr20", 6, 55, "+"],
              ["chr20", 51, 120, "chr20", 66, 135, "+"],
              # overlapping second chain on chr20 -> many-to-one collisions
              ["chr20", 11, 30, "chr20", 71, 90, "+"],
          ])
print("TSV fixtures written")


# ---- transcript projection prototype (mirrors Kotlin projector) ------------
class Tx:
    def __init__(self, tid, contig, strand, exon_starts, exon_ends,
                 cds_piece_starts, cds_piece_ends):
        self.tid = tid
        self.contig = contig
        self.strand = strand
        self.exons = list(zip(exon_starts, exon_ends))
        self.coding = list(zip(cds_piece_starts, cds_piece_ends))
        self.cds_s = cds_piece_starts[0]
        self.cds_e = cds_piece_ends[-1]

    def in_exon(self, g):
        return any(s <= g <= e for s, e in self.exons)

    def cds_pos(self, g):
        """1-based CDS coordinate of a coding genomic base, or None."""
        if self.strand == "+":
            if not any(s <= g <= e for s, e in self.coding):
                return None
            return sum(min(g, e) - s + 1 for s, e in self.coding if g >= s)
        if not any(s <= g <= e for s, e in self.coding):
            return None
        return sum(e - max(g, s) + 1 for s, e in reversed(self.coding) if e >= g)


TX_BETA = Tx("TX_BETA", "chr11", "+", [21, 121], [80, 180],
             [40, 121], [80, 169])
TX_DELTA = Tx("TX_DELTA", "chr17", "-", [121, 221], [180, 280],
              [121, 240], [169, 280])


def codon_at_cds(seq, tx, cds_index0):
    """Three coding bases for a 0-based codon index, splicing-aware."""
    out = []
    for k in range(3):
        cds1 = cds_index0 * 3 + k + 1
        if tx.strand == "+":
            g = _genomic_from_cds_plus(tx, cds1)
            out.append(seq[g - 1])
        else:
            g = _genomic_from_cds_minus(tx, cds1)
            out.append(COMP[seq[g - 1]])
    return "".join(out)


def _genomic_from_cds_plus(tx, cds1):
    remaining = cds1
    for s, e in tx.coding:
        length = e - s + 1
        if remaining <= length:
            return s + remaining - 1
        remaining -= length
    raise ValueError(cds1)


def _genomic_from_cds_minus(tx, cds1):
    remaining = cds1
    for s, e in reversed(tx.coding):
        length = e - s + 1
        if remaining <= length:
            return e - remaining + 1
        remaining -= length
    raise ValueError(cds1)


AA = {
    "ATG": "M", "GCT": "A", "TTC": "F", "GTT": "V", "CAA": "Q", "TCT": "S",
    "GGG": "G", "ACT": "T", "ACC": "T", "ATC": "I", "CTG": "L", "GCA": "A",
    "AGT": "S", "TTA": "L", "GGC": "G", "TAT": "Y", "CGT": "R", "GCC": "A",
    "TTG": "L", "GTC": "V", "TAC": "Y", "TGG": "W", "GCT": "A", "TAA": "*",
    "TAG": "*", "TGA": "*",
}


def aa_of(codon):
    return AA.get(codon, "?")


proj_expected = {}
# TX_BETA checks
checks_beta = [
    (40, "ATG", "M"),   # start codon, exon1
    (79, "TTA", "L"),   # codon14 = cds40(g79),cds41(g80),cds42(g121)
    (169, "TTG", "L"),  # codon29 g168,g169,g170
]
for g, codon, aa in checks_beta:
    cds = TX_BETA.cds_pos(g)
    got = codon_at_cds(B37_CHR11, TX_BETA, (cds - 1) // 3)
    proj_expected["beta_g%d" % g] = {"cds": cds, "codon": got, "aa": aa_of(got),
                                     "expect_codon": codon}
# TX_DELTA spot checks
for g in (280, 121):
    cds = TX_DELTA.cds_pos(g)
    got = codon_at_cds(B37_CHR17, TX_DELTA, (cds - 1) // 3)
    proj_expected["delta_g%d" % g] = {"cds": cds, "codon": got, "aa": aa_of(got)}
with open(os.path.join(FIX, "expected_proj.json"), "w") as fh:
    json.dump(proj_expected, fh, indent=2, sort_keys=True)
print(json.dumps(proj_expected, indent=2, sort_keys=True))


# ---- seed VCF records -------------------------------------------------------
def base(seq, pos1, length):
    return seq[pos1 - 1:pos1 - 1 + length]


vcf_rows = [
    # CHROM POS ID REF ALT ...
    # 1) chromosome-start substitution on chr20
    ("chr20", 1, "START_SUB", base(B37_CHR20, 1, 1), "T"),
    # 2) multi-allelic inside poly-A: deletion + insertion reported mid-run;
    #    both ALTs normalize to the same equivalence group as separate child
    #    evidence that must share one group without duplicating the parent.
    ("chr20", 14, "POLY_MULTI", "A", "AA,C"),
    # 3) left-alignment demo poly-T deletion reported far right
    ("chr20", 66, "POLYT_RIGHTMOST", "TA", "A"),
    # 4) reference mismatch: claim REF that does not match b37
    ("chr20", 25, "REF_MISMATCH", "ZZZ", "Z"),
    # 5) start codon hit on TX_BETA: g40 ATG -> ACG (M->T)
    ("chr11", 40, "BETA_START", base(B37_CHR11, 40, 1), "C"),
    # 6) stop codon loss: g167 is first base of terminal TAA -> AAA
    ("chr11", 167, "BETA_STOP", base(B37_CHR11, 167, 1), "A"),
    # 7) cross-exon codon change: g80 is last base of exon1 codon14
    ("chr11", 80, "BETA_EXON_BOUNDARY", base(B37_CHR11, 80, 1), "C"),
    # 8) intronic variant
    ("chr11", 100, "BETA_INTRON", base(B37_CHR11, 100, 1), "T"),
    # 9) 5'UTR base (g30 is UTR)
    ("chr11", 30, "BETA_UTR5", base(B37_CHR11, 30, 1), "G"),
    # 10) intergenic far away
    ("chr11", 5, "INTERGENIC", base(B37_CHR11, 5, 1), "T"),
    # 11) minus-strand start codon hit: g280 is CDS1 (first base of ATG read)
    ("chr17", 280, "DELTA_START", base(B37_CHR17, 280, 1), "C"),
    # 12) chr20 variant inside primary chain (liftover clean)
    ("chr20", 30, "LIFT_CLEAN", base(B37_CHR20, 30, 1), "G"),
    # 13) chr20 deletion whose REF span (g48-g51) crosses the 50|51 chain gap
    ("chr20", 48, "LIFT_GAP", base(B37_CHR20, 48, 4),
     base(B37_CHR20, 48, 1)),
    # 14) chr20 variant inside many-to-one collision region src 11-30
    ("chr20", 20, "LIFT_COLLISION", base(B37_CHR20, 20, 1), "T"),
    # 15) unique in-frame deletion of GAC at g160 (window TTGACTT)
    ("chr11", 160, "BETA_INFRAME_DEL", base(B37_CHR11, 160, 4),
     base(B37_CHR11, 160, 1)),
]
vcf = [
    "##fileformat=VCFv4.2",
    "##source=variant-normalizer-fixture",
    "##reference=b37",
    '##INFO=<ID=NOTE,Number=1,Type=String,Description="fixture case">',
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
]
for chrom, pp, vid, ref, alt in vcf_rows:
    vcf.append("%s\t%d\t%s\t%s\t%s\t.\tPASS\tNOTE=%s" % (chrom, pp, vid, ref, alt, vid))
with open(os.path.join(FIX, "seed_b37.vcf"), "w") as fh:
    fh.write("\n".join(vcf) + "\n")
print("seed VCF written:", len(vcf_rows), "records")
for row in vcf_rows:
    print(row[2], row[0], row[1], row[3], "->", row[4])
