package app

/**
 * VCF normalization, following the standard bcftools/rtg semantics:
 *
 *  - POS is 1-based; the first base of REF and ALT is the anchor base.
 *    A pure insertion of n bases is written length(REF)=1, length(ALT)=n+1.
 *    A pure deletion of n bases is written length(REF)=n+1, length(ALT)=1.
 *  - Minimal representation strips common suffix (then prefix) bases while
 *    keeping at least one base on each allele.
 *  - Left alignment repeatedly borrows the reference base immediately upstream
 *    and drops one common suffix base (VCF convention: variants are anchored to
 *    the leftmost position they can occupy).
 *  - In a simple repeat the same physical change has several equally valid VCF
 *    placements. We enumerate ALL of them and keep the leftmost-minimal one as
 *    the canonical/preferred placement; the ambiguity is never collapsed away.
 */
object Normalizer {

    const val PREFERENCE_RULE =
        "leftmost placement after minimal representation (VCF left-alignment tie-break); all equivalent placements retained"

    fun normalize(input: Variant, ref: RefGenome): NormalizedAllele {
        val steps = mutableListOf<NormStep>()
        val upper = input.copy(ref = input.ref.uppercase(), alt = input.alt.uppercase())
        steps += NormStep("parse", "accepted VCF 1-based record", upper.pos, upper.ref, upper.alt)

        // 1) REF validation against the active reference build.
        val expected = if (ref.contains(upper.chrom, upper.pos, upper.ref.length))
            ref.substring(upper.chrom, upper.pos, upper.ref.length) else null
        val matches = expected != null && expected == upper.ref
        if (expected == null) {
            steps += NormStep("ref_check", "variant lies outside reference contig bounds", upper.pos)
            return NormalizedAllele(upper, false, null, steps, emptyList(), null, PREFERENCE_RULE)
        }
        steps += NormStep(
            "ref_check",
            if (matches) "REF matches reference build ${ref.build} at ${upper.chrom}:${upper.pos}"
            else "REF mismatch: reference build ${ref.build} has $expected",
            upper.pos, upper.ref, upper.alt,
        )
        if (!matches) {
            return NormalizedAllele(upper, false, expected, steps, emptyList(), null, PREFERENCE_RULE)
        }

        // 2) minimal representation: strip shared suffix, then shared prefix.
        var p = upper
        p = stripCommonSuffix(p, steps)
        p = stripCommonPrefix(p, steps)

        // 3) enumerate every equivalent placement across the simple repeat.
        val placements = enumeratePlacements(p, ref, steps)
        val canonical = placements.minWith(compareBy<Placement> { it.pos }.thenBy { it.ref.length })
        steps += NormStep("canonical", "preferred placement: leftmost minimal among ${placements.size} equivalent positions",
            canonical.pos, canonical.ref, canonical.alt)
        return NormalizedAllele(upper, true, expected, steps, placements, canonical, PREFERENCE_RULE)
    }

    fun stripCommonSuffix(v: Variant, steps: MutableList<NormStep>): Variant {
        var cur = v
        while (cur.ref.length > 1 && cur.alt.length > 1 &&
            cur.ref.last() == cur.alt.last()) {
            cur = Variant(cur.chrom, cur.pos, cur.ref.dropLast(1), cur.alt.dropLast(1))
        }
        if (cur != v) steps += NormStep("trim_suffix", "removed shared trailing anchor base", cur.pos, cur.ref, cur.alt)
        return cur
    }

    fun stripCommonPrefix(v: Variant, steps: MutableList<NormStep>): Variant {
        var cur = v
        while (cur.ref.length > 1 && cur.alt.length > 1 &&
            cur.ref.first() == cur.alt.first()) {
            cur = Variant(cur.chrom, cur.pos + 1, cur.ref.drop(1), cur.alt.drop(1))
        }
        if (cur != v) steps += NormStep("trim_prefix", "removed shared leading anchor base; POS advanced", cur.pos, cur.ref, cur.alt)
        return cur
    }

    /** One VCF left-alignment step: borrow the upstream base, drop the shared suffix. */
    fun slideLeft(v: Variant, ref: RefGenome): Variant? {
        if (v.pos <= 1) return null
        if (v.ref.last() != v.alt.last()) return null
        val up = ref.base(v.chrom, v.pos - 1)
        return Variant(v.chrom, v.pos - 1, (up + v.ref).dropLast(1), (up + v.alt).dropLast(1))
    }

    /** Mirror step: borrow the downstream base, drop the shared prefix. */
    fun slideRight(v: Variant, ref: RefGenome): Variant? {
        val after = v.pos + v.ref.length
        if (!ref.contains(v.chrom, after, 1)) return null
        if (v.ref.first() != v.alt.first()) return null
        val down = ref.base(v.chrom, after)
        return Variant(v.chrom, v.pos + 1, (v.ref + down).drop(1), (v.alt + down).drop(1))
    }

    /**
     * Enumerate every distinct minimal placement reachable by sliding within a
     * simple repeat. Left-alignment terminates by itself when the suffix anchors
     * stop matching (the repeat boundary); the right walk stops at the mirrored
     * boundary. The leftmost placement is VCF-canonical and preferred, while
     * every other physical position is retained rather than collapsed away.
     */
    fun enumeratePlacements(minimal: Variant, ref: RefGenome, steps: MutableList<NormStep>): List<Placement> {
        var leftmost = minimal
        while (true) {
            leftmost = slideLeft(leftmost, ref) ?: break
        }
        val seen = LinkedHashMap<Triple<Int, String, String>, Placement>()
        seen[Triple(leftmost.pos, leftmost.ref, leftmost.alt)] =
            Placement(leftmost.chrom, leftmost.pos, leftmost.ref, leftmost.alt, 0)
        var at = leftmost
        var slide = 0
        while (true) {
            val nxt = slideRight(at, ref) ?: break
            slide += 1
            val k = Triple(nxt.pos, nxt.ref, nxt.alt)
            if (k in seen) break
            seen[k] = Placement(nxt.chrom, nxt.pos, nxt.ref, nxt.alt, slide)
            at = nxt
        }
        if (seen.size > 1) {
            steps += NormStep("repeat_ambiguity",
                "simple repeat yields ${seen.size} equivalent minimal placements; all retained",
                leftmost.pos, leftmost.ref, leftmost.alt)
        }
        return seen.values.toList()
    }
}
