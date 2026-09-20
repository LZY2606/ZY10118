package app.norm

import app.model.NormStep
import app.model.Normalized

/**
 * Minimal representation plus parsimonious left alignment (vt-style).
 * Positions are 1-based VCF positions.
 *
 * Left alignment walks the allele through repeat context:
 *   - deletion (REF longer): the window moves left while REF's first and last
 *     bases agree and the reference base immediately left matches them;
 *   - insertion (ALT longer): the analogous condition on ALT.
 *
 * Equivalent anchors:
 *   In a tandem repeat one molecule admits several VCF representations.
 *   [equivalentAnchors] enumerates every physical anchor that re-normalizes to
 *   the canonical record. The preferred representation is the leftmost anchor.
 *   Substitutions/MNVs always have a single physical anchor.
 */
class Normalizer(private val reference: String) {

    private fun trimRight(r: String, a: String): Pair<String, String> {
        val rb = StringBuilder(r); val ab = StringBuilder(a)
        while (rb.length > 1 && ab.length > 1 && rb.last() == ab.last()) {
            rb.deleteCharAt(rb.length - 1); ab.deleteCharAt(ab.length - 1)
        }
        return rb.toString() to ab.toString()
    }

    private fun trimLeft(pos: Int, r: String, a: String): Triple<Int, String, String> {
        var p = pos
        val rb = StringBuilder(r); val ab = StringBuilder(a)
        while (rb.length > 1 && ab.length > 1 && rb.first() == ab.first()) {
            rb.deleteCharAt(0); ab.deleteCharAt(0); p += 1
        }
        return Triple(p, rb.toString(), ab.toString())
    }

    /** Core alignment without audit steps and without anchor enumeration. */
    private fun align(startPos: Int, startRef: String, startAlt: String): Triple<Int, String, String> {
        var (r, a) = trimRight(startRef, startAlt)
        var t = trimLeft(startPos, r, a)
        var p = t.first; r = t.second; a = t.third
        var guard = 0
        while (true) {
            if (r.length > a.length) {
                if (!(r.first() == r.last() && p > 1 && reference[p - 2] == r.first())) break
                val lead = reference[p - 2]
                r = lead + r.dropLast(1)
                a = lead + a.dropLast(1)
                p -= 1
            } else if (a.length > r.length) {
                if (!(a.first() == a.last() && p > 1 && reference[p - 2] == a.first())) break
                val lead = reference[p - 2]
                r = lead + r.dropLast(1)
                a = lead + a.dropLast(1)
                p -= 1
            } else break
            if (++guard > reference.length + 4)
                throw IllegalStateException("left alignment did not terminate at $startPos")
        }
        val tr = trimRight(r, a); r = tr.first; a = tr.second
        val tl = trimLeft(p, r, a)
        return Triple(tl.first, tl.second, tl.third)
    }

    fun normalize(startPos: Int, startRef: String, startAlt: String): Normalized {
        val steps = ArrayList<NormStep>()
        var pos = startPos
        var r = startRef
        var a = startAlt

        val tr0 = trimRight(r, a)
        if (tr0.first != r || tr0.second != a) {
            r = tr0.first; a = tr0.second
            steps.add(NormStep("TRIM_RIGHT", pos, r, a, "shared trailing base removed"))
        }
        val tl0 = trimLeft(pos, r, a)
        if (tl0.second != r || tl0.third != a) {
            pos = tl0.first; r = tl0.second; a = tl0.third
            steps.add(NormStep("TRIM_LEFT", pos, r, a, "shared leading base removed"))
        }

        var guard = 0
        while (true) {
            val canShift: Boolean
            if (r.length > a.length) {
                canShift = r.first() == r.last() && pos > 1 && reference[pos - 2] == r.first()
            } else if (a.length > r.length) {
                canShift = a.first() == a.last() && pos > 1 && reference[pos - 2] == a.first()
            } else break
            if (!canShift) break
            val lead = reference[pos - 2]
            r = lead + r.dropLast(1)
            a = lead + a.dropLast(1)
            pos -= 1
            guard++
            steps.add(NormStep("LEFT_ALIGN", pos, r, a,
                "prepended reference base '$lead'; window shifted left to $pos"))
            if (guard > reference.length + 4)
                throw IllegalStateException("left alignment did not terminate at $startPos")
        }
        val tr1 = trimRight(r, a)
        if (tr1.first != r || tr1.second != a) {
            r = tr1.first; a = tr1.second
            steps.add(NormStep("TRIM_RIGHT", pos, r, a, "trailing base removed after walk"))
        }
        val tl1 = trimLeft(pos, r, a)
        if (tl1.second != r || tl1.third != a) {
            pos = tl1.first; r = tl1.second; a = tl1.third
            steps.add(NormStep("TRIM_LEFT", pos, r, a, "leading base removed after walk"))
        }

        val anchors = enumerateAnchors(pos, r, a)
        return Normalized(pos, r, a, steps, anchors, anchors.size == 1)
    }

    /** Enumerate physical anchors by repeatedly inverting one left shift.
     *  Deletion and insertion use different inverse window transforms. */
    private fun enumerateAnchors(canonPos: Int, canonRef: String, canonAlt: String): List<Int> {
        if (canonRef.length == canonAlt.length) return listOf(canonPos)
        val anchors = arrayListOf(canonPos)
        var pos = canonPos
        var r = canonRef
        var a = canonAlt
        var guard = 0
        while (pos - 1 + r.length < reference.length) {
            val tail = reference[pos - 1 + r.length]
            val r2: String
            val a2: String
            if (r.length > a.length) {        // deletion inverse
                r2 = r.drop(1) + tail
                a2 = a.drop(1) + tail
            } else {                          // insertion inverse
                r2 = tail + r.dropLast(1)
                a2 = tail + a.dropLast(1)
            }
            val aligned = align(pos + 1, r2, a2)
            if (aligned.first != canonPos || aligned.second != canonRef || aligned.third != canonAlt) break
            r = r2; a = a2; pos += 1
            anchors.add(pos)
            if (++guard > reference.length + 4) break
        }
        return anchors.sorted()
    }
}
