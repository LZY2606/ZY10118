package app

/** Tiny dependency-free JSON (de)serializer for the small values used by this app. */
object Json {
    fun stringify(v: Any?): String = when (v) {
        null -> "null"
        is Boolean -> if (v) "true" else "false"
        is Number -> v.toString()
        is String -> encodeString(v)
        is Map<*, *> -> v.entries.joinToString(",", "{", "}") { (k, val2) ->
            "${encodeString(k.toString())}:${stringify(val2)}"
        }
        is Iterable<*> -> v.joinToString(",", "[", "]") { stringify(it) }
        is Array<*> -> v.joinToString(",", "[", "]") { stringify(it) }
        else -> encodeString(v.toString())
    }

    private fun encodeString(s: String): String {
        val sb = StringBuilder("\"")
        for (c in s) when (c) {
            '"' -> sb.append("\\\"")
            '\\' -> sb.append("\\\\")
            '\n' -> sb.append("\\n")
            '\r' -> sb.append("\\r")
            '\t' -> sb.append("\\t")
            else -> if (c.code < 0x20) sb.append("\\u%04x".format(c.code)) else sb.append(c)
        }
        return sb.append('"').toString()
    }

    @Suppress("UNCHECKED_CAST")
    fun parse(text: String): Any? {
        val p = Parser(text)
        val v = p.readValue()
        p.ws()
        require(p.atEnd()) { "trailing characters in JSON" }
        return v
    }

    fun parseObject(text: String): Map<String, Any?> {
        val v = parse(text)
        require(v is Map<*, *>) { "expected JSON object" }
        return v as Map<String, Any?>
    }

    private class Parser(val s: String) {
        var i = 0
        fun atEnd() = i >= s.length
        fun ws() { while (i < s.length && s[i].isWhitespace()) i++ }
        private fun expect(c: Char) {
            ws(); require(i < s.length && s[i] == c) { "expected '$c' at $i" }; i++
        }
        fun readValue(): Any? {
            ws()
            require(i < s.length) { "unexpected end of JSON" }
            return when (s[i]) {
                '{' -> readObject()
                '[' -> readArray()
                '"' -> readString()
                't', 'f' -> readBool()
                'n' -> { require(s.startsWith("null", i)); i += 4; null }
                else -> readNumber()
            }
        }
        private fun readObject(): Map<String, Any?> {
            expect('{'); ws()
            val m = LinkedHashMap<String, Any?>()
            if (i < s.length && s[i] == '}') { i++; return m }
            while (true) {
                ws(); val k = readString(); ws(); expect(':'); val v = readValue(); m[k] = v
                ws()
                when (s[i]) { ',' -> { i++; continue }
                    '}' -> { i++; break }
                    else -> error("expected ',' or '}' at $i")
                }
            }
            return m
        }
        private fun readArray(): List<Any?> {
            expect('['); ws()
            val l = mutableListOf<Any?>()
            if (i < s.length && s[i] == ']') { i++; return l }
            while (true) {
                l += readValue(); ws()
                when (s[i]) {
                    ',' -> { i++; continue }
                    ']' -> { i++; break }
                    else -> error("expected ',' or ']' at $i")
                }
            }
            return l
        }
        private fun readString(): String {
            expect('"')
            val sb = StringBuilder()
            while (true) {
                val c = s[i++]
                when (c) {
                    '"' -> break
                    '\\' -> {
                        val e = s[i++]
                        when (e) {
                            '"' -> sb.append('"')
                            '\\' -> sb.append('\\')
                            '/' -> sb.append('/')
                            'n' -> sb.append('\n')
                            't' -> sb.append('\t')
                            'r' -> sb.append('\r')
                            'b' -> sb.append('\b')
                            'f' -> sb.append('\u000C')
                            'u' -> { sb.append(s.substring(i, i + 4).toInt(16).toChar()); i += 4 }
                            else -> error("bad escape \\$e")
                        }
                    }
                    else -> sb.append(c)
                }
            }
            return sb.toString()
        }
        private fun readBool(): Boolean {
            if (s.startsWith("true", i)) { i += 4; return true }
            if (s.startsWith("false", i)) { i += 5; return false }
            error("invalid literal at $i")
        }
        private fun readNumber(): Number {
            val start = i
            if (s[i] == '-') i++
            while (i < s.length && (s[i].isDigit() || s[i] in ".eE+-")) i++
            val t = s.substring(start, i)
            return if (t.any { it in ".eE" }) t.toDouble() else t.toLong()
        }
    }
}
