package app.service

/** Minimal dependency-free JSON encoder for the values this service emits. */
object Json {
    fun stringify(value: Any?): String = when (value) {
        null -> "null"
        is Boolean -> value.toString()
        is Number -> value.toString()
        is String -> encodeString(value)
        is Map<*, *> -> value.entries.joinToString(",", "{", "}") { (k, v) ->
            encodeString(k.toString()) + ":" + stringify(v)
        }
        is Iterable<*> -> value.joinToString(",", "[", "]") { stringify(it) }
        is Array<*> -> value.joinToString(",", "[", "]") { stringify(it) }
        else -> {
            val dc = encodeDataClass(value)
            dc ?: encodeString(value.toString())
        }
    }

    private fun encodeDataClass(value: Any): String? {
        val klass = value::class.java
        if (!klass.isAnnotationPresent(kotlin.Metadata::class.java)) return null
        // Kotlin data classes expose componentN operators; pull property names
        // from the declared constructor in declaration order.
        val ctor = klass.declaredConstructors.firstOrNull() ?: return null
        val names = ctor.parameters.map { it.name }
        val parts = names.mapIndexed { idx, name ->
            val fn = klass.methods.firstOrNull { it.name == "component${idx + 1}" } ?: return null
            fn.isAccessible = true
            encodeString(name) + ":" + stringify(fn.invoke(value))
        }
        return parts.joinToString(",", "{", "}")
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
}

/** Tiny JSON value parser sufficient for review request bodies. */
class JsonParser(private val s: String) {
    private var i = 0

    fun parse(): Any? {
        skipWs()
        val v = readValue()
        skipWs()
        require(i >= s.length) { "trailing JSON content" }
        return v
    }

    private fun skipWs() {
        while (i < s.length && s[i].isWhitespace()) i++
    }

    private fun readValue(): Any? {
        skipWs()
        return when (s[i]) {
            '{' -> readObject()
            '[' -> readArray()
            '"' -> readString()
            't', 'f' -> readBoolean()
            'n' -> { require(s.startsWith("null", i)); i += 4; null }
            else -> readNumber()
        }
    }

    private fun readObject(): Map<String, Any?> {
        val m = LinkedHashMap<String, Any?>()
        i++
        skipWs()
        if (s[i] == '}') { i++; return m }
        while (true) {
            val key = readString()
            skipWs(); require(s[i] == ':'); i++; skipWs()
            m[key] = readValue()
            skipWs()
            when (s[i]) {
                ',' -> { i++; skipWs() }
                '}' -> { i++; return m }
                else -> error("expected ',' or '}'")
            }
        }
    }

    private fun readArray(): List<Any?> {
        val list = ArrayList<Any?>()
        i++
        skipWs()
        if (s[i] == ']') { i++; return list }
        while (true) {
            list.add(readValue())
            skipWs()
            when (s[i]) {
                ',' -> { i++; skipWs() }
                ']' -> { i++; return list }
                else -> error("expected ',' or ']'")
            }
        }
    }

    private fun readString(): String {
        require(s[i] == '"'); i++
        val sb = StringBuilder()
        while (true) {
            val c = s[i++]
            when (c) {
                '"' -> return sb.toString()
                '\\' -> when (val e = s[i++]) {
                    '"' -> sb.append('"')
                    '\\' -> sb.append('\\')
                    '/' -> sb.append('/')
                    'n' -> sb.append('\n')
                    't' -> sb.append('\t')
                    'r' -> sb.append('\r')
                    'u' -> {
                        sb.append(s.substring(i, i + 4).toInt(16).toChar()); i += 4
                    }
                    else -> error("bad escape $e")
                }
                else -> sb.append(c)
            }
        }
    }

    private fun readBoolean(): Boolean {
        if (s.startsWith("true", i)) { i += 4; return true }
        if (s.startsWith("false", i)) { i += 5; return false }
        error("bad boolean")
    }

    private fun readNumber(): Any {
        val start = i
        if (s[i] == '-') i++
        while (i < s.length && (s[i].isDigit() || s[i] in ".eE+-")) i++
        val token = s.substring(start, i)
        return token.toDoubleOrNull()?.takeIf { token.any { it == '.' || it == 'e' || it == 'E' } }
            ?: token.toLong()
    }

    companion object {
        fun parse(body: String): Any? = JsonParser(body).parse()
    }
}
