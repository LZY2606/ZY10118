package app

object AppResources {
    fun read(path: String): String =
        AppResources::class.java.getResourceAsStream(path)
            ?.bufferedReader(Charsets.UTF_8)?.use { it.readText() }
            ?: error("missing resource $path")
}
