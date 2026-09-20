package app.service

class ConflictException(message: String) : RuntimeException(message)
class BadRequestException(message: String) : RuntimeException(message)
class NotFoundException(message: String) : RuntimeException(message)

data class BatchReceipt(
    val batchId: String,
    val build: String,
    val recordCount: Int,
    val acceptedAlleles: Int,
    val refMismatches: Int,
    val groupsCreated: Int,
    val duplicate: Boolean,
    val receivedAt: String,
) {
    fun toMap(): Map<String, Any?> = linkedMapOf(
        "batchId" to batchId, "build" to build, "recordCount" to recordCount,
        "acceptedAlleles" to acceptedAlleles, "refMismatches" to refMismatches,
        "groupsCreated" to groupsCreated, "duplicate" to duplicate,
        "receivedAt" to receivedAt,
    )
}

data class ImpactItem(
    val canonicalId: String,
    val status: String,
    val destContig: String?,
    val destPos: Int?,
    val segments: List<Int>,
    val note: String,
)

data class ImpactReport(
    val sourceBuild: String,
    val destBuild: String,
    val total: Int,
    val items: List<ImpactItem>,
    val summary: Map<String, Int>,
)
