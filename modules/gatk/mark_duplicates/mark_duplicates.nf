nextflow.enable.dsl = 2

params.gatk_path = "gatk"
params.java_opts = "-Xms4000m"
params.compression_level = 5

process GATK_MARK_DUPLICATES {
    tag "${sampleId}"
    label 'gatk4_container'

    input:
    //tuple val(sampleId), tuple val(lanes), tuple path(input_mapped_merged_bams)
    tuple val(sampleId), val(lanes), path(input_mapped_merged_bams)

    output:
    val(sampleId)
    path("${sampleId}_merged.deduped.bam")
    path("${sampleId}_merged.deduped.metrics.txt")

    script:
    def input_mapped_merged_bams2 = input_mapped_merged_bams.collect { "--INPUT $it " }.join(" ")
    

    """
    ${params.gatk_path} --java-options "-Dsamjdk.compression_level=${params.compression_level} ${params.java_opts}" \
                        MarkDuplicates \
                        ${input_mapped_merged_bams2} \
                        --OUTPUT ${sampleId}_merged.deduped.bam \
                        --METRICS_FILE ${sampleId}_merged.deduped.metrics.txt \
                        --VALIDATION_STRINGENCY SILENT \
                        --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
                        --ASSUME_SORT_ORDER "queryname" \
                        --CREATE_MD5_FILE true
    """

    stub:

    """
    touch ${sampleId}_merged.deduped.bam
    touch ${sampleId}_merged.deduped.metrics.txt
    """

}
