nextflow.enable.dsl = 2


//================================================================================
// Derive file names and location from the params.yaml
//================================================================================

fastq_files_list = file(params.input_fofn)

//================================================================================
// Include sub-workflows and (soft) override workflow-level parameters
//================================================================================

include { FORMAT_CONVERSION } from "./workflows/format_conversion/format_conversion.nf"
include { PREPROCESSING_MAPPING } from "./workflows/preprocessing_mapping/preprocessing_mapping.nf"
include { QUALITY_RECALIBRATION } from "./workflows/quality_recalibration/quality_recalibration.nf"
include { VARIANT_DISCOVERY } from "./workflows/variant_discovery/variant_discovery.nf"


//================================================================================
// Prepare channels
//================================================================================


// Convert input manifest to a channel.
fastq_params_ch = channel.fromPath(fastq_files_list)
        .splitText(keepHeader:true)
        .map { line ->
            cols = line.tokenize('\t')
            [
                    file(cols[2]), // fastq_1
                    file(cols[3]), // fastq_2
                    cols[7], // run_date
                    cols[1], // sample_name
                    // NEW
                    cols[4], // lane
                    cols[5], // library_name
                    cols[8], // platform_name
                    cols[6], // platform_unit
                    cols[0], // readgroup_name
                    cols[9] // sequencing_center

                    // OLD 
                    /*
                    cols[4], // library_name
                    cols[7], // platform_name
                    cols[5], // platform_name
                    cols[0], // readgroup_name
                    cols[8] // sequencing_center
                    */
            ]
        }


//================================================================================
// Main workflow
//================================================================================

workflow {

    FORMAT_CONVERSION({ fastq_params_ch })

    PREPROCESSING_MAPPING(FORMAT_CONVERSION.out)

    QUALITY_RECALIBRATION(PREPROCESSING_MAPPING.out)

    VARIANT_DISCOVERY(QUALITY_RECALIBRATION.out)

    VARIANT_DISCOVERY.out
            .map {
                sampleId, vcfFile -> "${sampleId}\ts3://${vcfFile}"
            }
            .collectFile(
                    name: 'merged_vcfs.tsv', newLine: true, storeDir: "${params.outdir}"
            ) 

}
