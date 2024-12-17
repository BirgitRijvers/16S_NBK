/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FASTQC                    } from '../modules/nf-core/fastqc/main'
include { PORECHOP_PORECHOP         } from '../modules/nf-core/porechop/porechop/main'
include { NANOPLOT                  } from '../modules/nf-core/nanoplot/main'
include { FILTLONG                  } from '../modules/nf-core/filtlong/main'
include { MINIMAP2_INDEX            } from '../modules/nf-core/minimap2/index/main'
include { MINIMAP2_ALIGN            } from '../modules/nf-core/minimap2/align/main'
include { SAMTOOLS_VIEW             } from '../modules/nf-core/samtools/view/main'
include { SAMTOOLS_FASTQ            } from '../modules/nf-core/samtools/fastq/main'
include { KRAKEN2_KRAKEN2           } from '../modules/nf-core/kraken2/kraken2/main'
include { KRAKENTOOLS_KREPORT2KRONA } from '../modules/nf-core/krakentools/kreport2krona/main'
include { KRONA_KTIMPORTTEXT        } from '../modules/nf-core/krona/ktimporttext/main'
include { MULTIQC                   } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap          } from 'plugin/nf-validation'
include { paramsSummaryMultiqc      } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML    } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText    } from '../subworkflows/local/utils_nfcore_nbksixteens_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE PARAMS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
if(params.fasta){
    ch_fasta = Channel.fromPath(params.fasta, checkIfExists: true).collect()
        .map{ it -> [[id:it[0].getSimpleName()], it[0]]}
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow NBKSIXTEENS {

    take:
    // channel: samplesheet read in from --input
    ch_samplesheet

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    //
    // MODULE: Run FastQC
    //
    FASTQC (
        ch_samplesheet
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    //
    // MODULE: Run Porechop
    //
    PORECHOP_PORECHOP (
        ch_samplesheet
    )
    ch_versions = ch_versions.mix(PORECHOP_PORECHOP.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(PORECHOP_PORECHOP.out.log.collect{it[1]})
    
    //
    // MODULE: Run filtlong
    //
    FILTLONG ( 
        PORECHOP_PORECHOP.out.reads.map { meta, reads -> [ meta, [], reads ] }
        // ch_samplesheet.map { meta, reads -> [ meta, [], reads ] }
    )
    ch_versions = ch_versions.mix(FILTLONG.out.versions)
    // throws error: Not a valid path value type: java.util.LinkedHashMap
    // ch_multiqc_files = ch_multiqc_files.mix(FILTLONG.out.log.collect{it[1]})

    //
    // MODULE: Run Nanoplot
    //
    NANOPLOT (
        FILTLONG.out.reads
    )
    ch_versions = ch_versions.mix(NANOPLOT.out.versions)
    // Nanoplot is not compatible with MultiQC

    // Check if Minimap2 index is provided
    if (!params.minimap2_index) {
        //
        // MODULE: Run Minimap2 index
        //
        MINIMAP2_INDEX(
            ch_fasta
        )
        ch_versions = ch_versions.mix(MINIMAP2_INDEX.out.versions)
        ch_minimap2_index = MINIMAP2_INDEX.out.index
    }
    else {
        // Use provided Minimap2 index if provided
        ch_minimap2_index = Channel.value([[id:'input_genome_index'], params.minimap2_index])
    }

    //
    // MODULE: Run Minimap2 align
    //
    MINIMAP2_ALIGN(
        FILTLONG.out.reads,
        ch_minimap2_index,
        true,
        "bai",
        false,
        false
    )
    ch_versions = ch_versions.mix(MINIMAP2_ALIGN.out.versions.first())
    // Create channel with mapping outputs, metadata and empty list for reference
    ch_host_reads = MINIMAP2_ALIGN.out.bam.map {meta, reads ->[meta, reads, []] }

    // TO BE ADDED
    // MODULE: Run Samtools stats to generate stats on mapping
    //


    //
    // MODULE: Run Samtools view to extract unmapped reads
    //
    SAMTOOLS_VIEW (
        ch_host_reads,
        // Empty list for reference fasta and meta
        [[],[]],
        // Empty list for qname
        []
    )
    ch_versions = ch_versions.mix(SAMTOOLS_VIEW.out.versions)

    //
    // MODULE: Run Samtools fastq to convert unmapped reads to FASTQ
    //
    SAMTOOLS_FASTQ (
        SAMTOOLS_VIEW.out.bam, 
        false 
        )
    ch_versions = ch_versions.mix(SAMTOOLS_FASTQ.out.versions)

    // 
    // MODULE: Run Kraken2
    //
    KRAKEN2_KRAKEN2 (
        SAMTOOLS_FASTQ.out.other,
        params.kraken2_db,
        false,
        false
    )
    ch_versions = ch_versions.mix(KRAKEN2_KRAKEN2.out.versions.first())
    // ch_multiqc_files = ch_multiqc_files.mix(KRAKEN2_KRAKEN2.out.report)

    //
    // MODULE: Run KrakenTools kreport2krona
    //
    KRAKENTOOLS_KREPORT2KRONA (
        KRAKEN2_KRAKEN2.out.report
    )
    ch_versions = ch_versions.mix(KRAKENTOOLS_KREPORT2KRONA.out.versions)

    //
    // MODULE: Run Krona ktimporttext
    //
    KRONA_KTIMPORTTEXT (
        KRAKENTOOLS_KREPORT2KRONA.out.txt
    )
    ch_versions = ch_versions.mix(KRONA_KTIMPORTTEXT.out.versions)

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_pipeline_software_mqc_versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))

    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
