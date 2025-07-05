#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RARE DISEASE VARIANT PRIORITIZATION DEMO
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Minimal working example demonstrating rare disease variant prioritization
    following GREGoR consortium guidelines.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Import modules
include { RARE_DISEASE_PRIORITIZE } from '../modules/prioritization/rare_disease_filter'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOW FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow RARE_DISEASE_DEMO {
    
    take:
    vcf_channel    // channel: [ val(meta), path(vcf) ]
    config_file    // path: prioritization config file (optional)
    
    main:
    
    ch_versions = Channel.empty()
    
    //
    // Prioritize variants for rare disease analysis
    //
    RARE_DISEASE_PRIORITIZE (
        vcf_channel,
        config_file
    )
    ch_versions = ch_versions.mix(RARE_DISEASE_PRIORITIZE.out.versions)
    
    emit:
    prioritized_variants = RARE_DISEASE_PRIORITIZE.out.prioritized_variants
    summary             = RARE_DISEASE_PRIORITIZE.out.summary
    versions            = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {
    
    //
    // Input validation and channel creation
    //
    if (!params.input) {
        error "Please provide an input VCF file with --input"
    }
    
    // Create meta map
    def meta = [
        id: params.sample_id ?: file(params.input).getBaseName()
    ]
    
    // Create input channel
    ch_vcf = Channel.of([meta, file(params.input)])
    
    // Config file (optional)
    ch_config = params.config ? file(params.config) : []
    
    //
    // Run workflow
    //
    RARE_DISEASE_DEMO (
        ch_vcf,
        ch_config
    )
    
    //
    // Print summary
    //
    RARE_DISEASE_DEMO.out.summary.view { it ->
        log.info "Sample ${it[0].id}: Prioritization complete. Summary saved to ${it[1]}"
    }
}
