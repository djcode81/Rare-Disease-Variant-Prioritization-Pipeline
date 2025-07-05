#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Import modules
include { RARE_DISEASE_PRIORITIZE } from '../modules/prioritization/rare_disease_filter'
include { INHERITANCE_ANALYSIS }    from '../modules/prioritization/inheritance_pattern'
include { CLINICAL_REPORT }         from '../modules/prioritization/clinical_report'

workflow RARE_DISEASE_ENHANCED {
    
    take:
    vcf_channel      // channel: [ val(meta), path(vcf) ]
    pedigree_file    // path: pedigree file (optional)
    gene_panels      // path: gene panels file (optional)
    config_file      // path: prioritization config file (optional)
    
    main:
    
    ch_versions = Channel.empty()
    
    //
    // Step 1: Prioritize variants for rare disease analysis
    //
    RARE_DISEASE_PRIORITIZE (
        vcf_channel,
        config_file
    )
    ch_versions = ch_versions.mix(RARE_DISEASE_PRIORITIZE.out.versions)
    
    //
    // Step 2: Analyze inheritance patterns
    //
    INHERITANCE_ANALYSIS (
        vcf_channel,
        pedigree_file,
        gene_panels
    )
    ch_versions = ch_versions.mix(INHERITANCE_ANALYSIS.out.versions)
    
    //
    // Step 3: Generate clinical report
    //
    CLINICAL_REPORT (
        RARE_DISEASE_PRIORITIZE.out.prioritized_variants,
        RARE_DISEASE_PRIORITIZE.out.summary,
        pedigree_file
    )
    ch_versions = ch_versions.mix(CLINICAL_REPORT.out.versions)
    
    emit:
    prioritized_variants    = RARE_DISEASE_PRIORITIZE.out.prioritized_variants
    inheritance_patterns    = INHERITANCE_ANALYSIS.out.inheritance_patterns
    compound_heterozygous   = INHERITANCE_ANALYSIS.out.compound_heterozygous
    de_novo_variants        = INHERITANCE_ANALYSIS.out.de_novo_variants
    clinical_html_report    = CLINICAL_REPORT.out.html_report
    clinical_pdf_report     = CLINICAL_REPORT.out.pdf_report
    summary_stats           = RARE_DISEASE_PRIORITIZE.out.summary
    inheritance_summary     = INHERITANCE_ANALYSIS.out.summary
    versions               = ch_versions
}

workflow {
    
    //
    // Input validation and channel creation
    //
    if (!params.input) {
        error "Please provide an input VCF file with --input"
    }
    
    // Create meta map
    def meta = [
        id: params.sample_id ?: file(params.input).getBaseName(),
        family_id: params.family_id ?: null,
        affected: params.affected ?: false
    ]
    
    // Create input channels
    ch_vcf = Channel.of([meta, file(params.input)])
    
    // Optional files
    ch_pedigree = params.pedigree ? file(params.pedigree) : []
    ch_gene_panels = params.gene_panels ? file(params.gene_panels) : []
    ch_config = params.config ? file(params.config) : []
    
    //
    // Run enhanced workflow
    //
    RARE_DISEASE_ENHANCED (
        ch_vcf,
        ch_pedigree,
        ch_gene_panels,
        ch_config
    )
    
    //
    // Print results summary
    //
    RARE_DISEASE_ENHANCED.out.summary_stats.view { sample_meta, summary_file ->
        log.info "Sample ${sample_meta.id}: Variant prioritization complete"
        log.info "Summary: ${summary_file}"
    }
    
    RARE_DISEASE_ENHANCED.out.prioritized_variants.view { sample_meta, variants_file ->
        log.info "Sample ${sample_meta.id}: Prioritized variants saved"
        log.info "Variants: ${variants_file}"
    }
    
    RARE_DISEASE_ENHANCED.out.inheritance_summary.view { sample_meta, inheritance_file ->
        log.info "Sample ${sample_meta.id}: Inheritance analysis complete"
        log.info "Results: ${inheritance_file}"
    }
    
    RARE_DISEASE_ENHANCED.out.clinical_html_report.view { sample_meta, report_file ->
        log.info "Sample ${sample_meta.id}: Clinical report generated"
        log.info "Report: ${report_file}"
    }
}
