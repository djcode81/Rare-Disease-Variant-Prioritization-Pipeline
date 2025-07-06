# Rare Disease Variant Prioritization Pipeline
NextFlow License: MIT

A production-ready NextFlow pipeline for rare disease variant analysis with dual-parser capability, clinical reporting, and inheritance pattern detection.

## Overview
This pipeline processes VCF files to identify and prioritize variants potentially associated with rare diseases. It features automatic VCF type detection and uses appropriate parsing strategies for different data sources, from clinical databases to patient genomic data.

## Key Features

### Dual-Parser Architecture
- **Automatic VCF type detection** - Routes to appropriate parser based on file content
- **ClinVar parser** - Processes expert-curated clinical significance data
- **VEP parser** - Handles computationally annotated patient genomic data
- **Population parser** - Manages allele frequency data from large cohorts
- **Graceful fallback** - Basic VCF handling for minimal annotation datasets


## Data Sources and Results

### ClinVar Clinical Database
- **Data Source:** NCBI ClinVar with expert-curated pathogenicity assessments
- **Results Location:** `results/real_clinvar/` 
- **Scoring Strategy:** Clinical significance-based prioritization
- **Use Case:** Validation against known pathogenic variants, quality control
- **Advantages:** Real clinical evidence, regulatory compliance, expert curation

### VEP-Annotated Patient Data
- **Data Source:** Patient WGS/WES with VEP computational annotations
- **Results Location:** `results/vep_annotated/`
- **Scoring Strategy:** Multi-tier computational prediction + frequency filtering
- **Use Case:** Novel variant discovery, rare disease gene identification
- **Advantages:** Comprehensive annotation, population context, prediction algorithms

### Population Genetics Data
- **Data Source:** 1000 Genomes, gnomAD, UK Biobank cohorts
- **Results Location:** `results/population_data/`
- **Scoring Strategy:** Allele frequency-based rarity assessment
- **Use Case:** Population-scale analysis, frequency validation
- **Advantages:** Large sample sizes, diverse populations, allele frequency precision

## Pipeline Components

### Variant Prioritization (rare_disease_filter.nf)
**Universal Parser with Auto-Detection:**
- Detects VCF type: ClinVar (CLNSIG), VEP (CSQ), Population (AF), Basic
- Routes to appropriate parsing strategy automatically
- Extracts relevant annotations based on data source
- Applies data-type-specific scoring algorithms

**Scoring Systems:**
- **ClinVar:** Clinical significance (Pathogenic=10pts) + consequence + gene panel
- **VEP:** CADD + REVEL + frequency + consequence + clinical evidence
- **Population:** Frequency rarity + basic consequence prediction
- **Unified output:** Consistent TSV format across all data types

**Gene Panel:** BRCA1, BRCA2, TP53, MLH1, MSH2, APC, PTEN, STK11, CDH1, CHEK2, ATM, PALB2, BRIP1, RAD51C, RAD51D, CDKN2A

### Clinical Reporting (clinical_report.nf)
- **Interactive HTML dashboards** with Plotly.js visualizations
- **Executive summary** with priority distribution and clinical interpretation
- **VCF type identification** displayed prominently in reports
- **Adaptive content** based on data source (clinical vs computational)
- **Professional styling** with Bootstrap framework for clinical review

### Inheritance Analysis (inheritance_pattern.nf)
- **PED file parsing** for family-based rare disease analysis
- **De novo variant identification** in parent-child trios
- **Compound heterozygous detection** within genes
- **Inheritance pattern classification** (AR, AD, X-linked)

### Workflow Orchestration (rare_disease_enhanced.nf)
- **Parallel processing** of prioritization, inheritance, and reporting
- **Error handling** with retry logic and comprehensive logging
- **Resource management** with configurable CPU/memory allocation
- **Multi-platform deployment** (local, cloud, HPC)

## Quick Start

### ClinVar Clinical Data Analysis
```bash
# Download ClinVar VCF
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
gunzip -c clinvar.vcf.gz | head -n 5000 > test_data/clinvar_test.vcf

# Run with clinical data
nextflow run workflows/rare_disease_enhanced.nf \
  --input test_data/clinvar_test.vcf \
  --sample_id "CLINVAR_ANALYSIS"
```

### Patient WGS/WES Analysis (VEP-annotated)
```bash
# Requires VEP-annotated VCF with CSQ field
nextflow run workflows/rare_disease_enhanced.nf \
  --input patient_data/sample_vep_annotated.vcf \
  --sample_id "PATIENT_WGS" \
  --pedigree family.ped
```

### Population Cohort Analysis
```bash
# Works with 1000G, gnomAD, or similar population VCFs
nextflow run workflows/rare_disease_enhanced.nf \
  --input population_data/cohort_chr22.vcf \
  --sample_id "POPULATION_ANALYSIS"
```

## Output Files
The pipeline generates results in work directories, organized by analysis type:

- `{SAMPLE_ID}.prioritized.tsv`: Ranked variants with comprehensive scoring
- `{SAMPLE_ID}.summary.json`: Analysis summary with VCF type and statistics
- `{SAMPLE_ID}.clinical_report.html`: Interactive clinical dashboard
- `{SAMPLE_ID}.inheritance.json`: Family-based inheritance analysis
- `{SAMPLE_ID}.compound_het.tsv`: Compound heterozygous variant pairs
- `{SAMPLE_ID}.de_novo.tsv`: Potential de novo variants in trios

## Scoring Methodology

### ClinVar Clinical Scoring
```
Priority Score = Clinical Significance + Consequence + Gene Panel + Frequency
- Clinical Significance: Pathogenic=10, Likely_pathogenic=8, VUS=2, Benign=0
- Consequence Impact: High=4 (frameshift, stop), Moderate=2 (missense)
- Gene Panel Bonus: Disease gene=2, Other=0
- Frequency Bonus: Unknown=3 (assumes rare), Known=frequency-based
```

### VEP Computational Scoring
```
Priority Score = Pathogenicity + Frequency + Consequence + Gene Panel
- Pathogenicity: CADD(≥30)=4 + REVEL(≥0.75)=3 + ClinVar evidence
- Frequency: gnomAD AF ≤0.00001=3, ≤0.0001=2, ≤0.001=1, >0.001=0
- Consequence: Same as clinical scoring
- Gene Panel: Same as clinical scoring
```

## Data Quality and Validation

### Real vs Simulated Data
- **ClinVar results:** 100% real clinical significance from expert curation
- **VEP results:** Real computational predictions (CADD, REVEL, SpliceAI)
- **Population results:** Real allele frequencies from large cohorts
- **No simulation:** Pipeline rejects datasets requiring simulated pathogenicity

### Quality Metrics
- **Sensitivity:** Correctly identifies known pathogenic ClinVar variants
- **Specificity:** Appropriately deprioritizes benign/likely benign variants
- **Clinical relevance:** Prioritizes variants in established disease genes
- **Reproducibility:** Consistent results across platform deployments

## Technical Implementation
- **Language:** NextFlow DSL2, Python 3.9+
- **Dependencies:** pandas, numpy, PyVCF3, plotly for visualizations
- **Container:** Python 3.11 slim with scientific computing libraries
- **Resource requirements:** 2-4 CPU cores, 4-8GB RAM per process
- **Scalability:** Parallel processing with configurable resource allocation

## Configuration and Deployment
The pipeline includes configurations for:
- **Local execution** with automatic dependency management
- **Docker containerization** for reproducible research environments
- **Google Cloud Platform** with Batch API and Cloud Storage integration
- **Amazon Web Services** with Batch and S3 data management
- **SLURM clusters** for high-performance computing environments

## File Structure
```
├── modules/prioritization/
│   ├── rare_disease_filter.nf      # Universal parser with auto-detection
│   ├── clinical_report.nf          # Interactive dashboard generation
│   └── inheritance_pattern.nf      # Family-based inheritance analysis
├── workflows/
│   └── rare_disease_enhanced.nf    # Main workflow orchestration
├── test_data/
│   ├── clinvar_test.vcf            # ClinVar clinical subset
│   └── test_subset.vcf             # Population genetics demo
├── results/
│   ├── real_clinvar/              # Clinical analysis results
│   ├── vep_annotated/             # Patient genomic analysis
│   └── population_data/           # Cohort analysis results
├── nextflow.config                 # Multi-platform configuration
└── README.md
```

## Use Cases and Applications

### Clinical Genomics Laboratory
- **Patient WGS/WES interpretation** for rare disease diagnosis
- **Variant validation** against ClinVar clinical database
- **Family-based inheritance analysis** for genetic counseling
- **Quality control** through known positive/negative controls

### Research Applications
- **Novel disease gene discovery** in rare disease cohorts
- **Phenotype expansion studies** for known disease genes
- **Population-scale variant effect prediction** and validation
- **Therapeutic target identification** through pathway analysis


## Limitations and Considerations
- **VEP dependency:** Optimal performance requires VEP-annotated input data
- **Clinical interpretation:** Results require expert clinical correlation
- **Population bias:** Frequency data may not represent all ancestry groups
- **Computational predictions:** CADD/REVEL scores supplement but don't replace clinical evidence

## Development and Citation
This pipeline demonstrates production-ready clinical genomics workflows suitable for:
- Academic research laboratories
- Clinical genetics laboratories  
- Pharmaceutical research and development
- Regulatory submission and compliance

The implementation showcases modern bioinformatics best practices including containerization, workflow orchestration, interactive reporting, and multi-platform deployment.

## License
MIT License - see LICENSE file for details.
