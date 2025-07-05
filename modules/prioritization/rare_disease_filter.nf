process RARE_DISEASE_PRIORITIZE {
    tag "$meta.id"
    label 'process_low'
    
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11' :
        'python:3.11-slim' }"

    input:
    tuple val(meta), path(vcf)
    path(config)

    output:
    tuple val(meta), path("*.prioritized.tsv"), emit: prioritized_variants
    tuple val(meta), path("*.summary.json"),    emit: summary
    path "versions.yml",                        emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def config_arg = config ? "-c $config" : ""
    
    """
    pip install pandas numpy PyVCF3
    
    python3 -c "
import pandas as pd
import numpy as np
import json
import logging
import vcf

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Configuration
config = {
    'max_gnomad_af': 0.001,
    'rare_disease_genes': ['BRCA1', 'BRCA2', 'TP53', 'MTOR', 'CDKN2A', 'CHEK2', 'EGFR'],
    'high_impact_consequences': ['stop_gained', 'frameshift_variant', 'splice_donor_variant'],
    'moderate_impact_consequences': ['missense_variant', 'splice_region_variant'],
    'clinvar_pathogenic': ['Pathogenic', 'Likely_pathogenic']
}

def load_and_parse_vcf(vcf_path):
    logger.info('Loading VCF file using PyVCF...')
    
    # Use PyVCF to parse VCF properly
    vcf_reader = vcf.Reader(open(vcf_path, 'r'))
    
    records = []
    for record in vcf_reader:
        # Extract basic info
        row = {
            'CHROM': record.CHROM,
            'POS': record.POS,
            'REF': record.REF,
            'ALT': ','.join([str(alt) for alt in record.ALT]) if record.ALT else '.',
            'INFO': record.INFO if hasattr(record, 'INFO') else {}
        }
        
        # Parse VEP CSQ if present
        row['SYMBOL'] = None
        row['Consequence'] = None
        row['gnomAD_AF'] = None
        row['CADD_PHRED'] = None
        row['REVEL'] = None
        row['CLIN_SIG'] = None
        
        if 'CSQ' in record.INFO:
            csq_data = record.INFO['CSQ']
            if isinstance(csq_data, list) and len(csq_data) > 0:
                # Take first annotation
                first_csq = csq_data[0]
                fields = first_csq.split('|')
                
                if len(fields) >= 4:
                    row['Consequence'] = fields[1] if len(fields) > 1 else None
                    row['SYMBOL'] = fields[3] if len(fields) > 3 else None
                
                # Extract scores (adjust indices as needed for your VEP format)
                if len(fields) >= 70:
                    try:
                        if len(fields) > 55 and fields[55]:
                            row['gnomAD_AF'] = float(fields[55])
                    except (ValueError, IndexError):
                        pass
                    try:
                        if len(fields) > 60 and fields[60]:
                            row['CLIN_SIG'] = fields[60]
                    except IndexError:
                        pass
                    try:
                        if len(fields) > 68 and fields[68]:
                            row['CADD_PHRED'] = float(fields[68])
                    except (ValueError, IndexError):
                        pass
                    try:
                        if len(fields) > 69 and fields[69]:
                            row['REVEL'] = float(fields[69])
                    except (ValueError, IndexError):
                        pass
        
        records.append(row)
    
    df = pd.DataFrame(records)
    logger.info(f'Loaded {len(df)} variants')
    return df

def prioritize_variants(df):
    logger.info('Calculating priority scores...')
    
    # Initialize scores
    df['pathogenicity_score'] = 0.0
    df['priority_score'] = 0.0
    
    # Frequency score (0-3 points)
    freq = pd.to_numeric(df['gnomAD_AF'], errors='coerce')
    df['priority_score'] += np.where(
        freq.isna(), 3,
        np.where(freq <= 0.00001, 3,
        np.where(freq <= 0.0001, 2,
        np.where(freq <= 0.001, 1, 0)))
    )
    
    # Consequence score (0-4 points)
    df['priority_score'] += df['Consequence'].apply(
        lambda x: 4 if pd.notna(x) and any(c in str(x) for c in config['high_impact_consequences'])
        else 2 if pd.notna(x) and any(c in str(x) for c in config['moderate_impact_consequences'])
        else 0
    )
    
    # CADD score (0-3 points)
    cadd = pd.to_numeric(df['CADD_PHRED'], errors='coerce')
    df['pathogenicity_score'] += np.where(cadd >= 30, 3, np.where(cadd >= 20, 2, np.where(cadd >= 15, 1, 0)))
    
    # REVEL score (0-2 points)
    revel = pd.to_numeric(df['REVEL'], errors='coerce')
    df['pathogenicity_score'] += np.where(revel >= 0.75, 2, np.where(revel >= 0.5, 1, 0))
    
    # ClinVar score (0-5 points)
    df['pathogenicity_score'] += df['CLIN_SIG'].apply(
        lambda x: 5 if pd.notna(x) and any(sig in str(x) for sig in config['clinvar_pathogenic']) else 0
    )
    
    # Gene bonus (0-2 points)
    df['priority_score'] += df['SYMBOL'].apply(
        lambda x: 2 if pd.notna(x) and str(x) in config['rare_disease_genes'] else 0
    )
    
    # Add pathogenicity to priority
    df['priority_score'] += df['pathogenicity_score']
    
    # Sort by priority
    df = df.sort_values('priority_score', ascending=False)
    
    logger.info(f'Prioritization complete: {len(df)} variants')
    return df

def generate_summary(df, sample_id):
    summary = {
        'sample_id': sample_id,
        'total_prioritized_variants': len(df),
        'high_priority_variants': len(df[df['priority_score'] >= 8]),
        'medium_priority_variants': len(df[(df['priority_score'] >= 5) & (df['priority_score'] < 8)]),
        'low_priority_variants': len(df[df['priority_score'] < 5]),
        'top_genes': df['SYMBOL'].dropna().value_counts().head(10).to_dict(),
        'consequence_distribution': df['Consequence'].value_counts().to_dict()
    }
    return summary

# Main execution
df = load_and_parse_vcf('${vcf}')
df = prioritize_variants(df)

# Save results
output_cols = ['CHROM', 'POS', 'REF', 'ALT', 'SYMBOL', 'Consequence', 
               'gnomAD_AF', 'CADD_PHRED', 'REVEL', 'CLIN_SIG', 
               'pathogenicity_score', 'priority_score']
available_cols = [col for col in output_cols if col in df.columns]
result_df = df[available_cols]

result_df.to_csv('${prefix}.prioritized.tsv', sep=chr(9), index=False)

# Generate summary
summary = generate_summary(df, '${prefix}')
with open('${prefix}.summary.json', 'w') as f:
    json.dump(summary, f, indent=2)

print(f'Processed {len(df)} variants')
"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.prioritized.tsv
    touch ${prefix}.summary.json
    touch versions.yml
    """
}
