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
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    pip install pandas numpy PyVCF3

    cat > parser.py << 'EOF'
import pandas as pd
import numpy as np
import json
import logging
import vcf

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def detect_vcf_type(vcf_path):
    vcf_reader = vcf.Reader(open(vcf_path, 'r'))
    
    if 'CLNSIG' in vcf_reader.infos:
        return 'clinvar'
    elif 'CSQ' in vcf_reader.infos:
        return 'vep_annotated'
    else:
        return 'basic'

def parse_clinvar(record):
    info = {}
    if 'CLNSIG' in record.INFO:
        info['CLIN_SIG'] = str(record.INFO['CLNSIG'])
    if 'GENEINFO' in record.INFO:
        geneinfo = record.INFO['GENEINFO']
        if isinstance(geneinfo, list):
            geneinfo = geneinfo[0]
        info['SYMBOL'] = str(geneinfo).split(':')[0]
    if 'MC' in record.INFO:
        mc = record.INFO['MC']
        if isinstance(mc, list):
            mc = mc[0]
        if '|' in str(mc):
            info['Consequence'] = str(mc).split('|')[1]
    return info

def load_vcf(vcf_path):
    vcf_type = detect_vcf_type(vcf_path)
    vcf_reader = vcf.Reader(open(vcf_path, 'r'))
    
    records = []
    for i, record in enumerate(vcf_reader):
        if i >= 10000:  # Reasonable limit
            break
            
        variant = {
            'CHROM': record.CHROM,
            'POS': record.POS,
            'REF': record.REF,
            'ALT': ','.join([str(alt) for alt in record.ALT]) if record.ALT else '.'
        }
        
        if vcf_type == 'clinvar':
            variant.update(parse_clinvar(record))
        
        # Set defaults
        variant.setdefault('SYMBOL', None)
        variant.setdefault('Consequence', None)
        variant.setdefault('CLIN_SIG', None)
        
        records.append(variant)
    
    return pd.DataFrame(records), vcf_type

def score_variants(df, vcf_type):
    df['priority_score'] = 0.0
    
    # Clinical significance scoring
    if vcf_type == 'clinvar':
        df['priority_score'] += df['CLIN_SIG'].apply(
            lambda x: 10 if pd.notna(x) and 'Pathogenic' in str(x) and 'Benign' not in str(x)
            else 8 if pd.notna(x) and 'Likely_pathogenic' in str(x)
            else 2 if pd.notna(x) and 'Uncertain_significance' in str(x)
            else 0
        )
    
    # Consequence scoring
    high_impact = ['stop_gained', 'frameshift_variant', 'splice_donor_variant']
    moderate_impact = ['missense_variant', 'splice_region_variant']
    
    df['priority_score'] += df['Consequence'].apply(
        lambda x: 4 if pd.notna(x) and any(c in str(x) for c in high_impact)
        else 2 if pd.notna(x) and any(c in str(x) for c in moderate_impact)
        else 1 if pd.notna(x)
        else 0
    )
    
    # Gene panel bonus
    disease_genes = ['BRCA1', 'BRCA2', 'TP53', 'MLH1', 'MSH2', 'APC']
    df['priority_score'] += df['SYMBOL'].apply(
        lambda x: 2 if pd.notna(x) and str(x) in disease_genes else 0
    )
    
    # Frequency bonus (unknown = rare)
    df['priority_score'] += 3  # Default high priority for rare variants
    
    return df.sort_values('priority_score', ascending=False)

# Main execution
df, vcf_type = load_vcf('${vcf}')
df = score_variants(df, vcf_type)

# Save results
output_cols = ['CHROM', 'POS', 'REF', 'ALT', 'SYMBOL', 'Consequence', 'CLIN_SIG', 'priority_score']
available_cols = [col for col in output_cols if col in df.columns]
df[available_cols].to_csv('${prefix}.prioritized.tsv', sep='\\t', index=False)

# Summary
summary = {
    'sample_id': '${prefix}',
    'vcf_type': vcf_type,
    'total_variants': len(df),
    'high_priority_variants': len(df[df['priority_score'] >= 8]),
    'medium_priority_variants': len(df[(df['priority_score'] >= 5) & (df['priority_score'] < 8)]),
    'low_priority_variants': len(df[df['priority_score'] < 5])
}

with open('${prefix}.summary.json', 'w') as f:
    json.dump(summary, f, indent=2)

logger.info(f'Processed {len(df)} variants ({vcf_type})')
EOF

    python3 parser.py

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
