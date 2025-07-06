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
    
    cat > prioritize_script.py << 'PYTHON_EOF'
import pandas as pd
import numpy as np
import json
import logging
import vcf
import re

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class AdaptiveVEPParser:
    def __init__(self):
        self.csq_format = None
        self.field_map = {}
        self.vcf_type = 'unknown'
        
    def detect_vcf_type(self, vcf_reader):
        \"\"\"Detect VCF annotation type and parse CSQ format if present\"\"\"
        
        # Check for ClinVar fields
        if any('CLNSIG' in str(info) for info in vcf_reader.infos.values()):
            self.vcf_type = 'clinvar'
            logger.info('Detected ClinVar VCF format')
            return
            
        # Check for VEP CSQ field
        if 'CSQ' in vcf_reader.infos:
            self.vcf_type = 'vep_annotated'
            self._parse_csq_format(vcf_reader.infos['CSQ'])
            logger.info(f'Detected VEP-annotated VCF format with {len(self.field_map)} fields')
            return
            
        # Check for population frequency fields
        if any(field in vcf_reader.infos for field in ['AF', 'AN', 'AC']):
            self.vcf_type = 'population'
            logger.info('Detected population genetics VCF format')
            return
            
        # Default to basic
        self.vcf_type = 'basic'
        logger.info('Detected basic VCF format (minimal annotations)')
        
    def _parse_csq_format(self, csq_info):
        \"\"\"Parse CSQ header to create field mapping\"\"\"
        description = csq_info.desc
        
        # Extract format from description
        format_match = re.search(r'Format:\s*([^\"]+)', description)
        if not format_match:
            logger.warning('Could not parse CSQ format from header')
            return
            
        fields = [field.strip() for field in format_match.group(1).split('|')]
        
        # Create mapping of field names to indices
        for i, field in enumerate(fields):
            self.field_map[field] = i
            
        logger.info(f'Parsed CSQ format: {len(fields)} fields')
        
        # Log important field positions
        important_fields = ['SYMBOL', 'Consequence', 'gnomAD_AF', 'CADD_PHRED', 'REVEL', 'CLIN_SIG']
        for field in important_fields:
            if field in self.field_map:
                logger.info(f'  {field}: position {self.field_map[field]}')
                
    def get_field_value(self, csq_fields, field_name, default=None):
        \"\"\"Safely extract field value from CSQ annotation\"\"\"
        if field_name not in self.field_map:
            return default
            
        index = self.field_map[field_name]
        if index >= len(csq_fields):
            return default
            
        value = csq_fields[index]
        return value if value and value != '' else default

class AdaptiveVariantAnalyzer:
    def __init__(self):
        self.config = {
            'max_gnomad_af': 0.001,
            'rare_disease_genes': ['BRCA1', 'BRCA2', 'TP53', 'MTOR', 'CDKN2A', 'CHEK2', 'EGFR', 'AGT'],
            'high_impact_consequences': ['stop_gained', 'frameshift_variant', 'splice_donor_variant'],
            'moderate_impact_consequences': ['missense_variant', 'splice_region_variant'],
            'clinvar_pathogenic': ['Pathogenic', 'Likely_pathogenic']
        }
        self.parser = AdaptiveVEPParser()
        
    def load_and_parse_vcf(self, vcf_path):
        logger.info('Loading VCF file with adaptive parsing...')
        
        vcf_reader = vcf.Reader(open(vcf_path, 'r'))
        self.parser.detect_vcf_type(vcf_reader)
        
        records = []
        for record in vcf_reader:
            row = {
                'CHROM': record.CHROM,
                'POS': record.POS,
                'REF': record.REF,
                'ALT': ','.join([str(alt) for alt in record.ALT]) if record.ALT else '.',
                'SYMBOL': None,
                'Consequence': None,
                'gnomAD_AF': None,
                'CADD_PHRED': None,
                'REVEL': None,
                'CLIN_SIG': None,
                'vcf_type': self.parser.vcf_type
            }
            
            # Parse based on VCF type
            if self.parser.vcf_type == 'clinvar':
                self._parse_clinvar_record(record, row)
            elif self.parser.vcf_type == 'vep_annotated':
                self._parse_vep_record(record, row)
            elif self.parser.vcf_type == 'population':
                self._parse_population_record(record, row)
            # Basic VCF keeps default None values
            
            records.append(row)
        
        df = pd.DataFrame(records)
        logger.info(f'Loaded {len(df)} variants as {self.parser.vcf_type} format')
        return df
        
    def _parse_clinvar_record(self, record, row):
        \"\"\"Parse ClinVar-specific fields\"\"\"
        if hasattr(record, 'INFO'):
            if 'CLNSIG' in record.INFO:
                clnsig = record.INFO['CLNSIG']
                row['CLIN_SIG'] = clnsig[0] if isinstance(clnsig, list) else clnsig
            if 'GENEINFO' in record.INFO:
                gene_info = record.INFO['GENEINFO']
                gene_str = gene_info[0] if isinstance(gene_info, list) else gene_info
                row['SYMBOL'] = gene_str.split(':')[0] if ':' in str(gene_str) else str(gene_str)
                
    def _parse_vep_record(self, record, row):
        \"\"\"Parse VEP CSQ annotations adaptively\"\"\"
        if 'CSQ' in record.INFO:
            csq_data = record.INFO['CSQ']
            if isinstance(csq_data, list) and len(csq_data) > 0:
                # Take first annotation (most severe)
                first_csq = csq_data[0]
                fields = first_csq.split('|')
                
                # Extract fields using adaptive mapping
                row['SYMBOL'] = self.parser.get_field_value(fields, 'SYMBOL')
                row['Consequence'] = self.parser.get_field_value(fields, 'Consequence')
                
                # Try multiple field names for frequency
                af_value = (self.parser.get_field_value(fields, 'gnomAD_AF') or 
                           self.parser.get_field_value(fields, 'AF') or
                           self.parser.get_field_value(fields, 'MAX_AF'))
                if af_value:
                    try:
                        row['gnomAD_AF'] = float(af_value)
                    except ValueError:
                        pass
                        
                # Extract prediction scores
                cadd_value = self.parser.get_field_value(fields, 'CADD_PHRED')
                if cadd_value:
                    try:
                        row['CADD_PHRED'] = float(cadd_value)
                    except ValueError:
                        pass
                        
                revel_value = self.parser.get_field_value(fields, 'REVEL')
                if revel_value:
                    try:
                        row['REVEL'] = float(revel_value)
                    except ValueError:
                        pass
                        
                # Clinical significance
                row['CLIN_SIG'] = self.parser.get_field_value(fields, 'CLIN_SIG')
                
    def _parse_population_record(self, record, row):
        \"\"\"Parse population genetics fields\"\"\"
        if hasattr(record, 'INFO'):
            if 'AF' in record.INFO:
                af = record.INFO['AF']
                try:
                    row['gnomAD_AF'] = float(af[0] if isinstance(af, list) else af)
                except (ValueError, TypeError):
                    pass
                    
    def prioritize_variants(self, df):
        logger.info('Calculating adaptive priority scores...')
        
        df['pathogenicity_score'] = 0.0
        df['priority_score'] = 0.0
        
        # Frequency score (0-3 points) - works for all VCF types
        freq = pd.to_numeric(df['gnomAD_AF'], errors='coerce')
        df['priority_score'] += np.where(
            freq.isna(), 3,  # Unknown frequency = assume rare
            np.where(freq <= 0.00001, 3,
            np.where(freq <= 0.0001, 2,
            np.where(freq <= 0.001, 1, 0)))
        )
        
        # Consequence score (0-4 points) - VEP and ClinVar
        df['priority_score'] += df['Consequence'].apply(
            lambda x: 4 if pd.notna(x) and any(c in str(x) for c in self.config['high_impact_consequences'])
            else 2 if pd.notna(x) and any(c in str(x) for c in self.config['moderate_impact_consequences'])
            else 0
        )
        
        # CADD score (0-3 points) - VEP only
        cadd = pd.to_numeric(df['CADD_PHRED'], errors='coerce')
        df['pathogenicity_score'] += np.where(cadd >= 30, 3, np.where(cadd >= 20, 2, np.where(cadd >= 15, 1, 0)))
        
        # REVEL score (0-2 points) - VEP only
        revel = pd.to_numeric(df['REVEL'], errors='coerce')
        df['pathogenicity_score'] += np.where(revel >= 0.75, 2, np.where(revel >= 0.5, 1, 0))
        
        # ClinVar clinical significance (0-10 points) - ClinVar and some VEP
        df['pathogenicity_score'] += df['CLIN_SIG'].apply(
            lambda x: 10 if pd.notna(x) and 'Pathogenic' in str(x) and 'Likely' not in str(x)
            else 8 if pd.notna(x) and 'Likely_pathogenic' in str(x)
            else 2 if pd.notna(x) and 'Uncertain' in str(x)
            else 0
        )
        
        # Gene panel bonus (0-2 points) - All types
        df['priority_score'] += df['SYMBOL'].apply(
            lambda x: 2 if pd.notna(x) and str(x) in self.config['rare_disease_genes'] else 0
        )
        
        # Add pathogenicity to priority
        df['priority_score'] += df['pathogenicity_score']
        
        # Sort by priority
        df = df.sort_values('priority_score', ascending=False)
        
        logger.info(f'Prioritization complete: {len(df)} variants')
        return df
        
    def generate_summary(self, df, sample_id):
        summary = {
            'sample_id': sample_id,
            'vcf_type': df['vcf_type'].iloc[0] if len(df) > 0 else 'unknown',
            'total_prioritized_variants': len(df),
            'high_priority_variants': len(df[df['priority_score'] >= 8]),
            'medium_priority_variants': len(df[(df['priority_score'] >= 5) & (df['priority_score'] < 8)]),
            'low_priority_variants': len(df[df['priority_score'] < 5]),
            'top_genes': df['SYMBOL'].dropna().value_counts().head(10).to_dict(),
            'consequence_distribution': df['Consequence'].value_counts().to_dict()
        }
        return summary

# Main execution
analyzer = AdaptiveVariantAnalyzer()
df = analyzer.load_and_parse_vcf('${vcf}')
df = analyzer.prioritize_variants(df)

# Save results
output_cols = ['CHROM', 'POS', 'REF', 'ALT', 'SYMBOL', 'Consequence', 
               'gnomAD_AF', 'CADD_PHRED', 'REVEL', 'CLIN_SIG', 
               'pathogenicity_score', 'priority_score']
available_cols = [col for col in output_cols if col in df.columns]
result_df = df[available_cols]

result_df.to_csv('${prefix}.prioritized.tsv', sep=chr(9), index=False)

# Generate summary
summary = analyzer.generate_summary(df, '${prefix}')
with open('${prefix}.summary.json', 'w') as f:
    json.dump(summary, f, indent=2)

logger.info(f'Analysis complete: {summary["vcf_type"]} format, {len(df)} variants processed')
PYTHON_EOF

python3 prioritize_script.py

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
