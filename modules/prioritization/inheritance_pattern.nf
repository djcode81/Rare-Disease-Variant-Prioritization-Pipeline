process INHERITANCE_ANALYSIS {
    tag "$meta.id"
    label 'process_medium'
    
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11' :
        'python:3.11-slim' }"

    input:
    tuple val(meta), path(vcf)
    path(pedigree)
    path(gene_panels)

    output:
    tuple val(meta), path("*.inheritance.tsv"),    emit: inheritance_patterns
    tuple val(meta), path("*.compound_het.tsv"),  emit: compound_heterozygous
    tuple val(meta), path("*.de_novo.tsv"),       emit: de_novo_variants
    tuple val(meta), path("*.inheritance.json"),  emit: summary
    path "versions.yml",                          emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def pedigree_arg = pedigree ? "-p $pedigree" : ""
    def panels_arg = gene_panels ? "-g $gene_panels" : ""
    
    """
    pip install pandas numpy PyVCF3 networkx

    python3 -c "
import pandas as pd
import numpy as np
import json
import logging
import vcf
from collections import defaultdict, Counter
import itertools

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class PedigreeParser:
    def __init__(self, pedigree_file=None):
        self.families = {}
        self.individuals = {}
        if pedigree_file:
            self.parse_pedigree(pedigree_file)
    
    def parse_pedigree(self, pedigree_file):
        try:
            df = pd.read_csv(pedigree_file, sep=chr(9))
            required_cols = ['FamilyID', 'IndividualID', 'FatherID', 'MotherID', 'Sex', 'Affected']
            
            for _, row in df.iterrows():
                family_id = row['FamilyID']
                individual_id = row['IndividualID']
                
                if family_id not in self.families:
                    self.families[family_id] = []
                
                individual = {
                    'id': individual_id,
                    'father': row['FatherID'] if pd.notna(row['FatherID']) and row['FatherID'] != '0' else None,
                    'mother': row['MotherID'] if pd.notna(row['MotherID']) and row['MotherID'] != '0' else None,
                    'sex': row['Sex'],
                    'affected': row['Affected'] == 1 or row['Affected'] == '1'
                }
                
                self.families[family_id].append(individual)
                self.individuals[individual_id] = individual
                
        except Exception as e:
            logger.warning(f'Could not parse pedigree: {e}')
    
    def get_parents(self, individual_id):
        if individual_id in self.individuals:
            ind = self.individuals[individual_id]
            return ind.get('father'), ind.get('mother')
        return None, None
    
    def is_affected(self, individual_id):
        if individual_id in self.individuals:
            return self.individuals[individual_id]['affected']
        return None

class InheritanceAnalyzer:
    def __init__(self, pedigree_parser=None):
        self.pedigree = pedigree_parser
        self.gene_inheritance = {
            'BRCA1': 'AD', 'BRCA2': 'AD', 'TP53': 'AD',
            'MLH1': 'AD', 'MSH2': 'AD', 'APC': 'AD',
            'PTEN': 'AD', 'STK11': 'AD', 'CDH1': 'AD',
            'CHEK2': 'AD', 'ATM': 'AD', 'PALB2': 'AD',
            'NBN': 'AR', 'FANCC': 'AR', 'FANCA': 'AR'
        }
    
    def load_vcf_with_genotypes(self, vcf_path):
        logger.info('Loading VCF with genotype information...')
        
        vcf_reader = vcf.Reader(open(vcf_path, 'r'))
        variants = []
        
        for record in vcf_reader:
            variant = {
                'CHROM': record.CHROM,
                'POS': record.POS,
                'REF': record.REF,
                'ALT': ','.join([str(alt) for alt in record.ALT]) if record.ALT else '.',
                'SYMBOL': None,
                'Consequence': None,
                'samples': {}
            }
            
            # Parse VEP annotations
            if 'CSQ' in record.INFO:
                csq_data = record.INFO['CSQ']
                if isinstance(csq_data, list) and len(csq_data) > 0:
                    first_csq = csq_data[0]
                    fields = first_csq.split('|')
                    if len(fields) >= 4:
                        variant['Consequence'] = fields[1] if len(fields) > 1 else None
                        variant['SYMBOL'] = fields[3] if len(fields) > 3 else None
            
            # Extract genotypes for all samples
            for sample in record.samples:
                gt = sample['GT'] if hasattr(sample, 'GT') and sample['GT'] else None
                variant['samples'][sample.sample] = {
                    'GT': gt,
                    'AD': getattr(sample, 'AD', None),
                    'DP': getattr(sample, 'DP', None),
                    'GQ': getattr(sample, 'GQ', None)
                }
            
            variants.append(variant)
        
        return pd.DataFrame(variants)
    
    def classify_genotype(self, gt_string):
        if not gt_string or gt_string in ['./.' , '.']:
            return 'missing'
        
        if '|' in gt_string:
            alleles = gt_string.split('|')
        elif '/' in gt_string:
            alleles = gt_string.split('/')
        else:
            return 'unknown'
        
        try:
            alleles = [int(a) for a in alleles if a != '.']
            if len(alleles) != 2:
                return 'missing'
            
            if alleles[0] == 0 and alleles[1] == 0:
                return 'homozygous_ref'
            elif alleles[0] != alleles[1] and (0 in alleles):
                return 'heterozygous'
            elif alleles[0] == alleles[1] and alleles[0] != 0:
                return 'homozygous_alt'
            elif alleles[0] != alleles[1] and (0 not in alleles):
                return 'compound_heterozygous'
            else:
                return 'unknown'
        except ValueError:
            return 'unknown'
    
    def find_de_novo_variants(self, variants_df):
        logger.info('Identifying de novo variants...')
        
        if not self.pedigree:
            logger.warning('No pedigree available for de novo analysis')
            return pd.DataFrame()
        
        de_novo_variants = []
        
        for _, variant in variants_df.iterrows():
            for sample_id, genotype_data in variant['samples'].items():
                if not genotype_data['GT']:
                    continue
                
                child_gt = self.classify_genotype(genotype_data['GT'])
                if child_gt not in ['heterozygous', 'homozygous_alt']:
                    continue
                
                father_id, mother_id = self.pedigree.get_parents(sample_id)
                if not father_id or not mother_id:
                    continue
                
                father_gt = None
                mother_gt = None
                
                if father_id in variant['samples']:
                    father_gt = self.classify_genotype(variant['samples'][father_id]['GT'])
                if mother_id in variant['samples']:
                    mother_gt = self.classify_genotype(variant['samples'][mother_id]['GT'])
                
                # De novo if both parents are homozygous reference
                if father_gt == 'homozygous_ref' and mother_gt == 'homozygous_ref':
                    de_novo_variants.append({
                        'CHROM': variant['CHROM'],
                        'POS': variant['POS'],
                        'REF': variant['REF'],
                        'ALT': variant['ALT'],
                        'SYMBOL': variant['SYMBOL'],
                        'Consequence': variant['Consequence'],
                        'child_id': sample_id,
                        'child_gt': genotype_data['GT'],
                        'father_id': father_id,
                        'father_gt': variant['samples'].get(father_id, {}).get('GT', 'missing'),
                        'mother_id': mother_id,
                        'mother_gt': variant['samples'].get(mother_id, {}).get('GT', 'missing'),
                        'quality_score': genotype_data.get('GQ', 'N/A')
                    })
        
        return pd.DataFrame(de_novo_variants)
    
    def find_compound_heterozygous(self, variants_df):
        logger.info('Identifying compound heterozygous variants...')
        
        compound_het_pairs = []
        
        # Group variants by gene and sample
        gene_variants = defaultdict(lambda: defaultdict(list))
        
        for idx, variant in variants_df.iterrows():
            if not variant['SYMBOL']:
                continue
            
            for sample_id, genotype_data in variant['samples'].items():
                if not genotype_data['GT']:
                    continue
                
                gt = self.classify_genotype(genotype_data['GT'])
                if gt == 'heterozygous':
                    gene_variants[variant['SYMBOL']][sample_id].append({
                        'variant_idx': idx,
                        'CHROM': variant['CHROM'],
                        'POS': variant['POS'],
                        'REF': variant['REF'],
                        'ALT': variant['ALT'],
                        'Consequence': variant['Consequence'],
                        'GT': genotype_data['GT']
                    })
        
        # Find pairs of heterozygous variants in the same gene
        for gene, samples in gene_variants.items():
            for sample_id, variants in samples.items():
                if len(variants) >= 2:
                    # Check all pairs
                    for var1, var2 in itertools.combinations(variants, 2):
                        # Verify they're on different chromosomes or distant positions
                        if (var1['CHROM'] != var2['CHROM'] or 
                            abs(var1['POS'] - var2['POS']) > 10000):
                            
                            compound_het_pairs.append({
                                'gene': gene,
                                'sample_id': sample_id,
                                'variant1_pos': f'{var1[\"CHROM\"]}:{var1[\"POS\"]}',
                                'variant1_change': f'{var1[\"REF\"]}>{var1[\"ALT\"]}',
                                'variant1_consequence': var1['Consequence'],
                                'variant2_pos': f'{var2[\"CHROM\"]}:{var2[\"POS\"]}',
                                'variant2_change': f'{var2[\"REF\"]}>{var2[\"ALT\"]}',
                                'variant2_consequence': var2['Consequence'],
                                'inheritance_pattern': 'compound_heterozygous'
                            })
        
        return pd.DataFrame(compound_het_pairs)
    
    def analyze_inheritance_patterns(self, variants_df):
        logger.info('Analyzing inheritance patterns...')
        
        inheritance_results = []
        
        for _, variant in variants_df.iterrows():
            if not variant['SYMBOL']:
                continue
            
            gene_inheritance = self.gene_inheritance.get(variant['SYMBOL'], 'unknown')
            
            for sample_id, genotype_data in variant['samples'].items():
                if not genotype_data['GT']:
                    continue
                
                gt = self.classify_genotype(genotype_data['GT'])
                if gt in ['heterozygous', 'homozygous_alt']:
                    
                    # Determine inheritance pattern compatibility
                    pattern_compatible = []
                    
                    if gt == 'heterozygous':
                        pattern_compatible.extend(['autosomal_dominant', 'X_linked_dominant'])
                        if gene_inheritance == 'AR':
                            pattern_compatible.append('autosomal_recessive_carrier')
                    
                    if gt == 'homozygous_alt':
                        pattern_compatible.extend(['autosomal_recessive', 'X_linked_recessive'])
                    
                    # Check affected status if pedigree available
                    affected_status = None
                    if self.pedigree:
                        affected_status = self.pedigree.is_affected(sample_id)
                    
                    inheritance_results.append({
                        'CHROM': variant['CHROM'],
                        'POS': variant['POS'],
                        'REF': variant['REF'],
                        'ALT': variant['ALT'],
                        'SYMBOL': variant['SYMBOL'],
                        'Consequence': variant['Consequence'],
                        'sample_id': sample_id,
                        'genotype': gt,
                        'genotype_string': genotype_data['GT'],
                        'known_inheritance': gene_inheritance,
                        'compatible_patterns': ','.join(pattern_compatible),
                        'affected_status': affected_status,
                        'quality_score': genotype_data.get('GQ', 'N/A')
                    })
        
        return pd.DataFrame(inheritance_results)

def generate_summary_stats(inheritance_df, de_novo_df, compound_het_df):
    summary = {
        'total_variants_analyzed': len(inheritance_df),
        'inheritance_patterns': {
            'autosomal_dominant_compatible': 0,
            'autosomal_recessive_compatible': 0,
            'x_linked_compatible': 0
        },
        'de_novo_variants': len(de_novo_df),
        'compound_heterozygous_pairs': len(compound_het_df),
        'affected_individuals': 0,
        'genes_with_variants': 0
    }
    
    if len(inheritance_df) > 0 and 'compatible_patterns' in inheritance_df.columns:
        summary['inheritance_patterns']['autosomal_dominant_compatible'] = len(inheritance_df[inheritance_df['compatible_patterns'].str.contains('autosomal_dominant', na=False)])
        summary['inheritance_patterns']['autosomal_recessive_compatible'] = len(inheritance_df[inheritance_df['compatible_patterns'].str.contains('autosomal_recessive', na=False)])
        summary['inheritance_patterns']['x_linked_compatible'] = len(inheritance_df[inheritance_df['compatible_patterns'].str.contains('X_linked', na=False)])
        
        if 'affected_status' in inheritance_df.columns:
            summary['affected_individuals'] = len(inheritance_df[inheritance_df['affected_status'] == True])
        
        if 'SYMBOL' in inheritance_df.columns:
            summary['genes_with_variants'] = inheritance_df['SYMBOL'].nunique()
    
    return summary

# Main execution
pedigree_parser = None
if '${pedigree}' and '${pedigree}' != '':
    pedigree_parser = PedigreeParser('${pedigree}')

analyzer = InheritanceAnalyzer(pedigree_parser)
variants_df = analyzer.load_vcf_with_genotypes('${vcf}')

# Perform analyses
inheritance_df = analyzer.analyze_inheritance_patterns(variants_df)
de_novo_df = analyzer.find_de_novo_variants(variants_df)
compound_het_df = analyzer.find_compound_heterozygous(variants_df)

# Generate summary
summary = generate_summary_stats(inheritance_df, de_novo_df, compound_het_df)

# Save results
inheritance_df.to_csv('${prefix}.inheritance.tsv', sep=chr(9), index=False)
de_novo_df.to_csv('${prefix}.de_novo.tsv', sep=chr(9), index=False)
compound_het_df.to_csv('${prefix}.compound_het.tsv', sep=chr(9), index=False)

with open('${prefix}.inheritance.json', 'w') as f:
    json.dump(summary, f, indent=2)

logger.info(f'Inheritance analysis complete: {len(inheritance_df)} variants analyzed')
"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
        PyVCF3: \$(python -c "import vcf; print('1.0.0')")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.inheritance.tsv
    touch ${prefix}.compound_het.tsv
    touch ${prefix}.de_novo.tsv
    touch ${prefix}.inheritance.json
    touch versions.yml
    """
}
