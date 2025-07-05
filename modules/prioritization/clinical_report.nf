process CLINICAL_REPORT {
    tag "$meta.id"
    label 'process_low'
    
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11' :
        'python:3.11-slim' }"

    input:
    tuple val(meta), path(prioritized_variants)
    tuple val(meta), path(summary_json)
    path(pedigree)

    output:
    tuple val(meta), path("*.clinical_report.html"), emit: html_report
    tuple val(meta), path("*.clinical_report.pdf"),  emit: pdf_report
    path "versions.yml",                             emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def pedigree_arg = pedigree ? "-p $pedigree" : ""
    
    """
    pip install pandas jinja2 matplotlib seaborn plotly

    python3 -c "
import pandas as pd
import json
import logging
from datetime import datetime
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
import base64
from io import BytesIO

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def load_data():
    variants_df = pd.read_csv('${prioritized_variants}', sep=chr(9))
    with open('${summary_json}', 'r') as f:
        summary = json.load(f)
    return variants_df, summary

def create_priority_distribution_plot(variants_df):
    fig = go.Figure()
    
    priority_bins = ['Low (0-4)', 'Medium (5-7)', 'High (8-10)', 'Critical (>10)']
    counts = [
        len(variants_df[variants_df['priority_score'] < 5]),
        len(variants_df[(variants_df['priority_score'] >= 5) & (variants_df['priority_score'] < 8)]),
        len(variants_df[(variants_df['priority_score'] >= 8) & (variants_df['priority_score'] <= 10)]),
        len(variants_df[variants_df['priority_score'] > 10])
    ]
    
    colors = ['#2E8B57', '#FF8C00', '#DC143C', '#8B0000']
    
    fig.add_trace(go.Bar(
        x=priority_bins,
        y=counts,
        marker_color=colors,
        text=counts,
        textposition='auto'
    ))
    
    fig.update_layout(
        title='Variant Priority Distribution',
        xaxis_title='Priority Category',
        yaxis_title='Number of Variants',
        template='plotly_white',
        height=400
    )
    
    return fig.to_html(include_plotlyjs='inline', div_id='priority_dist')

def create_consequence_plot(summary):
    consequences = summary.get('consequence_distribution', {})
    if not consequences:
        return '<p>No consequence data available</p>'
    
    fig = go.Figure(data=[go.Pie(
        labels=list(consequences.keys()),
        values=list(consequences.values()),
        hole=0.3
    )])
    
    fig.update_layout(
        title='Variant Consequence Distribution',
        template='plotly_white',
        height=400
    )
    
    return fig.to_html(include_plotlyjs='inline', div_id='consequence_dist')

def format_variant_table(variants_df):
    top_variants = variants_df.head(20).copy()
    
    for col in ['gnomAD_AF', 'CADD_PHRED', 'REVEL']:
        if col in top_variants.columns:
            top_variants[col] = pd.to_numeric(top_variants[col], errors='coerce')
            top_variants[col] = top_variants[col].apply(lambda x: f'{x:.4f}' if pd.notna(x) else 'N/A')
    
    top_variants['priority_score'] = top_variants['priority_score'].apply(lambda x: f'{x:.1f}')
    
    return top_variants.to_html(classes='table table-striped table-hover', table_id='variants-table', escape=False)

def generate_clinical_interpretation(variants_df, summary):
    high_priority = len(variants_df[variants_df['priority_score'] >= 8])
    total_variants = len(variants_df)
    
    interpretation = []
    
    if high_priority == 0:
        interpretation.append('No high-priority pathogenic variants identified in known rare disease genes.')
        interpretation.append('Consider broader gene panel analysis or research-based whole genome interpretation.')
    elif high_priority <= 3:
        interpretation.append(f'{high_priority} high-priority variant(s) identified requiring clinical follow-up.')
        interpretation.append('Recommend genetic counseling and family history assessment.')
    else:
        interpretation.append(f'{high_priority} high-priority variants identified.')
        interpretation.append('Multiple potential pathogenic variants suggest possible compound heterozygosity or phenotypic complexity.')
    
    top_genes = summary.get('top_genes', {})
    if any(gene in ['BRCA1', 'BRCA2', 'TP53'] for gene in top_genes.keys()):
        interpretation.append('Cancer predisposition genes identified. Recommend oncogenetics consultation.')
    
    return interpretation

def create_html_report(variants_df, summary):
    priority_plot = create_priority_distribution_plot(variants_df)
    consequence_plot = create_consequence_plot(summary)
    variant_table = format_variant_table(variants_df)
    interpretation = generate_clinical_interpretation(variants_df, summary)
    
    html_template = f'''
<!DOCTYPE html>
<html lang='en'>
<head>
    <meta charset='UTF-8'>
    <meta name='viewport' content='width=device-width, initial-scale=1.0'>
    <title>Clinical Genomics Report - {summary['sample_id']}</title>
    <link href='https://cdn.jsdelivr.net/npm/bootstrap@5.1.3/dist/css/bootstrap.min.css' rel='stylesheet'>
    <script src='https://cdn.plot.ly/plotly-latest.min.js'></script>
    <style>
        .header-section {{ background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); color: white; }}
        .priority-high {{ background-color: #dc3545; color: white; }}
        .priority-medium {{ background-color: #fd7e14; color: white; }}
        .priority-low {{ background-color: #28a745; color: white; }}
        .interpretation-box {{ background-color: #f8f9fa; border-left: 4px solid #007bff; }}
    </style>
</head>
<body>
    <div class='container-fluid'>
        <!-- Header -->
        <div class='row header-section p-4 mb-4'>
            <div class='col-md-8'>
                <h1>Clinical Genomics Report</h1>
                <h2>Sample: {summary['sample_id']}</h2>
            </div>
            <div class='col-md-4 text-end'>
                <p><strong>Report Date:</strong> {datetime.now().strftime('%Y-%m-%d %H:%M')}</p>
                <p><strong>Pipeline:</strong> Rare Disease Variant Prioritization v1.0</p>
            </div>
        </div>

        <!-- Executive Summary -->
        <div class='row mb-4'>
            <div class='col-12'>
                <h3>Executive Summary</h3>
                <div class='row'>
                    <div class='col-md-3'>
                        <div class='card priority-high text-center p-3'>
                            <h4>{summary['high_priority_variants']}</h4>
                            <p>High Priority Variants</p>
                        </div>
                    </div>
                    <div class='col-md-3'>
                        <div class='card priority-medium text-center p-3'>
                            <h4>{summary['medium_priority_variants']}</h4>
                            <p>Medium Priority Variants</p>
                        </div>
                    </div>
                    <div class='col-md-3'>
                        <div class='card priority-low text-center p-3'>
                            <h4>{summary['low_priority_variants']}</h4>
                            <p>Low Priority Variants</p>
                        </div>
                    </div>
                    <div class='col-md-3'>
                        <div class='card bg-info text-white text-center p-3'>
                            <h4>{summary['total_prioritized_variants']}</h4>
                            <p>Total Variants</p>
                        </div>
                    </div>
                </div>
            </div>
        </div>

        <!-- Clinical Interpretation -->
        <div class='row mb-4'>
            <div class='col-12'>
                <h3>Clinical Interpretation</h3>
                <div class='interpretation-box p-3'>
                    {'<br>'.join(interpretation)}
                </div>
            </div>
        </div>

        <!-- Visualizations -->
        <div class='row mb-4'>
            <div class='col-md-6'>
                {priority_plot}
            </div>
            <div class='col-md-6'>
                {consequence_plot}
            </div>
        </div>

        <!-- Top Prioritized Variants -->
        <div class='row mb-4'>
            <div class='col-12'>
                <h3>Top 20 Prioritized Variants</h3>
                <div class='table-responsive'>
                    {variant_table}
                </div>
            </div>
        </div>

        <!-- Methodology -->
        <div class='row mb-4'>
            <div class='col-12'>
                <h3>Methodology</h3>
                <div class='row'>
                    <div class='col-md-6'>
                        <h5>Prioritization Criteria</h5>
                        <ul>
                            <li><strong>Population Frequency:</strong> gnomAD AF < 0.1%</li>
                            <li><strong>Functional Impact:</strong> High/moderate consequence variants</li>
                            <li><strong>Pathogenicity:</strong> CADD, REVEL, ClinVar scores</li>
                            <li><strong>Gene Panels:</strong> Cancer predisposition genes</li>
                        </ul>
                    </div>
                    <div class='col-md-6'>
                        <h5>Scoring System</h5>
                        <ul>
                            <li><strong>Frequency Score:</strong> 0-3 points (rarer = higher)</li>
                            <li><strong>Consequence Score:</strong> 0-4 points</li>
                            <li><strong>Pathogenicity Score:</strong> 0-10 points</li>
                            <li><strong>Gene Panel Bonus:</strong> 0-2 points</li>
                        </ul>
                    </div>
                </div>
            </div>
        </div>

        <!-- Footer -->
        <div class='row bg-light p-3 mt-4'>
            <div class='col-12 text-center'>
                <p><small>This report is for research purposes. Clinical decisions should involve genetic counseling and additional testing.</small></p>
            </div>
        </div>
    </div>
</body>
</html>
    '''
    
    return html_template

# Main execution
variants_df, summary = load_data()
html_content = create_html_report(variants_df, summary)

with open('${prefix}.clinical_report.html', 'w', encoding='utf-8') as f:
    f.write(html_content)

try:
    import weasyprint
    weasyprint.HTML(string=html_content).write_pdf('${prefix}.clinical_report.pdf')
    logger.info('PDF report generated successfully')
except Exception as e:
    logger.warning(f'PDF generation failed: {e}')
    with open('${prefix}.clinical_report.pdf', 'w') as f:
        f.write('PDF generation not available')

logger.info('Clinical report generation complete')
"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
        plotly: \$(python -c "import plotly; print(plotly.__version__)")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.clinical_report.html
    touch ${prefix}.clinical_report.pdf
    touch versions.yml
    """
}
