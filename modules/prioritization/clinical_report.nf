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
    
    """
    pip install pandas plotly

    cat > report_generator.py << 'EOF'
import pandas as pd
import json
import logging
from datetime import datetime
import plotly.graph_objects as go

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def load_data():
    variants_df = pd.read_csv('${prioritized_variants}', sep='\\t')
    with open('${summary_json}', 'r') as f:
        summary = json.load(f)
    return variants_df, summary

def create_priority_distribution_plot(variants_df):
    priority_bins = ['Low (0-4)', 'Medium (5-7)', 'High (8+)']
    counts = [
        len(variants_df[variants_df['priority_score'] < 5]),
        len(variants_df[(variants_df['priority_score'] >= 5) & (variants_df['priority_score'] < 8)]),
        len(variants_df[variants_df['priority_score'] >= 8])
    ]
    
    colors = ['#28a745', '#fd7e14', '#dc3545']
    
    fig = go.Figure()
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

def format_variant_table(variants_df):
    top_variants = variants_df.head(20).copy()
    
    # Format numeric columns
    for col in ['priority_score']:
        if col in top_variants.columns:
            top_variants[col] = top_variants[col].apply(lambda x: f'{x:.1f}' if pd.notna(x) else 'N/A')
    
    return top_variants.to_html(classes='table table-striped table-hover', table_id='variants-table', escape=False)

def generate_clinical_interpretation(variants_df, summary):
    high_priority = len(variants_df[variants_df['priority_score'] >= 8])
    
    interpretation = []
    
    if high_priority == 0:
        interpretation.append('No high-priority variants identified.')
        interpretation.append('Consider broader analysis or alternative approaches.')
    elif high_priority <= 5:
        interpretation.append(f'{high_priority} high-priority variant(s) identified requiring follow-up.')
        interpretation.append('Recommend further validation and clinical correlation.')
    else:
        interpretation.append(f'{high_priority} high-priority variants identified.')
        interpretation.append('Multiple variants suggest complex inheritance or phenotypic overlap.')
    
    return interpretation

def create_html_report(variants_df, summary):
    priority_plot = create_priority_distribution_plot(variants_df)
    variant_table = format_variant_table(variants_df)
    interpretation = generate_clinical_interpretation(variants_df, summary)
    
    # Handle different summary field names
    total_variants = summary.get('total_variants', summary.get('total_prioritized_variants', len(variants_df)))
    high_priority = summary.get('high_priority_variants', len(variants_df[variants_df['priority_score'] >= 8]))
    medium_priority = summary.get('medium_priority_variants', len(variants_df[(variants_df['priority_score'] >= 5) & (variants_df['priority_score'] < 8)]))
    low_priority = summary.get('low_priority_variants', len(variants_df[variants_df['priority_score'] < 5]))
    
    html_template = f'''
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Clinical Genomics Report - {summary['sample_id']}</title>
    <link href="https://cdn.jsdelivr.net/npm/bootstrap@5.1.3/dist/css/bootstrap.min.css" rel="stylesheet">
    <style>
        .header-section {{ background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); color: white; }}
        .priority-high {{ background-color: #dc3545; color: white; }}
        .priority-medium {{ background-color: #fd7e14; color: white; }}
        .priority-low {{ background-color: #28a745; color: white; }}
        .interpretation-box {{ background-color: #f8f9fa; border-left: 4px solid #007bff; }}
    </style>
</head>
<body>
    <div class="container-fluid">
        <!-- Header -->
        <div class="row header-section p-4 mb-4">
            <div class="col-md-8">
                <h1>Clinical Genomics Report</h1>
                <h2>Sample: {summary['sample_id']}</h2>
                <p>VCF Type: {summary.get('vcf_type', 'Unknown')}</p>
            </div>
            <div class="col-md-4 text-end">
                <p><strong>Report Date:</strong> {datetime.now().strftime('%Y-%m-%d %H:%M')}</p>
                <p><strong>Pipeline:</strong> Rare Disease Variant Prioritization v2.0</p>
            </div>
        </div>

        <!-- Executive Summary -->
        <div class="row mb-4">
            <div class="col-12">
                <h3>Executive Summary</h3>
                <div class="row">
                    <div class="col-md-3">
                        <div class="card priority-high text-center p-3">
                            <h4>{high_priority}</h4>
                            <p>High Priority Variants</p>
                        </div>
                    </div>
                    <div class="col-md-3">
                        <div class="card priority-medium text-center p-3">
                            <h4>{medium_priority}</h4>
                            <p>Medium Priority Variants</p>
                        </div>
                    </div>
                    <div class="col-md-3">
                        <div class="card priority-low text-center p-3">
                            <h4>{low_priority}</h4>
                            <p>Low Priority Variants</p>
                        </div>
                    </div>
                    <div class="col-md-3">
                        <div class="card bg-info text-white text-center p-3">
                            <h4>{total_variants}</h4>
                            <p>Total Variants</p>
                        </div>
                    </div>
                </div>
            </div>
        </div>

        <!-- Clinical Interpretation -->
        <div class="row mb-4">
            <div class="col-12">
                <h3>Clinical Interpretation</h3>
                <div class="interpretation-box p-3">
                    {'<br>'.join(interpretation)}
                </div>
            </div>
        </div>

        <!-- Visualization -->
        <div class="row mb-4">
            <div class="col-12">
                {priority_plot}
            </div>
        </div>

        <!-- Top Variants -->
        <div class="row mb-4">
            <div class="col-12">
                <h3>Top 20 Prioritized Variants</h3>
                <div class="table-responsive">
                    {variant_table}
                </div>
            </div>
        </div>

        <!-- Footer -->
        <div class="row bg-light p-3 mt-4">
            <div class="col-12 text-center">
                <p><small>This report is for research purposes. Clinical decisions should involve genetic counseling.</small></p>
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

# Simple PDF placeholder
with open('${prefix}.clinical_report.pdf', 'w') as f:
    f.write('PDF generation not available in this environment')

logger.info('Clinical report generation complete')
EOF

    python3 report_generator.py

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
