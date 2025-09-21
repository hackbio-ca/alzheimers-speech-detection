import os
import json
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import plotly.express as px
from plotly.utils import PlotlyJSONEncoder
from flask import Flask, request, jsonify, render_template_string
from flask_cors import CORS
import requests
from scipy import stats
from scipy.stats import hypergeom
import math

app = Flask(__name__)
CORS(app)

# GO enrichment analysis functions
def get_go_terms(gene_list, species='human'):
    """
    Get GO terms for a list of genes using Ensembl API
    """
    # This is a simplified version - in practice you'd use proper GO databases
    # For demo purposes, we'll create mock GO terms
    go_terms = {
        'GO:0008150': {'name': 'biological_process', 'genes': gene_list[:len(gene_list)//3]},
        'GO:0005575': {'name': 'cellular_component', 'genes': gene_list[len(gene_list)//3:2*len(gene_list)//3]},
        'GO:0003674': {'name': 'molecular_function', 'genes': gene_list[2*len(gene_list)//3:]},
        'GO:0006412': {'name': 'translation', 'genes': gene_list[:len(gene_list)//4]},
        'GO:0006355': {'name': 'regulation_of_transcription', 'genes': gene_list[len(gene_list)//4:len(gene_list)//2]},
        'GO:0007049': {'name': 'cell_cycle', 'genes': gene_list[len(gene_list)//2:3*len(gene_list)//4]},
        'GO:0006915': {'name': 'apoptotic_process', 'genes': gene_list[3*len(gene_list)//4:]},
    }
    return go_terms

def perform_go_enrichment(gene_list, background_genes=None, p_value_threshold=0.05):
    """
    Perform GO enrichment analysis
    """
    if background_genes is None:
        # Use a larger background set for demo
        background_genes = [f"GENE_{i:04d}" for i in range(20000)]
    
    go_terms = get_go_terms(background_genes)
    enrichment_results = []
    
    for go_id, go_data in go_terms.items():
        # Count genes in our list that are in this GO term
        genes_in_term = set(go_data['genes'])
        genes_in_list = set(gene_list)
        overlap = genes_in_term.intersection(genes_in_list)
        
        if len(overlap) > 0:
            # Hypergeometric test
            N = len(background_genes)  # Total genes
            K = len(genes_in_term)     # Genes in GO term
            n = len(gene_list)         # Genes in our list
            x = len(overlap)           # Overlap
            
            # Calculate p-value using hypergeometric distribution
            p_value = hypergeom.sf(x-1, N, K, n)
            
            # Calculate fold enrichment
            expected = (K * n) / N
            fold_enrichment = x / expected if expected > 0 else 0
            
            enrichment_results.append({
                'go_id': go_id,
                'go_name': go_data['name'],
                'genes_in_term': K,
                'genes_in_list': n,
                'overlap': x,
                'expected': expected,
                'fold_enrichment': fold_enrichment,
                'p_value': p_value,
                'genes': list(overlap)
            })
    
    # Sort by p-value
    enrichment_results.sort(key=lambda x: x['p_value'])
    
    return enrichment_results

def create_volcano_plot(data):
    """
    Create a volcano plot from gene expression data
    """
    df = pd.DataFrame(data)
    
    # Calculate -log10(p-value)
    df['neg_log_pval'] = -np.log10(df['p_val'] + 1e-300)  # Add small value to avoid log(0)
    
    # Define significance thresholds
    pval_threshold = 0.05
    fc_threshold = 1.0  # log2 fold change threshold
    
    # Color points based on significance
    colors = []
    for _, row in df.iterrows():
        if row['p_val'] < pval_threshold and abs(row['LogFC']) > fc_threshold:
            if row['LogFC'] > 0:
                colors.append('red')  # Upregulated
            else:
                colors.append('blue')  # Downregulated
        else:
            colors.append('gray')  # Not significant
    
    # Create the plot
    fig = go.Figure()
    
    # Add scatter points
    fig.add_trace(go.Scatter(
        x=df['LogFC'],
        y=df['neg_log_pval'],
        mode='markers',
        marker=dict(
            color=colors,
            size=8,
            opacity=0.7,
            line=dict(width=1, color='black')
        ),
        text=df['gene'],
        hovertemplate='<b>%{text}</b><br>' +
                      'LogFC: %{x:.2f}<br>' +
                      'P-value: %{customdata:.2e}<br>' +
                      '-log10(P): %{y:.2f}<extra></extra>',
        customdata=df['p_val'],
        name='Genes'
    ))
    
    # Add significance lines
    fig.add_hline(y=-np.log10(pval_threshold), line_dash="dash", line_color="red", 
                  annotation_text=f"P = {pval_threshold}")
    fig.add_vline(x=fc_threshold, line_dash="dash", line_color="red")
    fig.add_vline(x=-fc_threshold, line_dash="dash", line_color="red")
    
    # Update layout
    fig.update_layout(
        title='Volcano Plot - Gene Expression Analysis',
        xaxis_title='Log2 Fold Change',
        yaxis_title='-log10(P-value)',
        width=800,
        height=600,
        showlegend=False,
        plot_bgcolor='white',
        xaxis=dict(gridcolor='lightgray'),
        yaxis=dict(gridcolor='lightgray')
    )
    
    return fig

def create_go_enrichment_plot(enrichment_results, top_n=20):
    """
    Create a bar plot for GO enrichment results
    """
    # Take top N results
    top_results = enrichment_results[:top_n]
    
    if not top_results:
        # Create empty plot if no results
        fig = go.Figure()
        fig.update_layout(
            title='GO Enrichment Analysis',
            xaxis_title='-log10(P-value)',
            yaxis_title='GO Terms',
            width=800,
            height=600
        )
        return fig
    
    # Prepare data
    go_names = [result['go_name'] for result in top_results]
    neg_log_pvals = [-np.log10(result['p_value'] + 1e-300) for result in top_results]
    fold_enrichments = [result['fold_enrichment'] for result in top_results]
    
    # Create color scale based on fold enrichment
    colors = px.colors.sequential.Viridis
    
    fig = go.Figure()
    
    fig.add_trace(go.Bar(
        x=neg_log_pvals,
        y=go_names,
        orientation='h',
        marker=dict(
            color=fold_enrichments,
            colorscale='Viridis',
            colorbar=dict(title="Fold Enrichment"),
            line=dict(color='black', width=1)
        ),
        text=[f"FC: {fe:.2f}" for fe in fold_enrichments],
        textposition='inside',
        hovertemplate='<b>%{y}</b><br>' +
                      '-log10(P): %{x:.2f}<br>' +
                      'Fold Enrichment: %{marker.color:.2f}<extra></extra>'
    ))
    
    fig.update_layout(
        title='GO Enrichment Analysis',
        xaxis_title='-log10(P-value)',
        yaxis_title='GO Terms',
        width=800,
        height=600,
        plot_bgcolor='white',
        yaxis=dict(autorange="reversed"),
        xaxis=dict(gridcolor='lightgray'),
        yaxis_gridcolor='lightgray'
    )
    
    return fig

@app.route("/api/analyze_expression", methods=['POST'])
def analyze_expression():
    """
    Main endpoint for gene expression analysis
    """
    try:
        data = request.get_json()
        
        if not data or 'genes' not in data:
            return jsonify({'error': 'No gene data provided'}), 400
        
        gene_data = data['genes']
        
        # Validate data format
        required_fields = ['gene', 'LogFC', 'p_val']
        for item in gene_data:
            for field in required_fields:
                if field not in item:
                    return jsonify({'error': f'Missing field: {field}'}), 400
        
        # Create volcano plot
        volcano_fig = create_volcano_plot(gene_data)
        volcano_json = json.dumps(volcano_fig, cls=PlotlyJSONEncoder)
        
        # Get significant genes for GO enrichment
        significant_genes = [
            item['gene'] for item in gene_data 
            if item['p_val'] < 0.05 and abs(item['LogFC']) > 1.0
        ]
        
        # Perform GO enrichment
        enrichment_results = perform_go_enrichment(significant_genes)
        
        # Create GO enrichment plot
        go_fig = create_go_enrichment_plot(enrichment_results)
        go_json = json.dumps(go_fig, cls=PlotlyJSONEncoder)
        
        # Prepare response
        response = {
            'volcano_plot': json.loads(volcano_json),
            'go_enrichment_plot': json.loads(go_json),
            'enrichment_results': enrichment_results,
            'significant_genes_count': len(significant_genes),
            'total_genes_count': len(gene_data)
        }
        
        return jsonify(response)
        
    except Exception as e:
        return jsonify({'error': str(e)}), 500

@app.route("/api/volcano_plot", methods=['POST'])
def volcano_plot():
    """
    Endpoint specifically for volcano plot generation
    """
    try:
        data = request.get_json()
        
        if not data or 'genes' not in data:
            return jsonify({'error': 'No gene data provided'}), 400
        
        gene_data = data['genes']
        
        # Create volcano plot
        fig = create_volcano_plot(gene_data)
        plot_json = json.dumps(fig, cls=PlotlyJSONEncoder)
        
        return jsonify({'plot': json.loads(plot_json)})
        
    except Exception as e:
        return jsonify({'error': str(e)}), 500

@app.route("/api/go_enrichment", methods=['POST'])
def go_enrichment():
    """
    Endpoint specifically for GO enrichment analysis
    """
    try:
        data = request.get_json()
        
        if not data or 'genes' not in data:
            return jsonify({'error': 'No gene data provided'}), 400
        
        gene_list = data['genes']
        
        # Perform GO enrichment
        enrichment_results = perform_go_enrichment(gene_list)
        
        # Create GO enrichment plot
        fig = create_go_enrichment_plot(enrichment_results)
        plot_json = json.dumps(fig, cls=PlotlyJSONEncoder)
        
        response = {
            'plot': json.loads(plot_json),
            'enrichment_results': enrichment_results
        }
        
        return jsonify(response)
        
    except Exception as e:
        return jsonify({'error': str(e)}), 500

@app.route("/api/example_data")
def example_data():
    """
    Provide example gene expression data for testing
    """
    # Generate example data
    np.random.seed(42)
    n_genes = 1000
    
    genes = [f"GENE_{i:04d}" for i in range(n_genes)]
    log_fc = np.random.normal(0, 2, n_genes)
    p_vals = np.random.uniform(0.001, 0.1, n_genes)
    
    # Make some genes significant
    significant_indices = np.random.choice(n_genes, size=50, replace=False)
    log_fc[significant_indices] = np.random.normal(0, 3, 50)
    p_vals[significant_indices] = np.random.uniform(0.001, 0.01, 50)
    
    example_data = [
        {
            'gene': gene,
            'LogFC': float(log_fc[i]),
            'p_val': float(p_vals[i])
        }
        for i, gene in enumerate(genes)
    ]
    
    return jsonify({'genes': example_data})

@app.route("/")
def index():
    """
    Simple HTML page for testing the API
    """
    html = """
    <!DOCTYPE html>
    <html>
    <head>
        <title>Gene Expression Analysis API</title>
        <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
        <style>
            body { font-family: Arial, sans-serif; margin: 40px; }
            .container { max-width: 1200px; margin: 0 auto; }
            .plot-container { margin: 20px 0; }
            button { padding: 10px 20px; margin: 10px; background: #007bff; color: white; border: none; border-radius: 5px; cursor: pointer; }
            button:hover { background: #0056b3; }
            .info { background: #f8f9fa; padding: 15px; border-radius: 5px; margin: 10px 0; }
        </style>
    </head>
    <body>
        <div class="container">
            <h1>Gene Expression Analysis API</h1>
            <div class="info">
                <h3>Available Endpoints:</h3>
                <ul>
                    <li><strong>POST /api/analyze_expression</strong> - Full analysis (volcano + GO enrichment)</li>
                    <li><strong>POST /api/volcano_plot</strong> - Volcano plot only</li>
                    <li><strong>POST /api/go_enrichment</strong> - GO enrichment only</li>
                    <li><strong>GET /api/example_data</strong> - Get example data</li>
                </ul>
            </div>
            
            <button onclick="loadExampleData()">Load Example Data</button>
            <button onclick="runFullAnalysis()">Run Full Analysis</button>
            
            <div id="volcano-plot" class="plot-container"></div>
            <div id="go-plot" class="plot-container"></div>
            <div id="results"></div>
        </div>

        <script>
            let currentData = null;

            async function loadExampleData() {
                try {
                    const response = await fetch('/api/example_data');
                    const data = await response.json();
                    currentData = data;
                    document.getElementById('results').innerHTML = 
                        `<div class="info">Loaded ${data.genes.length} genes from example data</div>`;
                } catch (error) {
                    console.error('Error loading example data:', error);
                }
            }

            async function runFullAnalysis() {
                if (!currentData) {
                    alert('Please load example data first');
                    return;
                }

                try {
                    const response = await fetch('/api/analyze_expression', {
                        method: 'POST',
                        headers: {
                            'Content-Type': 'application/json',
                        },
                        body: JSON.stringify(currentData)
                    });

                    const result = await response.json();
                    
                    if (result.error) {
                        throw new Error(result.error);
                    }

                    // Display volcano plot
                    Plotly.newPlot('volcano-plot', result.volcano_plot.data, result.volcano_plot.layout);
                    
                    // Display GO enrichment plot
                    Plotly.newPlot('go-plot', result.go_enrichment_plot.data, result.go_enrichment_plot.layout);
                    
                    // Display results summary
                    document.getElementById('results').innerHTML = `
                        <div class="info">
                            <h3>Analysis Results:</h3>
                            <p>Total genes: ${result.total_genes_count}</p>
                            <p>Significant genes: ${result.significant_genes_count}</p>
                            <p>GO terms enriched: ${result.enrichment_results.length}</p>
                        </div>
                    `;

                } catch (error) {
                    console.error('Error running analysis:', error);
                    document.getElementById('results').innerHTML = 
                        `<div class="info" style="color: red;">Error: ${error.message}</div>`;
                }
            }
        </script>
    </body>
    </html>
    """
    return html

if __name__ == "__main__":
    app.run(host="0.0.0.0", port=5001, debug=True)
