#!/usr/bin/env python3
"""
Example client for the Gene Expression Analysis API
Demonstrates how to use the backend from Python
"""

import requests
import json
import time

# Server configuration
BASE_URL = "http://localhost:5001"

def test_api_connection():
    """Test if the API server is running"""
    try:
        response = requests.get(f"{BASE_URL}/api/example_data", timeout=5)
        if response.status_code == 200:
            print("✅ API server is running")
            return True
        else:
            print(f"❌ API server returned status code: {response.status_code}")
            return False
    except requests.exceptions.RequestException as e:
        print(f"❌ Cannot connect to API server: {e}")
        print("Make sure the server is running with: python gene_expression_server.py")
        return False

def load_example_data():
    """Load example gene expression data from the API"""
    print("\n📊 Loading example data...")
    try:
        response = requests.get(f"{BASE_URL}/api/example_data")
        response.raise_for_status()
        data = response.json()
        print(f"✅ Loaded {len(data['genes'])} genes")
        return data
    except requests.exceptions.RequestException as e:
        print(f"❌ Error loading example data: {e}")
        return None

def run_volcano_plot_analysis(gene_data):
    """Run volcano plot analysis only"""
    print("\n🌋 Running volcano plot analysis...")
    try:
        response = requests.post(
            f"{BASE_URL}/api/volcano_plot",
            json=gene_data,
            headers={'Content-Type': 'application/json'}
        )
        response.raise_for_status()
        result = response.json()
        print("✅ Volcano plot generated successfully")
        return result
    except requests.exceptions.RequestException as e:
        print(f"❌ Error generating volcano plot: {e}")
        return None

def run_go_enrichment_analysis(gene_data):
    """Run GO enrichment analysis only"""
    print("\n🧬 Running GO enrichment analysis...")
    
    # Extract significant genes for GO enrichment
    significant_genes = [
        item['gene'] for item in gene_data['genes'] 
        if item['p_val'] < 0.05 and abs(item['LogFC']) > 1.0
    ]
    
    if not significant_genes:
        print("⚠️  No significant genes found for GO enrichment")
        return None
    
    print(f"Found {len(significant_genes)} significant genes for GO enrichment")
    
    try:
        response = requests.post(
            f"{BASE_URL}/api/go_enrichment",
            json={'genes': significant_genes},
            headers={'Content-Type': 'application/json'}
        )
        response.raise_for_status()
        result = response.json()
        print("✅ GO enrichment analysis completed")
        return result
    except requests.exceptions.RequestException as e:
        print(f"❌ Error running GO enrichment: {e}")
        return None

def run_full_analysis(gene_data):
    """Run complete gene expression analysis"""
    print("\n🔬 Running full gene expression analysis...")
    try:
        response = requests.post(
            f"{BASE_URL}/api/analyze_expression",
            json=gene_data,
            headers={'Content-Type': 'application/json'}
        )
        response.raise_for_status()
        result = response.json()
        
        print("✅ Full analysis completed successfully")
        print(f"   📈 Total genes: {result['total_genes_count']}")
        print(f"   🎯 Significant genes: {result['significant_genes_count']}")
        print(f"   🧬 GO terms enriched: {len(result['enrichment_results'])}")
        
        return result
    except requests.exceptions.RequestException as e:
        print(f"❌ Error running full analysis: {e}")
        return None

def display_enrichment_results(enrichment_results, top_n=5):
    """Display top GO enrichment results"""
    if not enrichment_results:
        print("No enrichment results to display")
        return
    
    print(f"\n📋 Top {min(top_n, len(enrichment_results))} GO Enrichment Results:")
    print("-" * 80)
    print(f"{'GO Term':<30} {'P-value':<12} {'Fold Enrichment':<15} {'Genes':<8}")
    print("-" * 80)
    
    for i, result in enumerate(enrichment_results[:top_n]):
        go_name = result['go_name'][:28] + ".." if len(result['go_name']) > 30 else result['go_name']
        p_val = f"{result['p_value']:.2e}"
        fold_enrichment = f"{result['fold_enrichment']:.2f}"
        overlap = result['overlap']
        
        print(f"{go_name:<30} {p_val:<12} {fold_enrichment:<15} {overlap:<8}")

def save_plots_to_html(result, filename="gene_expression_analysis.html"):
    """Save plots to an HTML file for viewing"""
    if not result:
        print("No results to save")
        return
    
    html_template = """
    <!DOCTYPE html>
    <html>
    <head>
        <title>Gene Expression Analysis Results</title>
        <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
        <style>
            body { font-family: Arial, sans-serif; margin: 40px; }
            .plot-container { margin: 20px 0; }
            .info { background: #f8f9fa; padding: 15px; border-radius: 5px; margin: 10px 0; }
        </style>
    </head>
    <body>
        <h1>Gene Expression Analysis Results</h1>
        <div class="info">
            <h3>Analysis Summary:</h3>
            <p>Total genes: {total_genes}</p>
            <p>Significant genes: {significant_genes}</p>
            <p>GO terms enriched: {go_terms}</p>
        </div>
        
        <div id="volcano-plot" class="plot-container"></div>
        <div id="go-plot" class="plot-container"></div>
        
        <script>
            // Volcano plot
            Plotly.newPlot('volcano-plot', {volcano_data}, {volcano_layout});
            
            // GO enrichment plot
            Plotly.newPlot('go-plot', {go_data}, {go_layout});
        </script>
    </body>
    </html>
    """
    
    try:
        with open(filename, 'w') as f:
            f.write(html_template.format(
                total_genes=result.get('total_genes_count', 'N/A'),
                significant_genes=result.get('significant_genes_count', 'N/A'),
                go_terms=len(result.get('enrichment_results', [])),
                volcano_data=json.dumps(result['volcano_plot']['data']),
                volcano_layout=json.dumps(result['volcano_plot']['layout']),
                go_data=json.dumps(result['go_enrichment_plot']['data']),
                go_layout=json.dumps(result['go_enrichment_plot']['layout'])
            ))
        print(f"✅ Plots saved to {filename}")
    except Exception as e:
        print(f"❌ Error saving plots: {e}")

def main():
    """Main function to run the example client"""
    print("🧬 Gene Expression Analysis API Client")
    print("=" * 50)
    
    # Test API connection
    if not test_api_connection():
        return
    
    # Load example data
    gene_data = load_example_data()
    if not gene_data:
        return
    
    # Run individual analyses
    volcano_result = run_volcano_plot_analysis(gene_data)
    go_result = run_go_enrichment_analysis(gene_data)
    
    # Run full analysis
    full_result = run_full_analysis(gene_data)
    
    # Display results
    if full_result and full_result.get('enrichment_results'):
        display_enrichment_results(full_result['enrichment_results'])
    
    # Save plots to HTML file
    if full_result:
        save_plots_to_html(full_result)
    
    print("\n🎉 Analysis complete!")
    print("You can also view the interactive interface at: http://localhost:5001")

if __name__ == "__main__":
    main()
