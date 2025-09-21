#!/usr/bin/env python3
"""
Example client for the Drug Gene Expression Analysis API
Demonstrates how to analyze drugs and their effects on gene expression
"""

import requests
import json
import time

# Server configuration
BASE_URL = "http://localhost:5001"

def test_api_connection():
    """Test if the API server is running"""
    try:
        response = requests.get(f"{BASE_URL}/api/available_drugs", timeout=5)
        if response.status_code == 200:
            print("✅ API server is running")
            return True
        else:
            print(f"❌ API server returned status code: {response.status_code}")
            return False
    except requests.exceptions.RequestException as e:
        print(f"❌ Cannot connect to API server: {e}")
        print("Make sure the server is running with: python drug_gene_expression_server.py")
        return False

def get_available_drugs():
    """Get list of example drugs that can be analyzed"""
    print("\n💊 Getting available drugs...")
    try:
        response = requests.get(f"{BASE_URL}/api/available_drugs")
        response.raise_for_status()
        data = response.json()
        print(f"✅ Found {len(data['drugs'])} example drugs")
        return data['drugs']
    except requests.exceptions.RequestException as e:
        print(f"❌ Error getting available drugs: {e}")
        return []

def analyze_drug(drug_name):
    """Analyze a drug's effects on gene expression"""
    print(f"\n🔬 Analyzing drug: {drug_name}")
    try:
        response = requests.post(
            f"{BASE_URL}/api/analyze_drug",
            json={'drug_name': drug_name},
            headers={'Content-Type': 'application/json'}
        )
        response.raise_for_status()
        result = response.json()
        
        if result.get('error'):
            raise Exception(result['error'])
        
        print("✅ Drug analysis completed successfully")
        print(f"   🧬 Drug: {result['drug_name']}")
        print(f"   📈 Total genes: {result['total_genes_count']}")
        print(f"   🎯 Significant genes: {result['significant_genes_count']}")
        print(f"   🧬 GO terms enriched: {len(result['enrichment_results'])}")
        
        return result
    except requests.exceptions.RequestException as e:
        print(f"❌ Error analyzing drug: {e}")
        return None

def display_gene_expression_summary(gene_data, top_n=10):
    """Display top gene expression changes"""
    if not gene_data:
        print("No gene expression data to display")
        return
    
    # Sort by absolute log2fc
    sorted_genes = sorted(gene_data, key=lambda x: abs(x['LogFC']), reverse=True)
    
    print(f"\n📊 Top {min(top_n, len(sorted_genes))} Gene Expression Changes:")
    print("-" * 80)
    print(f"{'Gene':<15} {'Log2FC':<10} {'P-value':<12} {'Effect':<20}")
    print("-" * 80)
    
    for gene_data in sorted_genes[:top_n]:
        gene = gene_data['gene']
        log2fc = gene_data['LogFC']
        p_val = gene_data['p_val']
        
        # Categorize effect
        if log2fc > 2.0:
            effect = "Strong Upregulation"
        elif log2fc > 1.0:
            effect = "Moderate Upregulation"
        elif log2fc > 0.5:
            effect = "Weak Upregulation"
        elif log2fc > -0.5:
            effect = "No Change"
        elif log2fc > -1.0:
            effect = "Weak Downregulation"
        elif log2fc > -2.0:
            effect = "Moderate Downregulation"
        else:
            effect = "Strong Downregulation"
        
        print(f"{gene:<15} {log2fc:<10.4f} {p_val:<12.2e} {effect:<20}")

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

def save_drug_analysis_to_html(result, filename="drug_analysis_results.html"):
    """Save drug analysis results to an HTML file for viewing"""
    if not result:
        print("No results to save")
        return
    
    html_template = """
    <!DOCTYPE html>
    <html>
    <head>
        <title>Drug Gene Expression Analysis - {drug_name}</title>
        <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
        <style>
            body {{ font-family: Arial, sans-serif; margin: 40px; }}
            .plot-container {{ margin: 20px 0; }}
            .info {{ background: #f8f9fa; padding: 15px; border-radius: 5px; margin: 10px 0; }}
            .drug-header {{ background: #e3f2fd; padding: 20px; border-radius: 5px; margin: 20px 0; }}
        </style>
    </head>
    <body>
        <div class="drug-header">
            <h1>🧬 Drug Gene Expression Analysis</h1>
            <h2>Drug: {drug_name}</h2>
        </div>
        
        <div class="info">
            <h3>Analysis Summary:</h3>
            <p>Total genes analyzed: {total_genes}</p>
            <p>Significantly affected genes: {significant_genes}</p>
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
                drug_name=result.get('drug_name', 'Unknown'),
                total_genes=result.get('total_genes_count', 'N/A'),
                significant_genes=result.get('significant_genes_count', 'N/A'),
                go_terms=len(result.get('enrichment_results', [])),
                volcano_data=json.dumps(result['volcano_plot']['data']),
                volcano_layout=json.dumps(result['volcano_plot']['layout']),
                go_data=json.dumps(result['go_enrichment_plot']['data']),
                go_layout=json.dumps(result['go_enrichment_plot']['layout'])
            ))
        print(f"✅ Analysis results saved to {filename}")
    except Exception as e:
        print(f"❌ Error saving results: {e}")

def run_drug_comparison(drugs):
    """Compare multiple drugs"""
    print(f"\n🔬 Comparing {len(drugs)} drugs...")
    
    results = {}
    for drug in drugs:
        result = analyze_drug(drug)
        if result:
            results[drug] = result
            time.sleep(1)  # Be nice to the server
    
    if not results:
        print("No successful drug analyses to compare")
        return
    
    print(f"\n📊 Drug Comparison Summary:")
    print("-" * 60)
    print(f"{'Drug':<15} {'Total Genes':<12} {'Significant':<12} {'GO Terms':<10}")
    print("-" * 60)
    
    for drug, result in results.items():
        print(f"{drug:<15} {result['total_genes_count']:<12} {result['significant_genes_count']:<12} {len(result['enrichment_results']):<10}")
    
    return results

def main():
    """Main function to run the drug analysis client"""
    print("🧬 Drug Gene Expression Analysis API Client")
    print("=" * 60)
    
    # Test API connection
    if not test_api_connection():
        return
    
    # Get available drugs
    available_drugs = get_available_drugs()
    if not available_drugs:
        return
    
    print(f"\nAvailable example drugs: {', '.join(available_drugs[:5])}...")
    
    # Analyze a single drug
    test_drug = "aspirin"  # Default test drug
    print(f"\n🎯 Analyzing single drug: {test_drug}")
    
    result = analyze_drug(test_drug)
    if result:
        # Display gene expression summary
        display_gene_expression_summary(result['gene_expression_data'])
        
        # Display enrichment results
        display_enrichment_results(result['enrichment_results'])
        
        # Save results to HTML
        save_drug_analysis_to_html(result, f"{test_drug}_analysis.html")
    
    # Compare multiple drugs
    print(f"\n🔄 Running drug comparison...")
    comparison_drugs = ["aspirin", "acetaminophen", "ibuprofen"]
    comparison_results = run_drug_comparison(comparison_drugs)
    
    print("\n🎉 Drug analysis complete!")
    print("You can also view the interactive interface at: http://localhost:5001")
    print("\nTo analyze a custom drug, use:")
    print("  curl -X POST http://localhost:5001/api/analyze_drug \\")
    print("       -H 'Content-Type: application/json' \\")
    print("       -d '{\"drug_name\": \"your_drug_name\"}'")

if __name__ == "__main__":
    main()
