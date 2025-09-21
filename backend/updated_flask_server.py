from flask import Flask, request, jsonify
from flask_cors import CORS
import sys
import os

# Add parent directory to path to import the drug analysis functions
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Import the drug analysis functions
try:
    from drug_gene_expression_server import (
        predict_drug_gene_expression, 
        create_volcano_plot, 
        create_go_enrichment_plot,
        perform_go_enrichment,
        find_model_file,
        DRUG_PREDICTION_AVAILABLE
    )
    print("✅ Drug analysis modules imported successfully")
except ImportError as e:
    print(f"❌ Error importing drug analysis modules: {e}")
    DRUG_PREDICTION_AVAILABLE = False

app = Flask(__name__)
CORS(app)

# Helper function to determine if input is a drug name
def is_drug_name(input_text):
    """
    Determine if the input text looks like a drug name
    """
    drug_keywords = [
        'aspirin', 'acetaminophen', 'ibuprofen', 'metformin', 'warfarin',
        'simvastatin', 'lisinopril', 'amlodipine', 'omeprazole', 'prednisone',
        'caffeine', 'morphine', 'diazepam', 'fluoxetine', 'sertraline',
        'paracetamol', 'tylenol', 'advil', 'motrin', 'aleve', 'naproxen',
        'glucophage', 'coumadin', 'zocor', 'prinivil', 'norvasc', 'prilosec',
        'deltasone', 'xanax', 'valium', 'prozac', 'zoloft'
    ]
    
    lower_input = input_text.lower().strip()
    
    # Check if it's a known drug name
    if lower_input in drug_keywords:
        return True
    
    # Check if it looks like a drug name (no spaces, mostly letters, reasonable length)
    if (len(lower_input) > 2 and len(lower_input) < 30 and 
        lower_input.replace(' ', '').replace('-', '').isalpha() and 
        ' ' not in lower_input and lower_input.count('-') <= 1):
        return True
    
    return False

# Your friend's original endpoint - now enhanced with drug analysis
@app.route("/api/find", methods=["POST"])
def find_side_effects():
    """
    Enhanced endpoint that handles both drug analysis and general search
    """
    try:
        data = request.get_json()
        query = data.get("query", "").strip()
        
        if not query:
            return jsonify({"error": "No query provided"}), 400
        
        # Check if the query looks like a drug name
        if is_drug_name(query):
            print(f"🔬 Detected drug name: {query}")
            # Use the drug analysis functionality
            return analyze_drug_internal(query)
        else:
            print(f"🔍 General search query: {query}")
            # Use the original functionality
            return jsonify({"result": f"You searched for: {query}"})
            
    except Exception as e:
        print(f"❌ Error in find_side_effects: {e}")
        return jsonify({"error": str(e)}), 500

# Internal drug analysis function
def analyze_drug_internal(drug_name):
    """
    Internal function to analyze drug effects on gene expression
    """
    if not DRUG_PREDICTION_AVAILABLE:
        return jsonify({
            'error': 'Drug prediction modules not available. Please ensure all required files are in place.',
            'result': f'Drug analysis not available for: {drug_name}'
        }), 500
    
    try:
        print(f"🧬 Starting drug analysis for: {drug_name}")
        
        # Use the drug analysis functionality
        gene_expression_data, actual_drug_name = predict_drug_gene_expression(drug_name)
        
        # Create volcano plot
        volcano_fig = create_volcano_plot(gene_expression_data)
        
        # Get significant genes for GO enrichment
        significant_genes = [
            item['gene'] for item in gene_expression_data 
            if item['p_val'] < 0.05 and abs(item['LogFC']) > 1.0
        ]
        
        # Perform GO enrichment
        enrichment_results = perform_go_enrichment(significant_genes)
        
        # Create GO enrichment plot
        go_fig = create_go_enrichment_plot(enrichment_results)
        
        # Prepare response with both original format and new data
        response = {
            # Original format for backward compatibility
            "result": f"Drug analysis completed for {actual_drug_name}: {len(gene_expression_data)} genes analyzed, {len(significant_genes)} significantly affected",
            
            # New drug analysis data
            "drug_name": actual_drug_name,
            "volcano_plot": volcano_fig.to_dict(),
            "go_enrichment_plot": go_fig.to_dict(),
            "enrichment_results": enrichment_results,
            "gene_expression_data": gene_expression_data,
            "significant_genes_count": len(significant_genes),
            "total_genes_count": len(gene_expression_data),
            "is_drug_analysis": True
        }
        
        print(f"✅ Drug analysis completed for {actual_drug_name}")
        return jsonify(response)
        
    except Exception as e:
        print(f"❌ Error in drug analysis: {e}")
        return jsonify({
            "error": str(e),
            "result": f"Error analyzing drug {drug_name}: {str(e)}"
        }), 500

# Direct drug analysis endpoint (alternative)
@app.route("/api/analyze_drug", methods=["POST"])
def analyze_drug_direct():
    """
    Direct endpoint for drug analysis (alternative to /api/find)
    """
    try:
        data = request.get_json()
        drug_name = data.get('drug_name', '').strip()
        
        if not drug_name:
            return jsonify({'error': 'No drug name provided'}), 400
        
        return analyze_drug_internal(drug_name)
        
    except Exception as e:
        return jsonify({'error': str(e)}), 500

# Health check endpoint
@app.route("/api/health")
def health_check():
    """
    Check if the server and drug analysis are working
    """
    try:
        model_path = find_model_file()
        model_exists = os.path.exists(model_path)
        
        return jsonify({
            'status': 'healthy',
            'drug_analysis_available': DRUG_PREDICTION_AVAILABLE,
            'model_exists': model_exists,
            'model_path': os.path.abspath(model_path),
            'endpoints': [
                '/api/find - Enhanced endpoint (drug analysis + general search)',
                '/api/analyze_drug - Direct drug analysis',
                '/api/health - This endpoint'
            ]
        })
    except Exception as e:
        return jsonify({
            'status': 'error',
            'error': str(e)
        }), 500

# Available drugs endpoint
@app.route("/api/available_drugs")
def available_drugs():
    """
    Return a list of example drugs that can be analyzed
    """
    example_drugs = [
        "aspirin", "acetaminophen", "ibuprofen", "metformin", "warfarin",
        "simvastatin", "lisinopril", "amlodipine", "omeprazole", "prednisone",
        "caffeine", "morphine", "diazepam", "fluoxetine", "sertraline",
        "paracetamol", "tylenol", "advil", "motrin", "aleve"
    ]
    
    return jsonify({
        'drugs': example_drugs,
        'note': 'These are example drugs. The system can analyze any drug that has PubChem data.'
    })

if __name__ == "__main__":
    print("🚀 Starting enhanced Flask server...")
    print("📊 Available endpoints:")
    print("   - POST /api/find - Enhanced endpoint (drug analysis + general search)")
    print("   - POST /api/analyze_drug - Direct drug analysis")
    print("   - GET /api/health - Health check")
    print("   - GET /api/available_drugs - Example drugs")
    print(f"🧬 Drug prediction available: {DRUG_PREDICTION_AVAILABLE}")
    
    app.run(host="0.0.0.0", port=5001, debug=True)
