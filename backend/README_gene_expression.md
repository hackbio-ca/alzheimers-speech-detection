# Drug Gene Expression Analysis Backend

A comprehensive Python Flask backend that analyzes drug effects on gene expression, generating volcano plots and GO enrichment plots. The system integrates with machine learning models to predict gene expression changes from drug molecular features.

## Features

- **Drug Analysis**: Predict gene expression changes from drug names using trained ML models
- **Volcano Plot Generation**: Creates interactive volcano plots from gene expression data (gene, LogFC, p-val)
- **GO Enrichment Analysis**: Performs Gene Ontology enrichment analysis with statistical testing
- **Interactive Plots**: Uses Plotly for interactive, web-ready visualizations
- **RESTful API**: Clean API endpoints for frontend integration
- **Built-in Test Interface**: HTML interface for testing the API
- **Drug Comparison**: Compare effects of multiple drugs on gene expression

## Installation

1. Install dependencies:
```bash
pip install -r requirements_gene_expression.txt
```

2. Run the server:
```bash
python gene_expression_server.py
```

The server will start on `http://localhost:5001`

## API Endpoints

### Drug Analysis Endpoints (NEW!)
- **POST** `/api/analyze_drug`
  - Analyzes drug effects on gene expression using ML predictions
  - Input: JSON with `drug_name` field
  - Output: Complete analysis with volcano plot, GO enrichment, and gene expression data
  - Example: `{"drug_name": "aspirin"}`

### Direct Gene Expression Analysis
- **POST** `/api/analyze_expression`
  - Performs complete analysis (volcano plot + GO enrichment)
  - Input: JSON with `genes` array containing `{gene, LogFC, p_val}` objects
  - Output: Both plots and enrichment results

### Individual Plot Endpoints
- **POST** `/api/volcano_plot` - Volcano plot only
- **POST** `/api/go_enrichment` - GO enrichment analysis only

### Utility Endpoints
- **GET** `/api/available_drugs` - Returns list of example drugs for analysis
- **GET** `/api/example_data` - Returns sample gene expression data for testing
- **GET** `/` - HTML interface for testing the API

## Input Data Format

### Drug Analysis Input
For drug analysis, simply provide the drug name:

```json
{
  "drug_name": "aspirin"
}
```

### Direct Gene Expression Input
For direct gene expression analysis:

```json
{
  "genes": [
    {
      "gene": "GENE_0001",
      "LogFC": 2.5,
      "p_val": 0.001
    },
    {
      "gene": "GENE_0002", 
      "LogFC": -1.8,
      "p_val": 0.05
    }
  ]
}
```

## Output Format

### Drug Analysis Output
The drug analysis endpoint returns:

```json
{
  "drug_name": "aspirin",
  "volcano_plot": { /* Plotly JSON */ },
  "go_enrichment_plot": { /* Plotly JSON */ },
  "enrichment_results": [ /* GO enrichment data */ ],
  "gene_expression_data": [ /* Full gene expression predictions */ ],
  "significant_genes_count": 50,
  "total_genes_count": 1000
}
```

### Direct Analysis Output
The direct gene expression analysis endpoint returns:

```json
{
  "volcano_plot": { /* Plotly JSON */ },
  "go_enrichment_plot": { /* Plotly JSON */ },
  "enrichment_results": [ /* GO enrichment data */ ],
  "significant_genes_count": 50,
  "total_genes_count": 1000
}
```

## Usage Examples

### Drug Analysis Example

```python
import requests
import json

# Analyze a drug
response = requests.post(
    'http://localhost:5001/api/analyze_drug',
    json={'drug_name': 'aspirin'},
    headers={'Content-Type': 'application/json'}
)

result = response.json()
print(f"Drug: {result['drug_name']}")
print(f"Found {result['significant_genes_count']} significantly affected genes")
```

### Direct Gene Expression Analysis Example

```python
import requests
import json

# Load example data
response = requests.get('http://localhost:5001/api/example_data')
data = response.json()

# Run full analysis
analysis_response = requests.post(
    'http://localhost:5001/api/analyze_expression',
    json=data,
    headers={'Content-Type': 'application/json'}
)

result = analysis_response.json()
print(f"Found {result['significant_genes_count']} significant genes")
```

### JavaScript/Frontend Integration

```javascript
// Analyze a drug
async function analyzeDrug(drugName) {
    const response = await fetch('/api/analyze_drug', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ drug_name: drugName })
    });
    
    const result = await response.json();
    
    // Display plots using Plotly
    Plotly.newPlot('volcano-plot', result.volcano_plot.data, result.volcano_plot.layout);
    Plotly.newPlot('go-plot', result.go_enrichment_plot.data, result.go_enrichment_plot.layout);
    
    console.log(`Analyzed ${result.drug_name}: ${result.significant_genes_count} significant genes`);
}

// Direct gene expression analysis
async function runDirectAnalysis() {
    const dataResponse = await fetch('/api/example_data');
    const data = await dataResponse.json();
    
    const analysisResponse = await fetch('/api/analyze_expression', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify(data)
    });
    
    const result = await analysisResponse.json();
    
    // Display plots using Plotly
    Plotly.newPlot('volcano-plot', result.volcano_plot.data, result.volcano_plot.layout);
    Plotly.newPlot('go-plot', result.go_enrichment_plot.data, result.go_enrichment_plot.layout);
}
```

## Plot Features

### Volcano Plot
- X-axis: Log2 Fold Change
- Y-axis: -log10(P-value)
- Color coding: Red (upregulated), Blue (downregulated), Gray (not significant)
- Significance thresholds: P < 0.05, |LogFC| > 1.0
- Interactive hover information

### GO Enrichment Plot
- Horizontal bar chart showing enriched GO terms
- X-axis: -log10(P-value)
- Color coding: Fold enrichment
- Top 20 most significant terms displayed
- Interactive hover with detailed statistics

## Technical Details

### Drug Prediction Pipeline
- **QSAR Features**: Extracts molecular descriptors using RDKit from PubChem SMILES
- **ML Model**: Uses trained neural network (ImprovedMLP) for gene expression prediction
- **Gene Coverage**: Predicts expression for all genes in training dataset
- **P-value Generation**: Mock p-values based on effect size (can be replaced with real statistical tests)

### Analysis Features
- **Statistical Testing**: Uses hypergeometric test for GO enrichment
- **Significance Thresholds**: P < 0.05, |LogFC| > 1.0 (configurable)
- **GO Database**: Simplified mock database (can be extended with real GO data)
- **Plot Rendering**: Plotly for interactive, web-ready plots
- **CORS Enabled**: Ready for frontend integration

## Extending the System

1. **Real GO Database**: Replace mock GO terms with actual Gene Ontology database
2. **Additional Plot Types**: Add MA plots, heatmaps, etc.
3. **Statistical Methods**: Add multiple testing correction, different enrichment methods
4. **Data Validation**: Enhanced input validation and error handling
5. **Caching**: Add caching for large datasets and repeated analyses
6. **Real P-values**: Replace mock p-values with actual statistical tests from experimental data
7. **Model Updates**: Integrate with updated ML models or different prediction algorithms
8. **Drug Database**: Expand drug database with more compounds and better name matching

## Dependencies

- Flask: Web framework
- Flask-CORS: Cross-origin resource sharing
- Pandas: Data manipulation
- NumPy: Numerical computing
- Plotly: Interactive plotting
- SciPy: Statistical functions
- Requests: HTTP library
- PyTorch: Machine learning framework for drug prediction
- Scikit-learn: Machine learning utilities
- RDKit: Chemical informatics and molecular descriptors
- PubChemPy: PubChem database access
