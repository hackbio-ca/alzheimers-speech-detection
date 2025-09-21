# Installation Guide for Your Friend

## Files to Copy

Copy these files to your friend's project:

### Backend Files:
1. `updated_flask_server.py` → Replace their existing Flask server
2. `drug_gene_expression_server.py` → Place in backend/ directory
3. `requirements_gene_expression.txt` → Place in backend/ directory

### ML Model Files:
4. `best_fast_model.pth` → Place in project root directory
5. `predict_all_genes.py` → Place in project root directory  
6. `standalone_drug_descriptors.py` → Place in project root directory

### Frontend Files:
7. `updated_react_component.jsx` → Replace their FindSideEffects component
8. `updated_styles.css` → Add to their CSS file

## Step-by-Step Installation

### 1. Backend Setup

```bash
# Install Python dependencies
pip install -r backend/requirements_gene_expression.txt

# Replace their server.py with updated_flask_server.py
# Make sure the file structure looks like this:
your_project/
├── backend/
│   ├── updated_flask_server.py (rename from their server.py)
│   ├── drug_gene_expression_server.py
│   └── requirements_gene_expression.txt
├── best_fast_model.pth
├── predict_all_genes.py
└── standalone_drug_descriptors.py
```

### 2. Frontend Setup

```bash
# Install React dependencies
npm install plotly.js-dist react-plotly.js

# Replace their FindSideEffects component with updated_react_component.jsx
# Add the CSS styles from updated_styles.css to their CSS file
```

### 3. File Structure

Your friend's final project structure should look like:

```
their_project/
├── backend/
│   ├── updated_flask_server.py (their main server)
│   ├── drug_gene_expression_server.py
│   └── requirements_gene_expression.txt
├── best_fast_model.pth
├── predict_all_genes.py
├── standalone_drug_descriptors.py
└── frontend/
    ├── src/
    │   └── components/
    │       └── FindSideEffects.jsx (updated component)
    └── styles/
        └── styles.css (with added styles)
```

## How to Run

### 1. Start the Backend
```bash
cd backend
python updated_flask_server.py
```

### 2. Start the Frontend
```bash
cd frontend
npm start
```

## What Your Friend Gets

### Enhanced API Endpoint
- **Same endpoint**: `/api/find` (no frontend changes needed for requests)
- **Automatic detection**: Backend detects if input is a drug name
- **Drug analysis**: Returns gene expression predictions and plots
- **Backward compatibility**: Non-drug queries still work

### New Features
1. **Interactive Volcano Plots**: Shows gene expression changes
2. **GO Enrichment Analysis**: Biological pathway analysis
3. **Gene Expression Table**: Top affected genes with statistics
4. **Responsive Design**: Works on mobile and desktop
5. **Error Handling**: Clear error messages and loading states

### Example Usage
- User enters "aspirin" → Gets full drug analysis with plots
- User enters "side effects" → Gets original search functionality
- User enters "ibuprofen" → Gets drug analysis
- User enters "what is cancer" → Gets original search

## Testing

### Test the Backend
```bash
# Test health check
curl http://localhost:5001/api/health

# Test drug analysis
curl -X POST http://localhost:5001/api/find \
  -H "Content-Type: application/json" \
  -d '{"query": "aspirin"}'

# Test general search
curl -X POST http://localhost:5001/api/find \
  -H "Content-Type: application/json" \
  -d '{"query": "hello world"}'
```

### Test the Frontend
1. Open the React app
2. Enter "aspirin" in the search box
3. Click "Analyze Drug"
4. Should see volcano plot, GO enrichment, and gene table

## Troubleshooting

### Model File Not Found
```bash
# Check if model file exists
ls -la best_fast_model.pth

# Should show: -rw-r--r-- 1 user user 1459206 date best_fast_model.pth
```

### Import Errors
```bash
# Install missing dependencies
pip install torch scikit-learn rdkit-pypi pubchempy

# Check Python path
python -c "import sys; print(sys.path)"
```

### Frontend Errors
```bash
# Install Plotly
npm install plotly.js-dist react-plotly.js

# Check if components are imported correctly
```

## API Response Format

### Drug Analysis Response
```json
{
  "result": "Drug analysis completed for aspirin: 1000 genes analyzed, 50 significantly affected",
  "drug_name": "aspirin",
  "volcano_plot": { /* Plotly data */ },
  "go_enrichment_plot": { /* Plotly data */ },
  "gene_expression_data": [ /* Array of gene data */ ],
  "significant_genes_count": 50,
  "total_genes_count": 1000,
  "is_drug_analysis": true
}
```

### General Search Response
```json
{
  "result": "You searched for: hello world"
}
```

## Support

If your friend has issues:
1. Check the console for error messages
2. Verify all files are in the correct locations
3. Ensure all dependencies are installed
4. Test the backend endpoints directly with curl
5. Check the browser developer tools for frontend errors
