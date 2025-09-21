import React, { useState } from 'react';
import Plotly from 'plotly.js-dist';
import createPlotlyComponent from 'react-plotly.js';
import Footer from './Footer'; // Assuming they have this component

const Plot = createPlotlyComponent(Plotly);

function FindSideEffects() {
  const [input, setInput] = useState("");
  const [result, setResult] = useState("");
  const [loading, setLoading] = useState(false);
  const [analysisData, setAnalysisData] = useState(null);
  const [error, setError] = useState("");

  const handleSubmit = async (e) => {
    e.preventDefault();
    setLoading(true);
    setResult("");
    setAnalysisData(null);
    setError("");
    
    try {
      const res = await fetch("http://localhost:5001/api/find", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ query: input }),
      });
      
      const data = await res.json();
      
      if (data.error) {
        setError(data.error);
        setResult(data.result || `Error: ${data.error}`);
      } else if (data.is_drug_analysis) {
        // This is a drug analysis result
        setAnalysisData(data);
        setResult(data.result);
      } else {
        // This is a regular search result
        setResult(data.result);
      }
    } catch (err) {
      console.error(err);
      setError("Network error - please check if the server is running");
      setResult("Error fetching results");
    } finally {
      setLoading(false);
    }
  };

  const renderGeneTable = () => {
    if (!analysisData || !analysisData.gene_expression_data) return null;
    
    const topGenes = analysisData.gene_expression_data
      .sort((a, b) => Math.abs(b.LogFC) - Math.abs(a.LogFC))
      .slice(0, 15); // Show top 15 genes
    
    return (
      <div className="gene-table-container">
        <h3>Top Gene Expression Changes</h3>
        <div className="table-wrapper">
          <table className="gene-table">
            <thead>
              <tr>
                <th>Rank</th>
                <th>Gene</th>
                <th>Log2 Fold Change</th>
                <th>P-value</th>
                <th>Effect</th>
              </tr>
            </thead>
            <tbody>
              {topGenes.map((gene, index) => (
                <tr key={index}>
                  <td>{index + 1}</td>
                  <td className="gene-name">{gene.gene}</td>
                  <td className="logfc">{gene.LogFC.toFixed(4)}</td>
                  <td className="pvalue">{gene.p_val.toExponential(2)}</td>
                  <td className={`effect ${getEffectClass(gene.LogFC)}`}>
                    {getEffectDescription(gene.LogFC)}
                  </td>
                </tr>
              ))}
            </tbody>
          </table>
        </div>
      </div>
    );
  };

  const renderGOEnrichment = () => {
    if (!analysisData || !analysisData.enrichment_results || analysisData.enrichment_results.length === 0) {
      return null;
    }
    
    const topGO = analysisData.enrichment_results.slice(0, 10);
    
    return (
      <div className="go-enrichment-container">
        <h3>Top GO Enrichment Results</h3>
        <div className="go-table-wrapper">
          <table className="go-table">
            <thead>
              <tr>
                <th>GO Term</th>
                <th>P-value</th>
                <th>Fold Enrichment</th>
                <th>Genes</th>
              </tr>
            </thead>
            <tbody>
              {topGO.map((go, index) => (
                <tr key={index}>
                  <td className="go-term">{go.go_name}</td>
                  <td className="pvalue">{go.p_value.toExponential(2)}</td>
                  <td className="fold-enrichment">{go.fold_enrichment.toFixed(2)}</td>
                  <td className="gene-count">{go.overlap}</td>
                </tr>
              ))}
            </tbody>
          </table>
        </div>
      </div>
    );
  };

  const getEffectDescription = (log2fc) => {
    if (log2fc > 2.0) return "Strong Upregulation";
    if (log2fc > 1.0) return "Moderate Upregulation";
    if (log2fc > 0.5) return "Weak Upregulation";
    if (log2fc > -0.5) return "No Change";
    if (log2fc > -1.0) return "Weak Downregulation";
    if (log2fc > -2.0) return "Moderate Downregulation";
    return "Strong Downregulation";
  };

  const getEffectClass = (log2fc) => {
    if (log2fc > 1.0) return "upregulated";
    if (log2fc < -1.0) return "downregulated";
    return "no-change";
  };

  const getSignificanceClass = (pval) => {
    if (pval < 0.001) return "highly-significant";
    if (pval < 0.01) return "very-significant";
    if (pval < 0.05) return "significant";
    return "not-significant";
  };

  return (
    <>
      <div id="findSideEffects">
        <form id="sideEffectsForm" onSubmit={handleSubmit}>
          <h1>Find Drug Effects</h1>
          <p className="subtitle">Enter a drug name to analyze its effects on gene expression</p>
          
          <div className="input-container">
            <input
              type="text"
              placeholder="Enter drug candidate (e.g. Aspirin, Ibuprofen, Metformin)"
              value={input}
              onChange={(e) => setInput(e.target.value)}
              required
              className="drug-input"
            />
            <button 
              id="submission" 
              type="submit" 
              disabled={loading}
              className="submit-button"
            >
              {loading ? "Analyzing..." : "Analyze Drug"}
            </button>
          </div>
          
          {/* Example drugs */}
          <div className="example-drugs">
            <p>Try these examples:</p>
            <div className="drug-tags">
              {['aspirin', 'ibuprofen', 'metformin', 'caffeine', 'morphine'].map(drug => (
                <span 
                  key={drug}
                  className="drug-tag" 
                  onClick={() => setInput(drug)}
                >
                  {drug}
                </span>
              ))}
            </div>
          </div>
          
          {/* Results */}
          {error && (
            <div className="error-message">
              <h4>❌ Error</h4>
              <p>{error}</p>
            </div>
          )}
          
          {result && (
            <div className="result">
              <p>{result}</p>
            </div>
          )}
          
          {/* Drug Analysis Results */}
          {analysisData && (
            <div className="analysis-results">
              <h2>🧬 Analysis Results for: {analysisData.drug_name}</h2>
              
              <div className="stats-grid">
                <div className="stat-card">
                  <h4>Total Genes</h4>
                  <p className="stat-number">{analysisData.total_genes_count}</p>
                </div>
                <div className="stat-card">
                  <h4>Significant Genes</h4>
                  <p className="stat-number">{analysisData.significant_genes_count}</p>
                </div>
                <div className="stat-card">
                  <h4>GO Terms</h4>
                  <p className="stat-number">{analysisData.enrichment_results.length}</p>
                </div>
                <div className="stat-card">
                  <h4>Success Rate</h4>
                  <p className="stat-number">
                    {((analysisData.significant_genes_count / analysisData.total_genes_count) * 100).toFixed(1)}%
                  </p>
                </div>
              </div>
              
              {/* Volcano Plot */}
              {analysisData.volcano_plot && (
                <div className="plot-container">
                  <h3>Volcano Plot - Gene Expression Changes</h3>
                  <div className="plot-wrapper">
                    <Plot
                      data={analysisData.volcano_plot.data}
                      layout={{
                        ...analysisData.volcano_plot.layout,
                        width: 800,
                        height: 600,
                        responsive: true
                      }}
                      config={{
                        displayModeBar: true,
                        displaylogo: false,
                        modeBarButtonsToRemove: ['pan2d', 'lasso2d', 'select2d']
                      }}
                    />
                  </div>
                </div>
              )}
              
              {/* GO Enrichment Plot */}
              {analysisData.go_enrichment_plot && analysisData.enrichment_results.length > 0 && (
                <div className="plot-container">
                  <h3>GO Enrichment Analysis</h3>
                  <div className="plot-wrapper">
                    <Plot
                      data={analysisData.go_enrichment_plot.data}
                      layout={{
                        ...analysisData.go_enrichment_plot.layout,
                        width: 800,
                        height: 600,
                        responsive: true
                      }}
                      config={{
                        displayModeBar: true,
                        displaylogo: false,
                        modeBarButtonsToRemove: ['pan2d', 'lasso2d', 'select2d']
                      }}
                    />
                  </div>
                </div>
              )}
              
              {/* Gene Expression Table */}
              {renderGeneTable()}
              
              {/* GO Enrichment Table */}
              {renderGOEnrichment()}
            </div>
          )}
        </form>
      </div>
      <Footer />
    </>
  );
}

export default FindSideEffects;
