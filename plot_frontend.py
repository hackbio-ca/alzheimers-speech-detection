fetch("http://localhost:8000/predict/", {
  method: "POST",
  headers: { "Content-Type": "application/json" },
  body: JSON.stringify({
    drug_name: "aspirin",
    gene_list: ["TP53", "EGFR", "BRCA1", "MYC", "PTEN"]
  })
})
  .then(res => res.json())
  .then(data => {
    Plotly.newPlot("volcanoDiv", JSON.parse(data.volcano));
    Plotly.newPlot("barplotDiv", JSON.parse(data.barplot));
  });
