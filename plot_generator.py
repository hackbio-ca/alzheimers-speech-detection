# plot_generator.py
import pandas as pd
import numpy as np
import plotly.express as px
from tahoe_mlp_optimized import OptimizedTahoeMLPTrainer

### Step 1: Generate Predictions DataFrame
def generate_predictions(trainer, drug_name: str, gene_list: list) -> pd.DataFrame:
    """
    Uses trained model to generate predictions for a drug across multiple genes.
    Returns DataFrame with columns: gene, logFC, pval
    """
    predictions = []
    rng = np.random.default_rng(42)  # temporary mock p-values

    for gene in gene_list:
        logFC = trainer.predict(drug_name, gene)
        predictions.append({
            "gene": gene,
            "logFC": logFC,
            "pval": rng.uniform(0.0001, 0.5)  # mock until model outputs real stats
        })

    return pd.DataFrame(predictions)

### Step 2a: Volcano Plot
def make_volcano_plot(pred_df: pd.DataFrame):
    pred_df = pred_df.copy()
    pred_df["-log10p"] = -np.log10(pred_df["pval"])

    fig = px.scatter(
        pred_df,
        x="logFC",
        y="-log10p",
        text="gene",
        title="Predicted Drug-Gene Volcano Plot",
        hover_data=["gene", "logFC", "pval"],
        color="logFC",
        color_continuous_scale="RdBu",
    )
    fig.update_traces(marker=dict(size=10, opacity=0.7, line=dict(width=1, color="DarkSlateGrey")))
    fig.update_layout(template="plotly_white")
    return fig.to_json()

### Step 2b: Ranked Barplot
def make_ranked_barplot(pred_df: pd.DataFrame, top_n=20):
    df_top = pred_df.reindex(pred_df["logFC"].abs().sort_values(ascending=False).index).head(top_n)

    fig = px.bar(
        df_top,
        x="gene",
        y="logFC",
        title=f"Top {top_n} Predicted Differential Genes",
        color="logFC",
        color_continuous_scale="RdBu"
    )
    fig.update_layout(template="plotly_white", xaxis_tickangle=-45)
    return fig.to_json()

### Step 3: Wrap both plots
def generate_plots(trainer, drug_name: str, gene_list: list):
    pred_df = generate_predictions(trainer, drug_name, gene_list)
    return {
        "volcano": make_volcano_plot(pred_df),
        "barplot": make_ranked_barplot(pred_df)
    }
