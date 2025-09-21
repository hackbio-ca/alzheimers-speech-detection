# backend_api.py
from fastapi import FastAPI
from pydantic import BaseModel
from typing import List
from tahoe_mlp_optimized import OptimizedTahoeMLPTrainer
from plot_generator import generate_plots

app = FastAPI()

# Load model once at startup
trainer = OptimizedTahoeMLPTrainer()
trainer.load_model("best_optimized_model.pth")  # path to your trained model

class PredictionRequest(BaseModel):
    drug_name: str
    gene_list: List[str]

@app.post("/predict/")
def predict(request: PredictionRequest):
    plots = generate_plots(trainer, request.drug_name, request.gene_list)
    return plots  # returns {"volcano": "...json...", "barplot": "...json..."}
