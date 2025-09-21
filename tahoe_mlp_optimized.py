#!/usr/bin/env python3
"""
Optimized Tahoe-100M Drug-Gene Expression Prediction Model
Enhanced version with better performance, larger batch sizes, and advanced techniques.
"""

import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import Dataset, DataLoader
from datasets import load_dataset
import warnings
from typing import Dict, List, Tuple, Optional
import logging
from pathlib import Path
import pickle
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error, r2_score
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter
import gc

warnings.filterwarnings('ignore')

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class OptimizedTahoeDrugDataset(Dataset):
    """
    Optimized PyTorch Dataset with better memory management and feature engineering
    """
    
    def __init__(self, data: pd.DataFrame, drug_features: Dict[str, np.ndarray], 
                 gene_encoder: LabelEncoder, gene_features: Optional[Dict[str, np.ndarray]] = None,
                 use_gene_embeddings: bool = True, gene_embedding_dim: int = 128):
        """
        Initialize optimized dataset
        
        Args:
            data: DataFrame with columns ['drug', 'gene_name', 'log2FoldChange', ...]
            drug_features: Dict mapping drug names to feature vectors
            gene_encoder: LabelEncoder for gene names
            gene_features: Optional dict mapping gene names to feature vectors
            use_gene_embeddings: Whether to use learned gene embeddings
            gene_embedding_dim: Dimension of gene embeddings
        """
        self.data = data.reset_index(drop=True)
        self.drug_features = drug_features
        self.gene_encoder = gene_encoder
        self.gene_features = gene_features
        self.use_gene_embeddings = use_gene_embeddings
        self.gene_embedding_dim = gene_embedding_dim
        
        # Pre-compute valid indices and gene indices
        self.valid_indices = []
        self.gene_indices = []
        
        for idx, row in self.data.iterrows():
            if row['drug'] in self.drug_features:
                self.valid_indices.append(idx)
                if use_gene_embeddings:
                    try:
                        gene_idx = self.gene_encoder.transform([row['gene_name']])[0]
                        self.gene_indices.append(gene_idx)
                    except ValueError:
                        # Handle unseen genes
                        self.gene_indices.append(0)  # Default to first gene
                else:
                    self.gene_indices.append(0)  # Placeholder
        
        # Pre-compute drug feature tensors for faster access
        self.drug_feature_tensors = {}
        for drug, features in self.drug_features.items():
            self.drug_feature_tensors[drug] = torch.FloatTensor(features)
        
        logger.info(f"Optimized dataset initialized with {len(self.valid_indices)} valid samples")
        logger.info(f"Gene embedding dimension: {gene_embedding_dim if use_gene_embeddings else 'N/A'}")
    
    def __len__(self):
        return len(self.valid_indices)
    
    def __getitem__(self, idx):
        actual_idx = self.valid_indices[idx]
        row = self.data.iloc[actual_idx]
        
        # Get pre-computed drug features
        drug_feat = self.drug_feature_tensors[row['drug']]
        
        if self.use_gene_embeddings:
            gene_idx = self.gene_indices[idx]
            return drug_feat, torch.LongTensor([gene_idx]), torch.FloatTensor([float(row['log2FoldChange'])])
        else:
            # Use one-hot encoding for genes
            gene_feat = torch.zeros(1000)
            gene_hash = hash(row['gene_name']) % 1000
            gene_feat[gene_hash] = 1.0
            features = torch.cat([drug_feat, gene_feat])
            return features, torch.FloatTensor([float(row['log2FoldChange'])])

class GeneEmbeddingMLP(nn.Module):
    """
    MLP with learned gene embeddings for better gene representation
    """
    
    def __init__(self, drug_dim: int, num_genes: int, gene_embedding_dim: int = 128,
                 hidden_dims: List[int] = [512, 256, 128], dropout: float = 0.3):
        super(GeneEmbeddingMLP, self).__init__()
        
        # Gene embedding layer
        self.gene_embedding = nn.Embedding(num_genes, gene_embedding_dim)
        
        # Drug feature processing
        self.drug_processor = nn.Sequential(
            nn.Linear(drug_dim, 256),
            nn.BatchNorm1d(256),
            nn.ReLU(),
            nn.Dropout(dropout)
        )
        
        # Combined feature processing
        combined_dim = 256 + gene_embedding_dim
        layers = []
        prev_dim = combined_dim
        
        for hidden_dim in hidden_dims:
            layers.extend([
                nn.Linear(prev_dim, hidden_dim),
                nn.BatchNorm1d(hidden_dim),
                nn.ReLU(),
                nn.Dropout(dropout)
            ])
            prev_dim = hidden_dim
        
        self.combined_network = nn.Sequential(*layers)
        self.output_layer = nn.Linear(prev_dim, 1)
        
    def forward(self, drug_features, gene_indices):
        # Process drug features
        drug_processed = self.drug_processor(drug_features)
        
        # Get gene embeddings
        gene_embeddings = self.gene_embedding(gene_indices).squeeze(1)
        
        # Combine features
        combined = torch.cat([drug_processed, gene_embeddings], dim=1)
        
        # Final prediction
        output = self.combined_network(combined)
        return self.output_layer(output).squeeze(-1)

class ResidualMLP(nn.Module):
    """
    MLP with residual connections for better gradient flow
    """
    
    def __init__(self, input_dim: int, hidden_dims: List[int] = [512, 512, 256, 128], 
                 dropout: float = 0.3):
        super(ResidualMLP, self).__init__()
        
        self.input_dim = input_dim
        self.hidden_dims = hidden_dims
        
        # Input projection
        self.input_proj = nn.Linear(input_dim, hidden_dims[0])
        
        # Residual blocks
        self.residual_blocks = nn.ModuleList()
        for i in range(len(hidden_dims) - 1):
            block = nn.Sequential(
                nn.Linear(hidden_dims[i], hidden_dims[i+1]),
                nn.BatchNorm1d(hidden_dims[i+1]),
                nn.ReLU(),
                nn.Dropout(dropout)
            )
            self.residual_blocks.append(block)
        
        # Output layer
        self.output_layer = nn.Linear(hidden_dims[-1], 1)
        
    def forward(self, x):
        x = self.input_proj(x)
        
        for block in self.residual_blocks:
            residual = x
            x = block(x)
            # Add residual connection if dimensions match
            if x.shape == residual.shape:
                x = x + residual
        
        return self.output_layer(x).squeeze(-1)

class AttentionFusionModel(nn.Module):
    """
    Advanced fusion model with attention mechanism
    """
    
    def __init__(self, drug_dim: int, gene_dim: int, hidden_dim: int = 256, 
                 num_heads: int = 8, dropout: float = 0.2):
        super(AttentionFusionModel, self).__init__()
        
        # Drug encoder
        self.drug_encoder = nn.Sequential(
            nn.Linear(drug_dim, hidden_dim),
            nn.BatchNorm1d(hidden_dim),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim, hidden_dim)
        )
        
        # Gene encoder
        self.gene_encoder = nn.Sequential(
            nn.Linear(gene_dim, hidden_dim),
            nn.BatchNorm1d(hidden_dim),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim, hidden_dim)
        )
        
        # Multi-head attention
        self.attention = nn.MultiheadAttention(hidden_dim, num_heads, dropout=dropout, batch_first=True)
        
        # Fusion layers
        self.fusion = nn.Sequential(
            nn.Linear(hidden_dim * 2, hidden_dim),
            nn.BatchNorm1d(hidden_dim),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim, hidden_dim // 2),
            nn.BatchNorm1d(hidden_dim // 2),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim // 2, 1)
        )
        
    def forward(self, drug_features, gene_features):
        # Encode features
        drug_encoded = self.drug_encoder(drug_features).unsqueeze(1)  # [batch, 1, hidden]
        gene_encoded = self.gene_encoder(gene_features).unsqueeze(1)  # [batch, 1, hidden]
        
        # Concatenate for attention
        combined = torch.cat([drug_encoded, gene_encoded], dim=1)  # [batch, 2, hidden]
        
        # Apply attention
        attended, _ = self.attention(combined, combined, combined)
        
        # Flatten and fuse
        fused = attended.view(attended.size(0), -1)  # [batch, 2*hidden]
        output = self.fusion(fused)
        
        return output.squeeze(-1)

class OptimizedTahoeMLPTrainer:
    """
    Optimized trainer with better performance and memory management
    """
    
    def __init__(self, model_type: str = 'gene_embedding', device: str = 'auto'):
        self.model_type = model_type
        self.device = torch.device('cuda' if torch.cuda.is_available() and device == 'auto' else 'cpu')
        logger.info(f"Using device: {self.device}")
        
        self.model = None
        self.scaler = StandardScaler()
        self.drug_features = {}
        self.gene_encoder = None
        self.gene_features = {}
        
    def load_drug_descriptors(self, csv_path: str) -> Tuple[Dict[str, np.ndarray], List[str]]:
        """Load and normalize drug descriptors"""
        logger.info(f"Loading drug descriptors from {csv_path}")
        df = pd.read_csv(csv_path)
        
        drug_features = {}
        feature_columns = [col for col in df.columns if col not in ['drug_name', 'smiles']]
        
        # Clean and validate features
        feature_matrix = df[feature_columns].values
        
        # Check for infinite values and replace with NaN
        feature_matrix = np.where(np.isfinite(feature_matrix), feature_matrix, np.nan)
        
        # Check for extremely large values and clip them
        feature_matrix = np.clip(feature_matrix, -1e10, 1e10)
        
        # Replace any remaining NaN values with median of each column
        for col_idx in range(feature_matrix.shape[1]):
            col_data = feature_matrix[:, col_idx]
            if np.any(np.isnan(col_data)):
                median_val = np.nanmedian(col_data)
                feature_matrix[:, col_idx] = np.where(np.isnan(col_data), median_val, col_data)
                logger.warning(f"Replaced NaN values in column {feature_columns[col_idx]} with median: {median_val}")
        
        # Final check for any remaining issues
        if not np.all(np.isfinite(feature_matrix)):
            logger.error("Still contains infinite values after cleaning")
            raise ValueError("Feature matrix contains infinite values that could not be cleaned")
        
        logger.info(f"Feature matrix shape: {feature_matrix.shape}")
        logger.info(f"Feature matrix range: [{np.min(feature_matrix):.3f}, {np.max(feature_matrix):.3f}]")
        
        # Normalize features
        normalized_features = self.scaler.fit_transform(feature_matrix)
        
        for i, row in df.iterrows():
            drug_name = row['drug_name']
            features = normalized_features[i].astype(np.float32)
            drug_features[drug_name] = features
        
        logger.info(f"Loaded and normalized features for {len(drug_features)} drugs")
        logger.info(f"Feature dimension: {len(feature_columns)}")
        
        return drug_features, feature_columns
    
    def stream_tahoe_data(self, max_samples: int = 1000000, batch_size: int = 1000,
                         use_stratified: bool = True, min_cells: int = 1, 
                         min_genes: int = 1000, min_drugs: int = 30) -> pd.DataFrame:
        """Fast streaming with lightweight diversity sampling"""
        if use_stratified:
            logger.info(f"Streaming Tahoe-100M data with fast diversity sampling")
            logger.info(f"Max samples: {max_samples}, Min cells: {min_cells}, Min genes: {min_genes}, Min drugs: {min_drugs}")
            
            # Fast diversity sampling - much more efficient
            return self._fast_diversity_sampling(max_samples, min_cells, min_genes, min_drugs)
        
        else:
            # Original simple streaming
            logger.info(f"Streaming Tahoe-100M data (max {max_samples} samples)")
            
            dataset = load_dataset("tahoebio/Tahoe-100M", 
                                  "pseudobulk_differential_expression", 
                                  streaming=True, 
                                  split="train")
            
            data_rows = []
            sample_count = 0
            
            for sample in dataset:
                if sample_count >= max_samples:
                    break
                    
                if pd.notna(sample['log2FoldChange']):
                    data_rows.append({
                        'drug': sample['drug'],
                        'gene_name': sample['gene_name'],
                        'log2FoldChange': sample['log2FoldChange'],
                        'baseMean': sample['baseMean'],
                        'lfcSE': sample['lfcSE'],
                        'pvalue': sample['pvalue'],
                        'padj': sample['padj'],
                        'concentration': sample['concentration'],
                        'Cell_Name_Vevo': sample['Cell_Name_Vevo']
                    })
                    sample_count += 1
                
                # Memory management
                if sample_count % 50000 == 0:
                    logger.info(f"Processed {sample_count} samples")
                    gc.collect()  # Force garbage collection
            
            df = pd.DataFrame(data_rows)
            logger.info(f"Streamed {len(df)} valid samples")
            
            return df
    
    def _fast_diversity_sampling(self, max_samples: int, min_cells: int, min_genes: int, min_drugs: int) -> pd.DataFrame:
        """Fast diversity sampling with early stopping"""
        logger.info("Using fast diversity sampling strategy")
        
        dataset = load_dataset("tahoebio/Tahoe-100M", 
                              "pseudobulk_differential_expression", 
                              streaming=True, 
                              split="train")
        
        data_rows = []
        sample_count = 0
        processed_count = 0
        
        # Track diversity
        unique_cells = set()
        unique_genes = set()
        unique_drugs = set()
        
        # Sampling strategy: collect samples in batches and then filter for diversity
        batch_size = 10000
        batch_data = []
        
        for sample in dataset:
            processed_count += 1
            
            # Skip samples with missing log2FoldChange
            if pd.isna(sample.get('log2FoldChange')):
                continue
            
            # Add to batch
            batch_data.append({
                'drug': sample['drug'],
                'gene_name': sample['gene_name'],
                'log2FoldChange': sample['log2FoldChange'],
                'baseMean': sample['baseMean'],
                'lfcSE': sample['lfcSE'],
                'pvalue': sample['pvalue'],
                'padj': sample['padj'],
                'concentration': sample['concentration'],
                'Cell_Name_Vevo': sample['Cell_Name_Vevo']
            })
            
            # Process batch when full
            if len(batch_data) >= batch_size:
                # Add all samples from batch (fast approach)
                data_rows.extend(batch_data)
                sample_count += len(batch_data)
                
                # Update diversity tracking
                for item in batch_data:
                    unique_cells.add(item['Cell_Name_Vevo'])
                    unique_genes.add(item['gene_name'])
                    unique_drugs.add(item['drug'])
                
                batch_data = []
                
                # Log progress
                if sample_count % 50000 == 0:
                    logger.info(f"Collected {sample_count:,} samples, processed {processed_count:,} total")
                    logger.info(f"  Unique cells: {len(unique_cells)}, genes: {len(unique_genes)}, drugs: {len(unique_drugs)}")
                
                # Early stopping if we have enough samples and diversity
                if (sample_count >= max_samples and 
                    len(unique_cells) >= min_cells and 
                    len(unique_genes) >= min_genes and 
                    len(unique_drugs) >= min_drugs):
                    logger.info("Diversity requirements met, stopping early")
                    break
                
                # Memory management
                gc.collect()
        
        # Add remaining batch data
        if batch_data:
            data_rows.extend(batch_data)
            sample_count += len(batch_data)
            for item in batch_data:
                unique_cells.add(item['Cell_Name_Vevo'])
                unique_genes.add(item['gene_name'])
                unique_drugs.add(item['drug'])
        
        df = pd.DataFrame(data_rows)
        
        # Final diversity check
        logger.info(f"\nFinal sampling results:")
        logger.info(f"  Total samples: {len(df):,}")
        logger.info(f"  Unique cells: {len(unique_cells)}")
        logger.info(f"  Unique genes: {len(unique_genes)}")
        logger.info(f"  Unique drugs: {len(unique_drugs)}")
        logger.info(f"  Total processed: {processed_count:,}")
        
        # If we don't have enough diversity, sample more strategically
        if (len(unique_cells) < min_cells or 
            len(unique_genes) < min_genes or 
            len(unique_drugs) < min_drugs):
            logger.warning("Diversity requirements not fully met, but proceeding with available data")
        
        return df
    
    def prepare_gene_encoder(self, data: pd.DataFrame) -> LabelEncoder:
        """Create gene encoder for embedding-based models"""
        unique_genes = data['gene_name'].unique()
        encoder = LabelEncoder()
        encoder.fit(unique_genes)
        logger.info(f"Created gene encoder for {len(unique_genes)} unique genes")
        return encoder
    
    def prepare_data(self, data: pd.DataFrame, drug_features: Dict[str, np.ndarray], 
                    test_size: float = 0.2, random_state: int = 42) -> Tuple[DataLoader, DataLoader]:
        """Prepare optimized data loaders"""
        logger.info("Preparing optimized data for training")
        
        # Filter data
        valid_drugs = set(drug_features.keys())
        filtered_data = data[data['drug'].isin(valid_drugs)].copy()
        
        logger.info(f"Filtered to {len(filtered_data)} samples with known drug features")
        
        # Create gene encoder
        self.gene_encoder = self.prepare_gene_encoder(filtered_data)
        
        # Split data
        train_data, test_data = train_test_split(
            filtered_data, test_size=test_size, random_state=random_state, 
            stratify=filtered_data['drug']
        )
        
        # Create optimized datasets
        use_embeddings = self.model_type in ['gene_embedding', 'attention_fusion']
        
        train_dataset = OptimizedTahoeDrugDataset(
            train_data, drug_features, self.gene_encoder, 
            self.gene_features, use_embeddings
        )
        test_dataset = OptimizedTahoeDrugDataset(
            test_data, drug_features, self.gene_encoder, 
            self.gene_features, use_embeddings
        )
        
        # Optimized data loaders with larger batch sizes
        batch_size = 1024 if torch.cuda.is_available() else 512
        num_workers = min(8, torch.get_num_threads())
        
        train_loader = DataLoader(
            train_dataset, batch_size=batch_size, shuffle=True, 
            num_workers=num_workers, pin_memory=True, persistent_workers=True
        )
        test_loader = DataLoader(
            test_dataset, batch_size=batch_size, shuffle=False, 
            num_workers=num_workers, pin_memory=True, persistent_workers=True
        )
        
        logger.info(f"Created data loaders with batch size {batch_size} and {num_workers} workers")
        
        return train_loader, test_loader
    
    def create_model(self, drug_dim: int, num_genes: int = None) -> nn.Module:
        """Create optimized model based on type"""
        if self.model_type == 'gene_embedding':
            return GeneEmbeddingMLP(drug_dim, num_genes).to(self.device)
        elif self.model_type == 'residual':
            input_dim = drug_dim + 1000  # drug + gene features
            return ResidualMLP(input_dim).to(self.device)
        elif self.model_type == 'attention_fusion':
            gene_dim = 1000  # Fixed gene feature dimension
            return AttentionFusionModel(drug_dim, gene_dim).to(self.device)
        else:
            raise ValueError(f"Unknown model type: {self.model_type}")
    
    def train_model(self, train_loader: DataLoader, test_loader: DataLoader, 
                   epochs: int = 50, learning_rate: float = 0.001) -> Dict:
        """Optimized training with better techniques"""
        logger.info(f"Training {self.model_type} model for {epochs} epochs")
        
        # Create model
        drug_dim = len(list(self.drug_features.values())[0])
        num_genes = len(self.gene_encoder.classes_) if self.gene_encoder else None
        
        self.model = self.create_model(drug_dim, num_genes)
        
        # Count parameters
        total_params = sum(p.numel() for p in self.model.parameters())
        trainable_params = sum(p.numel() for p in self.model.parameters() if p.requires_grad)
        logger.info(f"Model has {total_params:,} total parameters ({trainable_params:,} trainable)")
        
        # Setup training with better optimizer
        criterion = nn.MSELoss()
        optimizer = optim.AdamW(self.model.parameters(), lr=learning_rate, weight_decay=1e-4)
        scheduler = optim.lr_scheduler.CosineAnnealingWarmRestarts(optimizer, T_0=10, T_mult=2)
        
        # Training history
        history = {
            'train_loss': [],
            'test_loss': [],
            'train_r2': [],
            'test_r2': [],
            'learning_rate': []
        }
        
        best_test_loss = float('inf')
        patience = 10
        patience_counter = 0
        
        for epoch in range(epochs):
            # Training phase
            self.model.train()
            train_loss = 0.0
            train_predictions = []
            train_targets = []
            
            for batch_idx, batch in enumerate(train_loader):
                if self.model_type in ['gene_embedding', 'attention_fusion']:
                    drug_features, gene_indices, targets = batch
                    drug_features = drug_features.to(self.device, non_blocking=True)
                    gene_indices = gene_indices.to(self.device, non_blocking=True)
                    targets = targets.to(self.device, non_blocking=True)
                    
                    optimizer.zero_grad()
                    outputs = self.model(drug_features, gene_indices)
                    loss = criterion(outputs, targets.squeeze())
                    loss.backward()
                    
                    # Gradient clipping
                    torch.nn.utils.clip_grad_norm_(self.model.parameters(), max_norm=1.0)
                    
                    optimizer.step()
                    
                    train_loss += loss.item()
                    train_predictions.extend(outputs.detach().cpu().numpy())
                    train_targets.extend(targets.detach().cpu().numpy())
                
                else:  # residual model
                    features, targets = batch
                    features = features.to(self.device, non_blocking=True)
                    targets = targets.to(self.device, non_blocking=True)
                    
                    optimizer.zero_grad()
                    outputs = self.model(features)
                    loss = criterion(outputs, targets.squeeze())
                    loss.backward()
                    
                    torch.nn.utils.clip_grad_norm_(self.model.parameters(), max_norm=1.0)
                    
                    optimizer.step()
                    
                    train_loss += loss.item()
                    train_predictions.extend(outputs.detach().cpu().numpy())
                    train_targets.extend(targets.detach().cpu().numpy())
            
            # Validation phase
            self.model.eval()
            test_loss = 0.0
            test_predictions = []
            test_targets = []
            
            with torch.no_grad():
                for batch in test_loader:
                    if self.model_type in ['gene_embedding', 'attention_fusion']:
                        drug_features, gene_indices, targets = batch
                        drug_features = drug_features.to(self.device, non_blocking=True)
                        gene_indices = gene_indices.to(self.device, non_blocking=True)
                        targets = targets.to(self.device, non_blocking=True)
                        
                        outputs = self.model(drug_features, gene_indices)
                        loss = criterion(outputs, targets.squeeze())
                        
                        test_loss += loss.item()
                        test_predictions.extend(outputs.cpu().numpy())
                        test_targets.extend(targets.cpu().numpy())
                    
                    else:  # residual model
                        features, targets = batch
                        features = features.to(self.device, non_blocking=True)
                        targets = targets.to(self.device, non_blocking=True)
                        
                        outputs = self.model(features)
                        loss = criterion(outputs, targets.squeeze())
                        
                        test_loss += loss.item()
                        test_predictions.extend(outputs.cpu().numpy())
                        test_targets.extend(targets.cpu().numpy())
            
            # Calculate metrics
            train_loss /= len(train_loader)
            test_loss /= len(test_loader)
            
            train_r2 = r2_score(train_targets, train_predictions)
            test_r2 = r2_score(test_targets, test_predictions)
            
            # Update history
            history['train_loss'].append(train_loss)
            history['test_loss'].append(test_loss)
            history['train_r2'].append(train_r2)
            history['test_r2'].append(test_r2)
            history['learning_rate'].append(optimizer.param_groups[0]['lr'])
            
            # Learning rate scheduling
            scheduler.step()
            
            # Early stopping
            if test_loss < best_test_loss:
                best_test_loss = test_loss
                patience_counter = 0
                self.save_model('best_optimized_model.pth')
            else:
                patience_counter += 1
            
            if patience_counter >= patience:
                logger.info(f"Early stopping at epoch {epoch}")
                break
            
            if epoch % 5 == 0:
                logger.info(f"Epoch {epoch:3d}: Train Loss: {train_loss:.4f}, Test Loss: {test_loss:.4f}, "
                          f"Train R²: {train_r2:.4f}, Test R²: {test_r2:.4f}, LR: {optimizer.param_groups[0]['lr']:.6f}")
        
        return history
    
    def save_model(self, path: str):
        """Save model with all necessary components"""
        torch.save({
            'model_state_dict': self.model.state_dict(),
            'model_type': self.model_type,
            'drug_features': self.drug_features,
            'gene_encoder': self.gene_encoder,
            'scaler': self.scaler
        }, path)
        logger.info(f"Model saved to {path}")
    
    def load_model(self, path: str):
        """Load model with all components"""
        checkpoint = torch.load(path, map_location=self.device)
        self.model_type = checkpoint['model_type']
        self.drug_features = checkpoint['drug_features']
        self.gene_encoder = checkpoint['gene_encoder']
        self.scaler = checkpoint['scaler']
        
        # Recreate model
        drug_dim = len(list(self.drug_features.values())[0])
        num_genes = len(self.gene_encoder.classes_) if self.gene_encoder else None
        
        self.model = self.create_model(drug_dim, num_genes)
        self.model.load_state_dict(checkpoint['model_state_dict'])
        self.model.to(self.device)
        logger.info(f"Model loaded from {path}")
    
    def predict(self, drug_name: str, gene_name: str) -> float:
        """Make prediction for a specific drug-gene pair"""
        if drug_name not in self.drug_features:
            raise ValueError(f"Drug {drug_name} not found in features")
        
        self.model.eval()
        
        # Prepare features
        drug_feat = torch.FloatTensor(self.drug_features[drug_name]).unsqueeze(0).to(self.device)
        
        with torch.no_grad():
            if self.model_type in ['gene_embedding', 'attention_fusion']:
                try:
                    gene_idx = self.gene_encoder.transform([gene_name])[0]
                    gene_tensor = torch.LongTensor([gene_idx]).to(self.device)
                    prediction = self.model(drug_feat, gene_tensor).item()
                except ValueError:
                    # Handle unseen genes
                    gene_tensor = torch.LongTensor([0]).to(self.device)
                    prediction = self.model(drug_feat, gene_tensor).item()
            else:
                # For residual model, use one-hot encoding
                gene_feat = torch.zeros(1000)
                gene_hash = hash(gene_name) % 1000
                gene_feat[gene_hash] = 1.0
                features = torch.cat([drug_feat.squeeze(0), gene_feat]).unsqueeze(0).to(self.device)
                prediction = self.model(features).item()
        
        return prediction
    
    def plot_training_history(self, history: Dict, save_path: str = 'optimized_training_history.png'):
        """Enhanced plotting with learning rate"""
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 10))
        
        # Loss plot
        ax1.plot(history['train_loss'], label='Train Loss')
        ax1.plot(history['test_loss'], label='Test Loss')
        ax1.set_xlabel('Epoch')
        ax1.set_ylabel('MSE Loss')
        ax1.set_title('Training and Validation Loss')
        ax1.legend()
        ax1.grid(True)
        
        # R² plot
        ax2.plot(history['train_r2'], label='Train R²')
        ax2.plot(history['test_r2'], label='Test R²')
        ax2.set_xlabel('Epoch')
        ax2.set_ylabel('R² Score')
        ax2.set_title('Training and Validation R²')
        ax2.legend()
        ax2.grid(True)
        
        # Learning rate plot
        ax3.plot(history['learning_rate'], label='Learning Rate')
        ax3.set_xlabel('Epoch')
        ax3.set_ylabel('Learning Rate')
        ax3.set_title('Learning Rate Schedule')
        ax3.legend()
        ax3.grid(True)
        
        # Loss difference plot
        loss_diff = [t - v for t, v in zip(history['train_loss'], history['test_loss'])]
        ax4.plot(loss_diff, label='Train - Test Loss')
        ax4.set_xlabel('Epoch')
        ax4.set_ylabel('Loss Difference')
        ax4.set_title('Overfitting Indicator')
        ax4.legend()
        ax4.grid(True)
        
        plt.tight_layout()
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        plt.show()
        logger.info(f"Training history plot saved to {save_path}")

def main():
    """
    Main optimized training pipeline
    """
    # Enhanced configuration
    config = {
        'model_type': 'gene_embedding',  # 'gene_embedding', 'residual', 'attention_fusion'
        'max_samples': 500000,   # Reduced for faster sampling
        'epochs': 50,
        'learning_rate': 0.001,
        'test_size': 0.2,
        'drug_descriptors_path': 'drug_descriptors_scaled.csv',
        'use_stratified': True,  # Use fast diversity sampling
        'min_cells': 1,          # Minimum number of cell lines (relaxed)
        'min_genes': 1000,       # Minimum number of genes
        'min_drugs': 20          # Minimum number of drugs (relaxed)
    }
    
    logger.info("Starting Optimized Tahoe-100M Drug-Gene Expression Prediction Training")
    logger.info(f"Configuration: {config}")
    
    # Initialize optimized trainer
    trainer = OptimizedTahoeMLPTrainer(model_type=config['model_type'])
    
    # Load drug descriptors
    drug_features, feature_columns = trainer.load_drug_descriptors(config['drug_descriptors_path'])
    trainer.drug_features = drug_features
    
    # Stream Tahoe data with stratified sampling
    tahoe_data = trainer.stream_tahoe_data(
        max_samples=config['max_samples'],
        use_stratified=config['use_stratified'],
        min_cells=config['min_cells'],
        min_genes=config['min_genes'],
        min_drugs=config['min_drugs']
    )
    
    # Prepare data
    train_loader, test_loader = trainer.prepare_data(
        tahoe_data, drug_features, test_size=config['test_size']
    )
    
    # Train model
    history = trainer.train_model(
        train_loader, test_loader, 
        epochs=config['epochs'], 
        learning_rate=config['learning_rate']
    )
    
    # Plot results
    trainer.plot_training_history(history)
    
    # Example predictions
    try:
        sample_drugs = list(drug_features.keys())[:3]
        sample_genes = tahoe_data['gene_name'].unique()[:3]
        
        logger.info("Example predictions:")
        for drug in sample_drugs:
            for gene in sample_genes:
                try:
                    prediction = trainer.predict(drug, gene)
                    logger.info(f"{drug} + {gene}: {prediction:.4f}")
                except Exception as e:
                    logger.warning(f"Could not predict {drug} + {gene}: {e}")
    except Exception as e:
        logger.warning(f"Could not make example predictions: {e}")
    
    logger.info("Optimized training completed successfully!")

if __name__ == "__main__":
    main()
