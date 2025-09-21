#!/usr/bin/env python3
"""
Extract all unique drug names and their canonical SMILES from the Tahoe-100M dataset using streaming to conserve memory.
This script processes the drug_metadata subset in chunks to conserve memory while collecting
all unique drug names with their corresponding SMILES structures.
"""

from datasets import load_dataset
from tqdm import tqdm
import json
import sys

def extract_unique_drugs(batch_size=1000, max_samples=None):
    """
    Extract all unique drug names and their canonical SMILES from the Tahoe-100M drug_metadata subset.
    
    Args:
        batch_size (int): Number of samples to process at once
        max_samples (int): Maximum number of samples to process (None for all)
    
    Returns:
        set: Set of unique drug entries in format "drug_name,canonical_smiles"
    """
    print("Loading Tahoe-100M drug_metadata subset...")
    ds = load_dataset("tahoebio/Tahoe-100M", "drug_metadata", split="train", streaming=True)
    
    unique_drugs = set()
    processed_count = 0
    
    # Create an iterator over the dataset
    dataset_iter = iter(ds)
    
    # Use tqdm for progress tracking
    if max_samples:
        pbar = tqdm(total=max_samples, desc="Processing samples")
    else:
        pbar = tqdm(desc="Processing samples")
    
    try:
        while True:
            # Check if we've reached the maximum number of samples
            if max_samples and processed_count >= max_samples:
                break
                
            # Process samples in batches
            batch = []
            for _ in range(batch_size):
                try:
                    sample = next(dataset_iter)
                    batch.append(sample)
                    processed_count += 1
                    
                    if max_samples and processed_count >= max_samples:
                        break
                        
                except StopIteration:
                    # End of dataset
                    break
            
            if not batch:
                break
            
            # Extract drug names and canonical_smiles from the batch
            for sample in batch:
                if 'drug' in sample and sample['drug']:
                    drug_name = sample['drug']
                    canonical_smiles = sample.get('canonical_smiles', '')
                    # Store as "drug_name,canonical_smiles"
                    drug_entry = f"{drug_name},{canonical_smiles}"
                    unique_drugs.add(drug_entry)
            
            # Update progress bar
            pbar.update(len(batch))
            
            # Print progress every 10,000 samples
            if processed_count % 10000 == 0:
                print(f"\nProcessed {processed_count:,} samples, found {len(unique_drugs):,} unique drugs")
    
    except KeyboardInterrupt:
        print(f"\nInterrupted by user. Processed {processed_count:,} samples.")
    
    finally:
        pbar.close()
    
    return unique_drugs, processed_count

def main():
    """Main function to extract unique drugs and their SMILES structures and save results."""
    print("Starting extraction of unique drug names and SMILES from Tahoe-100M drug_metadata subset...")
    
    # You can limit the number of samples for testing by setting max_samples
    # For example: max_samples=100000 for first 100k samples
    max_samples = None  # Set to None to process all samples
    
    try:
        unique_drugs, total_processed = extract_unique_drugs(max_samples=max_samples)
        
        print(f"\nExtraction complete!")
        print(f"Total samples processed: {total_processed:,}")
        print(f"Unique drugs found: {len(unique_drugs):,}")
        
        # Convert set to sorted list for better readability
        unique_drugs_list = sorted(list(unique_drugs))
        
        # Save to JSON file
        output_file = "unique_drugs_tahoe100m.json"
        with open(output_file, 'w', encoding='utf-8') as f:
            json.dump({
                'total_samples_processed': total_processed,
                'unique_drugs_count': len(unique_drugs),
                'unique_drugs': unique_drugs_list
            }, f, indent=2, ensure_ascii=False)
        
        print(f"Results saved to: {output_file}")
        
        # Also save as a simple text file with one drug per line
        txt_output_file = "unique_drugs_tahoe100m.txt"
        with open(txt_output_file, 'w', encoding='utf-8') as f:
            for drug_entry in unique_drugs_list:
                f.write(f"{drug_entry}\n")
        
        print(f"Drug list also saved to: {txt_output_file}")
        
        # Print first 20 drugs as a preview
        print(f"\nFirst 20 unique drugs found:")
        for i, drug_entry in enumerate(unique_drugs_list[:20], 1):
            drug_name, canonical_smiles = drug_entry.split(',', 1)
            print(f"{i:2d}. {drug_name}")
            if canonical_smiles:
                print(f"    SMILES: {canonical_smiles}")
            else:
                print(f"    SMILES: (not available)")
        
        if len(unique_drugs_list) > 20:
            print(f"... and {len(unique_drugs_list) - 20:,} more")
            
    except Exception as e:
        print(f"Error occurred: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
