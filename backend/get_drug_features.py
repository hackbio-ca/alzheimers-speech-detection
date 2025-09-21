#!/usr/bin/env python3
"""
Simple script that outputs just the feature vector for a drug.
Usage: python get_drug_features.py <drug_name>
Output: Just the feature values (no drug name, SMILES, etc.)
"""

import sys
from standalone_drug_descriptors import get_drug_descriptors

def main():
    if len(sys.argv) < 2:
        print("Usage: python get_drug_features.py <drug_name>")
        print("Example: python get_drug_features.py aspirin")
        sys.exit(1)
    
    drug_name = sys.argv[1]
    
    # Get descriptors with no verbose output and no scaling to see raw values
    result = get_drug_descriptors(drug_name, scale_descriptors_flag=False, verbose=False)
    
    if result is None:
        sys.exit(1)
    
    # Extract just the feature vector (skip drug_name, smiles, cid, synonyms)
    feature_vector = result[4:]
    
    # Output just the feature values
    print("Raw (unscaled) feature values:")
    print(feature_vector)
    
    # Also show scaled version for comparison
    result_scaled = get_drug_descriptors(drug_name, scale_descriptors_flag=True, verbose=False)
    if result_scaled:
        feature_vector_scaled = result_scaled[4:]
        print("\nScaled feature values:")
        print(feature_vector_scaled)

if __name__ == "__main__":
    main()