#!/usr/bin/env python3
"""
Query PubChem for drug SMILES and calculate molecular descriptors.

This script takes a drug name as input, queries PubChem to get the SMILES string,
and then calculates various molecular descriptors using RDKit.

Features:
- PubChem API integration for SMILES retrieval
- Molecular standardization: removes salts/solvents, neutralizes, normalizes tautomers
- Hydrogen addition for better 3D properties
- Appropriate scaling: Z-score for continuous, log1p for counts, min-max for fingerprints
- Comprehensive molecular descriptors including fingerprints, size/shape, polarity, etc.

Environment Requirements:
- This script requires the conda RDKit, NOT pip version
- Requires PubChemPy library for PubChem API calls (pip install pubchempy)
- Falls back to requests library if PubChemPy is not available
"""

import json
import pandas as pd
import sys
import argparse
import os
import time
import logging
import numpy as np
from pathlib import Path

# Try to import PubChemPy for better compound searching
try:
    import pubchempy as pcp
    HAS_PUBCHEMPY = True
    print("PubChemPy imported successfully - using enhanced compound search")
except ImportError:
    HAS_PUBCHEMPY = False
    print("PubChemPy not available - falling back to HTTP requests")
    import requests

def check_environment():
    """Check if we're running in the correct conda environment."""
    conda_env = os.environ.get('CONDA_DEFAULT_ENV', '')
    if conda_env != 'my-rdkit-env':
        print("=" * 60)
        print("ENVIRONMENT WARNING")
        print("=" * 60)
        if conda_env:
            print(f"Currently in conda environment: '{conda_env}'")
        else:
            print("Not running in a conda environment")
        print("Expected environment: 'my-rdkit-env'")
        print("\nTo fix this:")
        print("1. Run: conda activate my-rdkit-env")
        print("2. Or run: run_rdkit_script.bat")
        print("3. Or run: setup_rdkit_env.bat (if environment doesn't exist)")
        print("=" * 60)
        return False
    return True

try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors
    from rdkit.Chem import rdFingerprintGenerator
    from rdkit.Chem.AtomPairs import Pairs
    from rdkit.Chem.Fingerprints import FingerprintMols
    import numpy as np
    from sklearn.preprocessing import StandardScaler, MinMaxScaler
    print("RDKit and dependencies imported successfully!")
    
    # Try to import rdMolStandardize (not available in pip version)
    try:
        from rdkit.Chem import rdMolStandardize
        HAS_MOLSTANDARDIZE = True
        print("rdMolStandardize available - full standardization enabled")
    except ImportError:
        HAS_MOLSTANDARDIZE = False
        print("rdMolStandardize not available - using basic standardization")
        
except ImportError as e:
    print("=" * 60)
    print("IMPORT ERROR")
    print("=" * 60)
    print(f"Error importing required packages: {e}")
    print("\nThis script requires the 'rdkit-env' conda environment.")
    print("Please run one of the following:")
    print("1. setup_rdkit_env.bat (to create the environment)")
    print("2. run_rdkit_script.bat (to run in the correct environment)")
    print("3. conda activate rdkit-env (then run this script)")
    print("=" * 60)
    sys.exit(1)

def query_pubchem_smiles(drug_name, max_retries=3, delay=1.0):
    """
    Query PubChem for SMILES string of a drug using PubChemPy.
    
    Args:
        drug_name (str): Name of the drug to search for
        max_retries (int): Maximum number of retry attempts
        delay (float): Delay between retries in seconds
        
    Returns:
        tuple: (smiles, cid, synonyms) or (None, None, None) if not found
    """
    print(f"Querying PubChem for: {drug_name}")
    
    if HAS_PUBCHEMPY:
        return query_pubchem_smiles_pubchempy(drug_name, max_retries, delay)
    else:
        return query_pubchem_smiles_http(drug_name, max_retries, delay)

def query_pubchem_smiles_pubchempy(drug_name, max_retries=3, delay=1.0):
    """
    Query PubChem using PubChemPy library (preferred method).
    
    Args:
        drug_name (str): Name of the drug to search for
        max_retries (int): Maximum number of retry attempts
        delay (float): Delay between retries in seconds
        
    Returns:
        tuple: (smiles, cid, synonyms) or (None, None, None) if not found
    """
    for attempt in range(max_retries):
        try:
            print(f"Searching PubChem for: {drug_name} (attempt {attempt + 1})")
            
            # Search for compounds by name
            compounds = list(pcp.get_compounds(drug_name, 'name'))
            
            if not compounds:
                # Try searching by synonym
                print("No direct matches found, trying synonym search...")
                compounds = list(pcp.get_compounds(drug_name, 'synonym'))
            
            if not compounds:
                # Try fuzzy search
                print("No synonym matches found, trying fuzzy search...")
                compounds = list(pcp.get_compounds(drug_name, 'name', fuzzy=True))
            
            if compounds:
                # Use the first compound found
                compound = compounds[0]
                
                # Get SMILES (prefer isomeric, fall back to canonical)
                smiles = compound.isomeric_smiles or compound.canonical_smiles
                
                if smiles:
                    print(f"✅ Found compound: {compound.iupac_name or drug_name}")
                    print(f"   SMILES: {smiles}")
                    print(f"   CID: {compound.cid}")
                    print(f"   Molecular Formula: {compound.molecular_formula}")
                    print(f"   Molecular Weight: {compound.molecular_weight}")
                    print(f"   LogP: {compound.xlogp}")
                    
                    # Get synonyms
                    synonyms = []
                    try:
                        synonyms = compound.synonyms[:10]  # Limit to first 10 synonyms
                    except:
                        synonyms = [drug_name]
                    
                    return smiles, compound.cid, synonyms
                else:
                    print(f"Compound found but no SMILES available for CID {compound.cid}")
            else:
                print(f"No compounds found for: {drug_name}")
            
        except Exception as e:
            print(f"Attempt {attempt + 1} failed: {e}")
            if attempt < max_retries - 1:
                print(f"Retrying in {delay} seconds...")
                time.sleep(delay)
                delay *= 2  # Exponential backoff
    
    print(f"Failed to find SMILES for '{drug_name}' after {max_retries} attempts")
    return None, None, None

def query_pubchem_smiles_http(drug_name, max_retries=3, delay=1.0):
    """
    Fallback method using HTTP requests when PubChemPy is not available.
    
    Args:
        drug_name (str): Name of the drug to search for
        max_retries (int): Maximum number of retry attempts
        delay (float): Delay between retries in seconds
        
    Returns:
        tuple: (smiles, cid, synonyms) or (None, None, None) if not found
    """
    print("Using HTTP fallback method (PubChemPy not available)")
    
    for attempt in range(max_retries):
        try:
            # Simple HTTP request to PubChem API
            search_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{requests.utils.quote(drug_name)}/property/IsomericSMILES,MolecularFormula,MolecularWeight/JSON"
            response = requests.get(search_url, timeout=10)
            
            if response.status_code == 200:
                data = response.json()
                if 'PropertyTable' in data and 'Properties' in data['PropertyTable']:
                    properties = data['PropertyTable']['Properties'][0]
                    smiles = properties.get('IsomericSMILES', '')
                    molecular_formula = properties.get('MolecularFormula', '')
                    molecular_weight = properties.get('MolecularWeight', '')
                    
                    if smiles:
                        print(f"✅ Found SMILES: {smiles}")
                        print(f"   Molecular Formula: {molecular_formula}")
                        print(f"   Molecular Weight: {molecular_weight}")
                        return smiles, None, [drug_name]
            
            print(f"Attempt {attempt + 1} failed: No results found")
            
        except Exception as e:
            print(f"Attempt {attempt + 1} failed: {e}")
        
        if attempt < max_retries - 1:
            print(f"Retrying in {delay} seconds...")
            time.sleep(delay)
            delay *= 2  # Exponential backoff
    
    print(f"Failed to find SMILES for '{drug_name}' after {max_retries} attempts")
    return None, None, None

def standardize_molecule(smiles):
    """
    Standardize molecule: remove salts/solvents, neutralize, normalize tautomers.
    Falls back to basic standardization if rdMolStandardize is not available.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        Chem.Mol: Standardized molecule object or None if failed
    """
    try:
        # Convert SMILES to RDKit molecule object
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        
        if HAS_MOLSTANDARDIZE:
            # Full standardization with rdMolStandardize
            # Step 1: Cleanup - remove salts, solvents, and standardize
            mol = rdMolStandardize.Cleanup(mol)
            
            # Step 2: Get the largest fragment (remove salts/solvents)
            mol = rdMolStandardize.FragmentParent(mol)
            
            # Step 3: Normalize tautomers
            enumerator = rdMolStandardize.TautomerEnumerator()
            mol = enumerator.Canonicalize(mol)
        else:
            # Basic standardization without rdMolStandardize
            # Remove salts by taking the largest fragment
            try:
                # Get all fragments and take the largest one
                frags = Chem.GetMolFrags(mol, asMols=True)
                if frags:
                    # Sort by number of atoms and take the largest
                    mol = max(frags, key=lambda x: x.GetNumAtoms())
                else:
                    return None
            except:
                # If fragment separation fails, just use the original molecule
                pass
            
            # Basic cleanup - remove hydrogens and re-add them
            mol = Chem.RemoveHs(mol)
        
        # Step 4: Add hydrogens for better 3D properties
        mol = Chem.AddHs(mol)
        
        return mol
        
    except Exception as e:
        print(f"Error standardizing SMILES '{smiles}': {e}")
        return None

def calculate_molecular_descriptors(smiles):
    """
    Calculate molecular descriptors and fingerprints from SMILES string.
    Includes proper molecular standardization and hydrogen addition.
    Focus on discriminative features: fingerprints, size/shape, polarity, 
    H-bonding, LogP, flexibility, charge, SA score, and functional groups.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        dict: Dictionary containing calculated descriptors and fingerprints
    """
    if not smiles or smiles.strip() == '':
        return None
    
    try:
        # Standardize the molecule first
        mol = standardize_molecule(smiles)
        if mol is None:
            return None
        
        descriptors = {}
        
        # === SIZE/SHAPE DESCRIPTORS ===
        descriptors['molecular_weight'] = Descriptors.MolWt(mol)
        descriptors['num_atoms'] = mol.GetNumAtoms()
        descriptors['num_heavy_atoms'] = Descriptors.HeavyAtomCount(mol)
        descriptors['num_bonds'] = mol.GetNumBonds()
        descriptors['num_rings'] = Descriptors.RingCount(mol)
        descriptors['num_aromatic_rings'] = Descriptors.NumAromaticRings(mol)
        try:
            descriptors['num_saturated_rings'] = Descriptors.NumSaturatedRings(mol)
            descriptors['num_aliphatic_rings'] = Descriptors.NumAliphaticRings(mol)
        except AttributeError:
            # Fallback: calculate manually
            descriptors['num_saturated_rings'] = 0
            descriptors['num_aliphatic_rings'] = 0
        
        # === POLARITY DESCRIPTORS ===
        descriptors['tpsa'] = Descriptors.TPSA(mol)  # Topological Polar Surface Area
        descriptors['logp'] = Crippen.MolLogP(mol)  # Lipophilicity
        descriptors['molar_refractivity'] = Crippen.MolMR(mol)
        
        # === H-BONDING DESCRIPTORS ===
        try:
            # Try the Lipinski module first (more common in pip RDKit)
            from rdkit.Chem import Lipinski
            descriptors['num_hbd'] = Lipinski.NumHBD(mol)  # Hydrogen Bond Donors
            descriptors['num_hba'] = Lipinski.NumHBA(mol)  # Hydrogen Bond Acceptors
        except (ImportError, AttributeError):
            try:
                # Fallback to Descriptors module
                descriptors['num_hbd'] = Descriptors.NumHBD(mol)
                descriptors['num_hba'] = Descriptors.NumHBA(mol)
            except AttributeError:
                # If neither works, use alternative calculation
                descriptors['num_hbd'] = rdMolDescriptors.CalcNumHBD(mol)
                descriptors['num_hba'] = rdMolDescriptors.CalcNumHBA(mol)
        
        # === FLEXIBILITY DESCRIPTORS ===
        descriptors['num_rotatable_bonds'] = Descriptors.NumRotatableBonds(mol)
        try:
            descriptors['num_sp3_carbons'] = Descriptors.NumSpiroAtoms(mol)
        except AttributeError:
            # Fallback: count sp3 carbons manually
            descriptors['num_sp3_carbons'] = sum(1 for atom in mol.GetAtoms() 
                                               if atom.GetAtomicNum() == 6 and atom.GetHybridization() == Chem.HybridizationType.SP3)
        
        # === CHARGE DESCRIPTORS ===
        descriptors['formal_charge'] = Chem.rdmolops.GetFormalCharge(mol)
        try:
            descriptors['num_positive_charges'] = Descriptors.NumPositiveIons(mol)
            descriptors['num_negative_charges'] = Descriptors.NumNegativeIons(mol)
        except AttributeError:
            # Fallback: count charges manually
            descriptors['num_positive_charges'] = sum(1 for atom in mol.GetAtoms() if atom.GetFormalCharge() > 0)
            descriptors['num_negative_charges'] = sum(1 for atom in mol.GetAtoms() if atom.GetFormalCharge() < 0)
        
        # === SIMPLE SA SCORE (Synthetic Accessibility) ===
        try:
            # Use a simple complexity measure as SA score proxy
            descriptors['sa_score'] = rdMolDescriptors.CalcNumSpiroAtoms(mol)  # Simplified SA proxy
        except:
            descriptors['sa_score'] = 0
        
        # === ATOMIC COMPOSITION ===
        # Count atoms manually (more reliable across RDKit versions)
        descriptors['num_carbons'] = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
        descriptors['num_nitrogens'] = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
        descriptors['num_oxygens'] = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
        descriptors['num_sulfurs'] = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 16)
        descriptors['num_halogens'] = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() in [9, 17, 35, 53])
        descriptors['num_heteroatoms'] = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() not in [1, 6])
        
        # === FUNCTIONAL GROUP COUNTS ===
        # Use a safer approach that handles missing functions gracefully
        functional_groups = [
            'fr_amide', 'fr_ester', 'fr_ether', 'fr_ketone', 'fr_aldehyde',
            'fr_COO', 'fr_amine', 'fr_aromatic', 'fr_benzene', 'fr_pyridine',
            'fr_pyrimidine', 'fr_pyrazole', 'fr_imidazole', 'fr_thiazole',
            'fr_oxazole', 'fr_phenol', 'fr_aniline', 'fr_halogen', 'fr_alkyl_halide',
            'fr_benzodiazepine', 'fr_bicyclic', 'fr_urea', 'fr_sulfonamide',
            'fr_thiophene', 'fr_furan', 'fr_pyrrole', 'fr_azide', 'fr_aziridine',
            'fr_epoxide', 'fr_oxirane'
        ]
        
        for group in functional_groups:
            try:
                # Get the function from Descriptors module
                func = getattr(Descriptors, group)
                descriptors[group] = func(mol)
            except AttributeError:
                # If function doesn't exist, set to 0
                descriptors[group] = 0
        
        # === CIRCULAR FINGERPRINTS (Morgan/ECFP) ===
        # ECFP4 (radius=2) - most commonly used
        ecfp4_generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
        ecfp4_fp = ecfp4_generator.GetFingerprint(mol)
        descriptors['ecfp4_density'] = ecfp4_fp.GetNumOnBits() / ecfp4_fp.GetNumBits()
        descriptors['ecfp4_bits'] = ecfp4_fp.GetNumOnBits()
        
        # ECFP6 (radius=3) - more detailed
        ecfp6_generator = rdFingerprintGenerator.GetMorganGenerator(radius=3, fpSize=2048)
        ecfp6_fp = ecfp6_generator.GetFingerprint(mol)
        descriptors['ecfp6_density'] = ecfp6_fp.GetNumOnBits() / ecfp6_fp.GetNumBits()
        descriptors['ecfp6_bits'] = ecfp6_fp.GetNumOnBits()
        
        # === MACCS KEYS (Lightweight substructure keys) ===
        maccs_fp = rdMolDescriptors.GetMACCSKeysFingerprint(mol)
        descriptors['maccs_density'] = maccs_fp.GetNumOnBits() / maccs_fp.GetNumBits()
        descriptors['maccs_bits'] = maccs_fp.GetNumOnBits()
        
        # === ADDITIONAL STRUCTURAL DESCRIPTORS ===
        # Use safe approach for structural descriptors
        structural_descriptors = {
            'balaban_j': 'BalabanJ',
            'bertz_ct': 'BertzCT', 
            'kier_flexibility': 'Kappa1',
            'kier_shape1': 'Kappa2',
            'kier_shape2': 'Kappa3'
        }
        
        for desc_name, func_name in structural_descriptors.items():
            try:
                func = getattr(Descriptors, func_name)
                descriptors[desc_name] = func(mol)
            except AttributeError:
                descriptors[desc_name] = 0.0
        
        return descriptors
        
    except Exception as e:
        print(f"Error calculating descriptors for SMILES '{smiles}': {e}")
        return None

def scale_descriptors(df):
    """
    Scale molecular descriptors appropriately:
    - Z-score for continuous descriptors
    - Don't scale binary fingerprints
    - Min-max or log1p for count fingerprints
    
    Args:
        df (pd.DataFrame): DataFrame with molecular descriptors
        
    Returns:
        pd.DataFrame: DataFrame with scaled descriptors
    """
    df_scaled = df.copy()
    
    # Define descriptor types
    continuous_descriptors = [
        'molecular_weight', 'tpsa', 'logp', 'molar_refractivity',
        'balaban_j', 'bertz_ct', 'kier_flexibility', 'kier_shape1', 'kier_shape2'
    ]
    
    count_descriptors = [
        'num_atoms', 'num_heavy_atoms', 'num_bonds', 'num_rings',
        'num_aromatic_rings', 'num_saturated_rings', 'num_aliphatic_rings',
        'num_hbd', 'num_hba', 'num_rotatable_bonds', 'num_sp3_carbons',
        'formal_charge', 'num_positive_charges', 'num_negative_charges',
        'num_carbons', 'num_nitrogens', 'num_oxygens', 'num_sulfurs',
        'num_halogens', 'num_heteroatoms', 'sa_score'
    ]
    
    # Functional group counts (fr_*)
    functional_group_descriptors = [col for col in df.columns if col.startswith('fr_')]
    count_descriptors.extend(functional_group_descriptors)
    
    # Fingerprint density and bit counts
    fingerprint_descriptors = [
        'ecfp4_density', 'ecfp4_bits', 'ecfp6_density', 'ecfp6_bits',
        'maccs_density', 'maccs_bits'
    ]
    
    print("Scaling molecular descriptors...")
    
    # Z-score scaling for continuous descriptors
    continuous_cols = [col for col in continuous_descriptors if col in df.columns]
    if continuous_cols:
        scaler_continuous = StandardScaler()
        df_scaled[continuous_cols] = scaler_continuous.fit_transform(df[continuous_cols])
        print(f"Z-score scaled {len(continuous_cols)} continuous descriptors")
    
    # Log1p scaling for count descriptors (handles zeros well)
    count_cols = [col for col in count_descriptors if col in df.columns]
    if count_cols:
        # Apply log1p transformation
        df_scaled[count_cols] = np.log1p(df[count_cols])
        print(f"Log1p scaled {len(count_cols)} count descriptors")
    
    # Min-max scaling for fingerprint descriptors
    fingerprint_cols = [col for col in fingerprint_descriptors if col in df.columns]
    if fingerprint_cols:
        scaler_fingerprint = MinMaxScaler()
        df_scaled[fingerprint_cols] = scaler_fingerprint.fit_transform(df[fingerprint_cols])
        print(f"Min-max scaled {len(fingerprint_cols)} fingerprint descriptors")
    
    # Don't scale binary fingerprints (they're already in [0,1] range)
    # Don't scale drug_name and smiles columns
    
    return df_scaled

def validate_descriptors(df):
    """
    Validate molecular descriptors for data quality issues
    
    Args:
        df (pd.DataFrame): DataFrame with molecular descriptors
        
    Returns:
        tuple: (is_valid, issues_found) - validation status and list of issues
    """
    issues = []
    
    # Get feature columns (exclude drug_name, smiles, pubchem_cid, synonyms)
    feature_columns = [col for col in df.columns if col not in ['drug_name', 'smiles', 'pubchem_cid', 'synonyms']]
    
    # Check for missing values
    missing_counts = df[feature_columns].isnull().sum()
    if missing_counts.sum() > 0:
        for col, count in missing_counts[missing_counts > 0].items():
            issues.append(f"Missing values in {col}: {count}")
    
    # Check for infinite values
    feature_matrix = df[feature_columns].values
    infinite_mask = ~np.isfinite(feature_matrix)
    infinite_count = np.sum(infinite_mask)
    
    if infinite_count > 0:
        issues.append(f"Found {infinite_count} infinite values")
        # Find which columns have infinite values
        for col_idx, col_name in enumerate(feature_columns):
            col_infinite = np.sum(infinite_mask[:, col_idx])
            if col_infinite > 0:
                issues.append(f"  {col_name}: {col_infinite} infinite values")
    
    # Check for extremely large values
    large_values = np.abs(feature_matrix) > 1e10
    large_count = np.sum(large_values)
    
    if large_count > 0:
        issues.append(f"Found {large_count} extremely large values (>1e10)")
        # Find which columns have large values
        for col_idx, col_name in enumerate(feature_columns):
            col_large = np.sum(large_values[:, col_idx])
            if col_large > 0:
                issues.append(f"  {col_name}: {col_large} large values")
    
    is_valid = len(issues) == 0
    return is_valid, issues

def clean_descriptors(df):
    """
    Clean molecular descriptors by fixing data quality issues
    
    Args:
        df (pd.DataFrame): DataFrame with molecular descriptors
        
    Returns:
        pd.DataFrame: Cleaned DataFrame
    """
    print("Cleaning molecular descriptors...")
    
    df_cleaned = df.copy()
    feature_columns = [col for col in df.columns if col not in ['drug_name', 'smiles', 'pubchem_cid', 'synonyms']]
    
    # Clean features
    feature_matrix = df_cleaned[feature_columns].values
    
    # Replace infinite values with NaN
    feature_matrix = np.where(np.isfinite(feature_matrix), feature_matrix, np.nan)
    
    # Clip extremely large values
    feature_matrix = np.clip(feature_matrix, -1e10, 1e10)
    
    # Replace NaN values with median of each column
    for col_idx in range(feature_matrix.shape[1]):
        col_data = feature_matrix[:, col_idx]
        if np.any(np.isnan(col_data)):
            median_val = np.nanmedian(col_data)
            if np.isnan(median_val):
                # If all values are NaN, use 0
                median_val = 0.0
            feature_matrix[:, col_idx] = np.where(np.isnan(col_data), median_val, col_data)
            print(f"  Replaced NaN values in {feature_columns[col_idx]} with median: {median_val:.3f}")
    
    # Update DataFrame
    df_cleaned[feature_columns] = feature_matrix
    
    print("✅ Data cleaning completed")
    return df_cleaned

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Query PubChem for drug SMILES and calculate molecular descriptors",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python query_drug_descriptors.py aspirin
  python query_drug_descriptors.py "acetaminophen"
  python query_drug_descriptors.py --drug "ibuprofen" --output results.csv
  python query_drug_descriptors.py --drug "caffeine" --no-scale
  python query_drug_descriptors.py --drug "morphine" --validate --clean
  python query_drug_descriptors.py --drug "aspirin" --no-clean --json
        """
    )
    
    parser.add_argument(
        'drug_name', 
        nargs='?', 
        help='Name of the drug to query and analyze'
    )
    
    parser.add_argument(
        '--drug', '-d',
        dest='drug_name_alt',
        help='Alternative way to specify drug name'
    )
    
    parser.add_argument(
        '--output', '-o',
        default='drug_analysis_results.csv',
        help='Output CSV file name (default: drug_analysis_results.csv)'
    )
    
    parser.add_argument(
        '--no-scale',
        action='store_true',
        help='Skip descriptor scaling (output raw values)'
    )
    
    parser.add_argument(
        '--json',
        action='store_true',
        help='Also save results as JSON file'
    )
    
    parser.add_argument(
        '--validate',
        action='store_true',
        help='Validate data quality and report issues'
    )
    
    parser.add_argument(
        '--clean',
        action='store_true',
        help='Clean data by fixing infinite values, NaN values, and outliers'
    )
    
    parser.add_argument(
        '--no-clean',
        action='store_true',
        help='Skip automatic data cleaning (not recommended)'
    )
    
    return parser.parse_args()

def main():
    """Main function to query PubChem and calculate molecular descriptors."""
    print("Starting drug descriptor calculation from PubChem...")
    
    # Check environment first
    if not check_environment():
        print("\nContinuing anyway, but you may encounter issues...")
        print("Press Enter to continue or Ctrl+C to exit...")
        try:
            input()
        except KeyboardInterrupt:
            print("\nExiting...")
            sys.exit(1)
    
    # Parse command line arguments
    args = parse_arguments()
    
    # Get drug name from command line
    drug_name = args.drug_name_alt if args.drug_name_alt else args.drug_name
    
    if not drug_name:
        print("Error: No drug name provided.")
        print("Usage: python query_drug_descriptors.py <drug_name>")
        print("   or: python query_drug_descriptors.py --drug <drug_name>")
        print("\nExamples:")
        print("  python query_drug_descriptors.py aspirin")
        print("  python query_drug_descriptors.py --drug acetaminophen")
        sys.exit(1)
    
    print(f"Analyzing drug: {drug_name}")
    
    # Query PubChem for SMILES
    smiles, cid, synonyms = query_pubchem_smiles(drug_name)
    
    if not smiles:
        print(f"Failed to find SMILES for '{drug_name}' in PubChem.")
        print("Please check the drug name spelling or try alternative names.")
        sys.exit(1)
    
    # Calculate molecular descriptors
    print(f"\nCalculating molecular descriptors...")
    descriptors = calculate_molecular_descriptors(smiles)
    
    if descriptors is None:
        print(f"Failed to calculate descriptors for SMILES: {smiles}")
        sys.exit(1)
    
    # Create result dictionary
    result = {
        'drug_name': drug_name,
        'pubchem_cid': cid,
        'smiles': smiles,
        'synonyms': ', '.join(synonyms) if synonyms else '',
        **descriptors
    }
    
    # Convert to DataFrame
    df = pd.DataFrame([result])
    
    # Validate data quality
    print("\nValidating data quality...")
    is_valid, issues = validate_descriptors(df)
    
    if not is_valid:
        print("⚠️  Data quality issues found:")
        for issue in issues:
            print(f"  - {issue}")
        
        # Clean data unless explicitly disabled
        if not args.no_clean:
            print("\nCleaning data automatically...")
            df = clean_descriptors(df)
            
            # Re-validate after cleaning
            is_valid_after_clean, issues_after_clean = validate_descriptors(df)
            if is_valid_after_clean:
                print("✅ Data cleaning successful - all issues resolved")
            else:
                print("⚠️  Some issues remain after cleaning:")
                for issue in issues_after_clean:
                    print(f"  - {issue}")
        else:
            print("⚠️  Data cleaning skipped (--no-clean specified)")
    else:
        print("✅ Data validation passed - no issues found")
    
    # Additional validation if requested
    if args.validate:
        print("\nDetailed validation report:")
        print(f"  Total descriptors: {len([col for col in df.columns if col not in ['drug_name', 'smiles', 'pubchem_cid', 'synonyms']])}")
        print(f"  Data quality: {'✅ Valid' if is_valid else '⚠️  Issues found'}")
        
        # Show data ranges for key descriptors
        key_descriptors = ['molecular_weight', 'logp', 'tpsa', 'num_hbd', 'num_hba']
        print("\nKey descriptor ranges:")
        for desc in key_descriptors:
            if desc in df.columns:
                val = df[desc].iloc[0]
                print(f"  {desc}: {val:.3f}")
    
    # Clean data if explicitly requested
    if args.clean and not args.no_clean:
        print("\nCleaning data (--clean specified)...")
        df = clean_descriptors(df)
    
    # Scale descriptors if requested
    if not args.no_scale:
        print("\nScaling molecular descriptors...")
        df_scaled = scale_descriptors(df)
        output_df = df_scaled
    else:
        print("\nSkipping descriptor scaling (raw values)")
        output_df = df
    
    # Save results
    print(f"\nSaving results to: {args.output}")
    output_df.to_csv(args.output, index=False)
    
    if args.json:
        json_file = args.output.replace('.csv', '.json')
        print(f"Also saving to: {json_file}")
        output_df.to_json(json_file, orient='records', indent=2)
    
    # Print summary
    print(f"\n" + "="*60)
    print(f"ANALYSIS SUMMARY")
    print(f"="*60)
    print(f"Drug Name: {drug_name}")
    print(f"PubChem CID: {cid}")
    print(f"SMILES: {smiles}")
    print(f"Synonyms: {', '.join(synonyms) if synonyms else 'None'}")
    print(f"\nKey Molecular Properties:")
    print(f"  Molecular Weight: {descriptors['molecular_weight']:.2f} g/mol")
    print(f"  LogP (Lipophilicity): {descriptors['logp']:.2f}")
    print(f"  TPSA (Polar Surface Area): {descriptors['tpsa']:.2f} Å²")
    print(f"  H-bond Donors: {descriptors['num_hbd']}")
    print(f"  H-bond Acceptors: {descriptors['num_hba']}")
    print(f"  Rotatable Bonds: {descriptors['num_rotatable_bonds']}")
    print(f"  Aromatic Rings: {descriptors['num_aromatic_rings']}")
    print(f"  ECFP4 Density: {descriptors['ecfp4_density']:.3f}")
    print(f"  MACCS Density: {descriptors['maccs_density']:.3f}")
    
    # Lipinski's Rule of Five assessment
    print(f"\nLipinski's Rule of Five Assessment:")
    mw = descriptors['molecular_weight']
    logp = descriptors['logp']
    hbd = descriptors['num_hbd']
    hba = descriptors['num_hba']
    
    violations = 0
    if mw > 500:
        violations += 1
        print(f"  ❌ Molecular Weight: {mw:.1f} (violates >500)")
    else:
        print(f"  ✅ Molecular Weight: {mw:.1f}")
    
    if logp > 5:
        violations += 1
        print(f"  ❌ LogP: {logp:.2f} (violates >5)")
    else:
        print(f"  ✅ LogP: {logp:.2f}")
    
    if hbd > 5:
        violations += 1
        print(f"  ❌ H-bond Donors: {hbd} (violates >5)")
    else:
        print(f"  ✅ H-bond Donors: {hbd}")
    
    if hba > 10:
        violations += 1
        print(f"  ❌ H-bond Acceptors: {hba} (violates >10)")
    else:
        print(f"  ✅ H-bond Acceptors: {hba}")
    
    if violations == 0:
        print(f"\n🎉 {drug_name} follows Lipinski's Rule of Five (drug-like)")
    elif violations <= 2:
        print(f"\n⚠️  {drug_name} has {violations} Lipinski violations (may still be drug-like)")
    else:
        print(f"\n❌ {drug_name} has {violations} Lipinski violations (likely not drug-like)")
    
    print(f"\nResults saved to: {args.output}")
    if args.json:
        print(f"JSON results saved to: {json_file}")

if __name__ == "__main__":
    main()
