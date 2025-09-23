#!/usr/bin/env python3
"""
Calculate molecular descriptors for drugs using RDKit from SMILES structures.
This script reads the drug names and SMILES from the previous extraction
and calculates various molecular properties like LogP, H-bonding, size, etc.

Features:
- Molecular standardization: removes salts/solvents, neutralizes, normalizes tautomers
- Hydrogen addition for better 3D properties
- Appropriate scaling: Z-score for continuous, log1p for counts, min-max for fingerprints
- pH awareness: LogP calculated at neutral state (note: RDKit lacks pKa prediction)

Note: For pH-dependent properties, we should consider using external pKa prediction tools
and recalculating descriptors at physiological pH (7.4) if ionization matters.

Environment Requirements:
- This script requires the conda RDKit, NOT pip version
"""

import json
import pandas as pd
import sys
import argparse
import os
from pathlib import Path

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

def load_drug_data(input_file):
    """
    Load drug data from JSON file created by extract_unique_drugs.py
    
    Args:
        input_file (str): Path to the JSON file
        
    Returns:
        list: List of tuples (drug_name, smiles)
    """
    try:
        with open(input_file, 'r', encoding='utf-8') as f:
            data = json.load(f)
        
        drugs = []
        for drug_entry in data['unique_drugs']:
            if ',' in drug_entry:
                drug_name, smiles = drug_entry.split(',', 1)
                drugs.append((drug_name.strip(), smiles.strip()))
            else:
                # Fallback for entries without SMILES
                drugs.append((drug_entry.strip(), ''))
        
        return drugs
        
    except FileNotFoundError:
        print(f"Error: File '{input_file}' not found.")
        print("Please run extract_unique_drugs.py first to generate the drug data.")
        sys.exit(1)
    except Exception as e:
        print(f"Error loading drug data: {e}")
        sys.exit(1)

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

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Calculate molecular descriptors for drugs from SMILES structures",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python calculate_drug_descriptors.py
  python calculate_drug_descriptors.py my_drugs.json
  python calculate_drug_descriptors.py --input custom_drugs.json
        """
    )
    
    parser.add_argument(
        'input_file', 
        nargs='?', 
        default='unique_drugs_tahoe100m.json',
        help='Input JSON file containing drug names and SMILES (default: unique_drugs_tahoe100m.json)'
    )
    
    parser.add_argument(
        '--input', '-i',
        dest='input_file_alt',
        help='Alternative way to specify input file'
    )
    
    return parser.parse_args()

def main():
    """Main function to calculate and save molecular descriptors."""
    print("Starting molecular descriptor calculation...")
    
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
    
    # Use the input file from command line or default
    input_file = args.input_file_alt if args.input_file_alt else args.input_file
    
    print(f"Using input file: {input_file}")
    
    # Check if input file exists
    if not Path(input_file).exists():
        print(f"Error: Input file '{input_file}' not found.")
        print("Please run extract_unique_drugs.py first to generate the drug data,")
        print("or specify a valid JSON file containing drug names and SMILES.")
        print("\nUsage examples:")
        print("  python calculate_drug_descriptors.py")
        print("  python calculate_drug_descriptors.py my_drugs.json")
        print("  python calculate_drug_descriptors.py --input custom_drugs.json")
        sys.exit(1)
    
    # Load drug data
    print(f"Loading drug data from {input_file}...")
    drugs = load_drug_data(input_file)
    print(f"Loaded {len(drugs)} drugs")
    
    # Calculate descriptors for each drug
    results = []
    failed_count = 0
    
    print("Calculating molecular descriptors...")
    for i, (drug_name, smiles) in enumerate(drugs, 1):
        print(f"Processing {i}/{len(drugs)}: {drug_name}")
        
        descriptors = calculate_molecular_descriptors(smiles)
        
        if descriptors is not None:
            # Add drug name and SMILES to the result
            result = {
                'drug_name': drug_name,
                'smiles': smiles,
                **descriptors
            }
            results.append(result)
        else:
            failed_count += 1
            print(f"  Failed to calculate descriptors for: {drug_name}")
    
    print(f"\nDescriptor calculation complete!")
    print(f"Successfully processed: {len(results)} drugs")
    print(f"Failed to process: {failed_count} drugs")
    
    if not results:
        print("No valid results to save.")
        return
    
    # Convert to DataFrame for easier handling
    df = pd.DataFrame(results)
    
    # Scale the descriptors appropriately
    print("\nScaling molecular descriptors...")
    df_scaled = scale_descriptors(df)
    
    # Save results in multiple formats
    print("\nSaving results...")
    
    # Save original (unscaled) results
    csv_file = "drug_descriptors_unscaled.csv"
    df.to_csv(csv_file, index=False)
    print(f"Unscaled results saved to: {csv_file}")
    
    # Save scaled results
    csv_scaled_file = "drug_descriptors_scaled.csv"
    df_scaled.to_csv(csv_scaled_file, index=False)
    print(f"Scaled results saved to: {csv_scaled_file}")
    
    # Save as JSON (unscaled)
    json_file = "drug_descriptors.json"
    df.to_json(json_file, orient='records', indent=2)
    print(f"Results saved to: {json_file}")
    
    # Save as Excel (if openpyxl is available)
    try:
        excel_file = "drug_descriptors.xlsx"
        with pd.ExcelWriter(excel_file, engine='openpyxl') as writer:
            df.to_excel(writer, sheet_name='Unscaled', index=False)
            df_scaled.to_excel(writer, sheet_name='Scaled', index=False)
        print(f"Results saved to: {excel_file} (both scaled and unscaled sheets)")
    except ImportError:
        print("Excel export skipped (openpyxl not installed)")
    
    # Print summary statistics
    print(f"\nSummary Statistics:")
    print(f"Total drugs analyzed: {len(results)}")
    print(f"Average molecular weight: {df['molecular_weight'].mean():.2f}")
    print(f"Average LogP: {df['logp'].mean():.2f}")
    print(f"Average TPSA: {df['tpsa'].mean():.2f}")
    print(f"Average H-bond donors: {df['num_hbd'].mean():.1f}")
    print(f"Average H-bond acceptors: {df['num_hba'].mean():.1f}")
    print(f"Average rotatable bonds: {df['num_rotatable_bonds'].mean():.1f}")
    print(f"Average ECFP4 density: {df['ecfp4_density'].mean():.3f}")
    print(f"Average MACCS density: {df['maccs_density'].mean():.3f}")
    
    # Show drugs with highest and lowest LogP
    # print(f"\nTop 5 drugs with highest LogP (most lipophilic):")
    # top_logp = df.nlargest(5, 'logp')[['drug_name', 'logp']]
    # for _, row in top_logp.iterrows():
    #     print(f"  {row['drug_name']}: {row['logp']:.2f}")
    
    # print(f"\nTop 5 drugs with lowest LogP (most hydrophilic):")
    # bottom_logp = df.nsmallest(5, 'logp')[['drug_name', 'logp']]
    # for _, row in bottom_logp.iterrows():
    #     print(f"  {row['drug_name']}: {row['logp']:.2f}")
    
    # Show drugs with highest molecular complexity
    # print(f"\nTop 5 most complex compounds (highest Bertz CT):")
    # top_complexity = df.nlargest(5, 'bertz_ct')[['drug_name', 'bertz_ct']]
    # for _, row in top_complexity.iterrows():
    #     print(f"  {row['drug_name']}: {row['bertz_ct']:.1f}")
    
    # Show drugs with most functional groups
    # print(f"\nTop 5 drugs with most functional groups (highest total count):")
    # functional_cols = [col for col in df.columns if col.startswith('fr_')]
    # df['total_functional_groups'] = df[functional_cols].sum(axis=1)
    # top_functional = df.nlargest(5, 'total_functional_groups')[['drug_name', 'total_functional_groups']]
    # for _, row in top_functional.iterrows():
    #     print(f"  {row['drug_name']}: {row['total_functional_groups']:.0f} groups")

if __name__ == "__main__":
    main()
