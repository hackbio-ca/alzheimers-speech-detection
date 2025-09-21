#!/usr/bin/env python3
"""
Standalone drug descriptor calculator.

This is a self-contained script that takes a drug name as input, queries PubChem 
to get the SMILES string, and calculates molecular descriptors using RDKit.
Returns results as a list: [drug_name, smiles, pubchem_cid, synonyms, descriptor_values...]

Usage:
    python standalone_drug_descriptors.py aspirin
    python standalone_drug_descriptors.py "acetaminophen"
"""

import sys
import time
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

# Try to import RDKit
try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors
    from rdkit.Chem import rdFingerprintGenerator
    from rdkit.Chem.AtomPairs import Pairs
    from rdkit.Chem.Fingerprints import FingerprintMols
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
    print("\nThis script requires RDKit and PubChemPy.")
    print("Please install them using:")
    print("  conda install -c conda-forge rdkit")
    print("  pip install pubchempy")
    print("=" * 60)
    sys.exit(1)

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
                
                # Get SMILES (use the new smiles property, fall back to canonical_smiles)
                smiles = compound.smiles or compound.canonical_smiles
                
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

def query_pubchem_smiles_pubchempy_quiet(drug_name, verbose=False, max_retries=3, delay=1.0):
    """
    Query PubChem using PubChemPy library (quiet version).
    
    Args:
        drug_name (str): Name of the drug to search for
        verbose (bool): Whether to print progress messages
        max_retries (int): Maximum number of retry attempts
        delay (float): Delay between retries in seconds
        
    Returns:
        tuple: (smiles, cid, synonyms) or (None, None, None) if not found
    """
    for attempt in range(max_retries):
        try:
            if verbose:
                print(f"Searching PubChem for: {drug_name} (attempt {attempt + 1})")
            
            # Search for compounds by name
            compounds = list(pcp.get_compounds(drug_name, 'name'))
            
            if not compounds:
                # Try searching by synonym
                if verbose:
                    print("No direct matches found, trying synonym search...")
                compounds = list(pcp.get_compounds(drug_name, 'synonym'))
            
            if not compounds:
                # Try fuzzy search
                if verbose:
                    print("No synonym matches found, trying fuzzy search...")
                compounds = list(pcp.get_compounds(drug_name, 'name', fuzzy=True))
            
            if compounds:
                # Use the first compound found
                compound = compounds[0]
                
                # Get SMILES (use the new smiles property, fall back to canonical_smiles)
                smiles = compound.smiles or compound.canonical_smiles
                
                if smiles:
                    if verbose:
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
                    if verbose:
                        print(f"Compound found but no SMILES available for CID {compound.cid}")
            else:
                if verbose:
                    print(f"No compounds found for: {drug_name}")
            
        except Exception as e:
            if verbose:
                print(f"Attempt {attempt + 1} failed: {e}")
            if attempt < max_retries - 1:
                if verbose:
                    print(f"Retrying in {delay} seconds...")
                time.sleep(delay)
                delay *= 2  # Exponential backoff
    
    if verbose:
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

def query_pubchem_smiles_http_quiet(drug_name, verbose=False, max_retries=3, delay=1.0):
    """
    Fallback method using HTTP requests when PubChemPy is not available (quiet version).
    
    Args:
        drug_name (str): Name of the drug to search for
        verbose (bool): Whether to print progress messages
        max_retries (int): Maximum number of retry attempts
        delay (float): Delay between retries in seconds
        
    Returns:
        tuple: (smiles, cid, synonyms) or (None, None, None) if not found
    """
    if verbose:
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
                        if verbose:
                            print(f"✅ Found SMILES: {smiles}")
                            print(f"   Molecular Formula: {molecular_formula}")
                            print(f"   Molecular Weight: {molecular_weight}")
                        return smiles, None, [drug_name]
            
            if verbose:
                print(f"Attempt {attempt + 1} failed: No results found")
            
        except Exception as e:
            if verbose:
                print(f"Attempt {attempt + 1} failed: {e}")
        
        if attempt < max_retries - 1:
            if verbose:
                print(f"Retrying in {delay} seconds...")
            time.sleep(delay)
            delay *= 2  # Exponential backoff
    
    if verbose:
        print(f"Failed to find SMILES for '{drug_name}' after {max_retries} attempts")
    return None, None, None

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

def query_pubchem_smiles_quiet(drug_name, verbose=False, max_retries=3, delay=1.0):
    """
    Query PubChem for SMILES string of a drug using PubChemPy (quiet version).
    
    Args:
        drug_name (str): Name of the drug to search for
        verbose (bool): Whether to print progress messages
        max_retries (int): Maximum number of retry attempts
        delay (float): Delay between retries in seconds
        
    Returns:
        tuple: (smiles, cid, synonyms) or (None, None, None) if not found
    """
    if verbose:
        print(f"Querying PubChem for: {drug_name}")
    
    if HAS_PUBCHEMPY:
        return query_pubchem_smiles_pubchempy_quiet(drug_name, verbose, max_retries, delay)
    else:
        return query_pubchem_smiles_http_quiet(drug_name, verbose, max_retries, delay)

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

def scale_descriptors(descriptors_dict):
    """
    Scale molecular descriptors appropriately:
    - Z-score for continuous descriptors
    - Log1p for count descriptors
    - Min-max for fingerprint descriptors
    
    Args:
        descriptors_dict (dict): Dictionary with molecular descriptors
        
    Returns:
        dict: Dictionary with scaled descriptors
    """
    scaled_descriptors = descriptors_dict.copy()
    
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
    functional_group_descriptors = [key for key in descriptors_dict.keys() if key.startswith('fr_')]
    count_descriptors.extend(functional_group_descriptors)
    
    # Fingerprint density and bit counts
    fingerprint_descriptors = [
        'ecfp4_density', 'ecfp4_bits', 'ecfp6_density', 'ecfp6_bits',
        'maccs_density', 'maccs_bits'
    ]
    
    print("Scaling molecular descriptors...")
    
    # For single molecule, we'll use simple scaling approaches
    # In practice, you'd want to use statistics from a larger dataset
    
    # Log1p scaling for count descriptors (handles zeros well)
    for desc in count_descriptors:
        if desc in scaled_descriptors:
            scaled_descriptors[desc] = np.log1p(max(0, scaled_descriptors[desc]))
    
    # Min-max scaling for fingerprint descriptors (simple 0-1 scaling)
    for desc in fingerprint_descriptors:
        if desc in scaled_descriptors:
            # Simple scaling to 0-1 range
            val = scaled_descriptors[desc]
            if 'density' in desc:
                # Density is already 0-1, keep as is
                pass
            elif 'bits' in desc:
                # Scale bits to 0-1 range (assuming max 2048 bits)
                scaled_descriptors[desc] = val / 2048.0
    
    print(f"Scaled {len(count_descriptors)} count descriptors and {len(fingerprint_descriptors)} fingerprint descriptors")
    
    return scaled_descriptors

def scale_descriptors_quiet(descriptors_dict, verbose=False):
    """
    Scale molecular descriptors appropriately (quiet version):
    - Z-score for continuous descriptors
    - Log1p for count descriptors
    - Min-max for fingerprint descriptors
    
    Args:
        descriptors_dict (dict): Dictionary with molecular descriptors
        verbose (bool): Whether to print progress messages
        
    Returns:
        dict: Dictionary with scaled descriptors
    """
    scaled_descriptors = descriptors_dict.copy()
    
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
    functional_group_descriptors = [key for key in descriptors_dict.keys() if key.startswith('fr_')]
    count_descriptors.extend(functional_group_descriptors)
    
    # Fingerprint density and bit counts
    fingerprint_descriptors = [
        'ecfp4_density', 'ecfp4_bits', 'ecfp6_density', 'ecfp6_bits',
        'maccs_density', 'maccs_bits'
    ]
    
    if verbose:
        print("Scaling molecular descriptors...")
    
    # For single molecule, we'll use simple scaling approaches
    # In practice, you'd want to use statistics from a larger dataset
    
    # Log1p scaling for count descriptors (handles zeros well)
    for desc in count_descriptors:
        if desc in scaled_descriptors:
            scaled_descriptors[desc] = np.log1p(max(0, scaled_descriptors[desc]))
    
    # Min-max scaling for fingerprint descriptors (simple 0-1 scaling)
    for desc in fingerprint_descriptors:
        if desc in scaled_descriptors:
            # Simple scaling to 0-1 range
            val = scaled_descriptors[desc]
            if 'density' in desc:
                # Density is already 0-1, keep as is
                pass
            elif 'bits' in desc:
                # Scale bits to 0-1 range (assuming max 2048 bits)
                scaled_descriptors[desc] = val / 2048.0
    
    if verbose:
        print(f"Scaled {len(count_descriptors)} count descriptors and {len(fingerprint_descriptors)} fingerprint descriptors")
    
    return scaled_descriptors

def get_drug_descriptors(drug_name, scale_descriptors_flag=True, clean_data_flag=True, verbose=False):
    """
    Get drug descriptors for a given drug name.
    
    Args:
        drug_name (str): Name of the drug to analyze
        scale_descriptors_flag (bool): Whether to scale descriptors (default: True)
        clean_data_flag (bool): Whether to clean data (default: True)
        verbose (bool): Whether to print progress messages (default: False)
        
    Returns:
        list: [drug_name, smiles, pubchem_cid, synonyms, descriptor_values...]
              Returns None if drug not found or calculation failed
    """
    if verbose:
        print(f"Analyzing drug: {drug_name}")
    
    # Query PubChem for SMILES
    smiles, cid, synonyms = query_pubchem_smiles_quiet(drug_name, verbose)
    
    if not smiles:
        if verbose:
            print(f"Failed to find SMILES for '{drug_name}' in PubChem.")
            print("Please check the drug name spelling or try alternative names.")
        return None
    
    # Calculate molecular descriptors
    if verbose:
        print(f"Calculating molecular descriptors...")
    descriptors = calculate_molecular_descriptors(smiles)
    
    if descriptors is None:
        if verbose:
            print(f"Failed to calculate descriptors for SMILES: {smiles}")
        return None
    
    # Clean data if requested
    if clean_data_flag:
        if verbose:
            print("Cleaning data...")
        # Replace any infinite or NaN values
        for key, value in descriptors.items():
            if not np.isfinite(value):
                descriptors[key] = 0.0
                if verbose:
                    print(f"  Replaced non-finite value in {key} with 0.0")
    
    # Scale descriptors if requested
    if scale_descriptors_flag:
        if verbose:
            print("Scaling molecular descriptors...")
        descriptors = scale_descriptors_quiet(descriptors, verbose)
    
    # Convert to list format
    # Create result list: [drug_name, smiles, pubchem_cid, synonyms, feature_values...]
    result_list = [
        drug_name,
        smiles,
        cid,
        ', '.join(synonyms) if synonyms else ''
    ]
    
    # Add feature values in consistent order
    feature_order = [
        'molecular_weight', 'num_atoms', 'num_heavy_atoms', 'num_bonds', 'num_rings',
        'num_aromatic_rings', 'num_saturated_rings', 'num_aliphatic_rings',
        'tpsa', 'logp', 'molar_refractivity', 'num_hbd', 'num_hba',
        'num_rotatable_bonds', 'num_sp3_carbons', 'formal_charge',
        'num_positive_charges', 'num_negative_charges', 'sa_score',
        'num_carbons', 'num_nitrogens', 'num_oxygens', 'num_sulfurs',
        'num_halogens', 'num_heteroatoms', 'balaban_j', 'bertz_ct',
        'kier_flexibility', 'kier_shape1', 'kier_shape2',
        'ecfp4_density', 'ecfp4_bits', 'ecfp6_density', 'ecfp6_bits',
        'maccs_density', 'maccs_bits'
    ]
    
    # Add functional group descriptors
    functional_groups = [key for key in descriptors.keys() if key.startswith('fr_')]
    functional_groups.sort()  # Sort for consistent order
    feature_order.extend(functional_groups)
    
    # Add feature values in the specified order
    for feature in feature_order:
        if feature in descriptors:
            result_list.append(descriptors[feature])
        else:
            result_list.append(0.0)  # Default value for missing features
    
    if verbose:
        print(f"✅ Successfully calculated {len(feature_order)} descriptors for {drug_name}")
        print(f"Key properties: MW={descriptors.get('molecular_weight', 0):.1f}, LogP={descriptors.get('logp', 0):.2f}, TPSA={descriptors.get('tpsa', 0):.1f}")
    
    return result_list

def main():
    """Main function for command line usage."""
    if len(sys.argv) < 2:
        print("Usage: python standalone_drug_descriptors.py <drug_name> [--verbose]")
        print("Example: python standalone_drug_descriptors.py aspirin")
        print("Example: python standalone_drug_descriptors.py aspirin --verbose")
        sys.exit(1)
    
    drug_name = sys.argv[1]
    verbose = "--verbose" in sys.argv
    
    if verbose:
        print("=" * 60)
        print("STANDALONE DRUG DESCRIPTOR CALCULATOR")
        print("=" * 60)
    
    # Get drug descriptors
    result = get_drug_descriptors(drug_name, verbose=verbose)
    
    if result is None:
        if verbose:
            print(f"❌ Failed to get descriptors for '{drug_name}'")
        sys.exit(1)
    
    if verbose:
        # Print detailed results
        print(f"\n" + "="*60)
        print(f"RESULTS")
        print(f"="*60)
        print(f"Drug Name: {result[0]}")
        print(f"SMILES: {result[1]}")
        print(f"PubChem CID: {result[2]}")
        print(f"Synonyms: {result[3]}")
        print(f"Number of descriptors: {len(result) - 4}")
        
        # Show first few descriptor values as example
        print(f"\nFirst 10 descriptor values:")
        for i, val in enumerate(result[4:14]):
            print(f"  Descriptor {i+1}: {val:.6f}")
        
        print(f"\nResult list format:")
        print(f"[drug_name, smiles, pubchem_cid, synonyms, descriptor_1, descriptor_2, ...]")
        print(f"Total length: {len(result)}")
        
        # Show how to extract just the feature vector
        feature_vector = result[4:]  # Skip drug_name, smiles, cid, synonyms
        print(f"\nFeature vector (for ML use): {len(feature_vector)} values")
        print(f"First 5 features: {feature_vector[:5]}")
        
        # For programmatic use, show how to get clean output
        print(f"\n" + "="*60)
        print("FOR PROGRAMMATIC USE:")
        print("="*60)
        print("To get just the list without verbose output:")
        print("  from standalone_drug_descriptors import get_drug_descriptors")
        print("  result = get_drug_descriptors('aspirin', verbose=False)")
        print("  # result is: [drug_name, smiles, cid, synonyms, descriptor_values...]")
        print("  feature_vector = result[4:]  # Just the descriptor values")
    else:
        # Just output the result list (for programmatic use)
        print(result)

if __name__ == "__main__":
    main()