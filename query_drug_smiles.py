#!/usr/bin/env python3
"""
Simple drug SMILES query script.

Takes a drug name as input and outputs the SMILES string from PubChem.
This is a simplified version focused only on SMILES retrieval.

Usage:
    python query_drug_smiles.py aspirin
    python query_drug_smiles.py "acetaminophen"
"""

import sys
import argparse

# Try to import PubChemPy for better compound searching
try:
    import pubchempy as pcp
    HAS_PUBCHEMPY = True
except ImportError:
    HAS_PUBCHEMPY = False
    import requests

def query_smiles_pubchempy(drug_name):
    """
    Query PubChem for SMILES string using PubChemPy library.
    
    Args:
        drug_name (str): Name of the drug to search for
        
    Returns:
        str: SMILES string or None if not found
    """
    try:
        # Search for compounds by name
        compounds = list(pcp.get_compounds(drug_name, 'name'))
        
        if not compounds:
            # Try searching by synonym
            compounds = list(pcp.get_compounds(drug_name, 'synonym'))
        
        if not compounds:
            # Try fuzzy search
            compounds = list(pcp.get_compounds(drug_name, 'name', fuzzy=True))
        
        if compounds:
            # Use the first compound found
            compound = compounds[0]
            
            # Get SMILES (use the new smiles property, fall back to canonical_smiles)
            smiles = compound.smiles or compound.canonical_smiles
            
            if smiles:
                return smiles
        
        return None
        
    except Exception as e:
        print(f"Error querying PubChem: {e}", file=sys.stderr)
        return None

def query_smiles_http(drug_name):
    """
    Fallback method using HTTP requests when PubChemPy is not available.
    
    Args:
        drug_name (str): Name of the drug to search for
        
    Returns:
        str: SMILES string or None if not found
    """
    try:
        # Simple HTTP request to PubChem API
        search_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{requests.utils.quote(drug_name)}/property/IsomericSMILES/JSON"
        response = requests.get(search_url, timeout=10)
        
        if response.status_code == 200:
            data = response.json()
            if 'PropertyTable' in data and 'Properties' in data['PropertyTable']:
                properties = data['PropertyTable']['Properties'][0]
                smiles = properties.get('IsomericSMILES', '')
                return smiles if smiles else None
        
        return None
        
    except Exception as e:
        print(f"Error querying PubChem: {e}", file=sys.stderr)
        return None

def query_drug_smiles(drug_name):
    """
    Query PubChem for SMILES string of a drug.
    
    Args:
        drug_name (str): Name of the drug to search for
        
    Returns:
        str: SMILES string or None if not found
    """
    if HAS_PUBCHEMPY:
        return query_smiles_pubchempy(drug_name)
    else:
        return query_smiles_http(drug_name)

def main():
    """Main function."""
    parser = argparse.ArgumentParser(
        description="Query PubChem for drug SMILES string",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python query_drug_smiles.py aspirin
  python query_drug_smiles.py "acetaminophen"
  python query_drug_smiles.py ibuprofen
        """
    )
    
    parser.add_argument(
        'drug_name',
        help='Name of the drug to search for'
    )
    
    parser.add_argument(
        '--verbose', '-v',
        action='store_true',
        help='Show verbose output including search method'
    )
    
    args = parser.parse_args()
    
    if args.verbose:
        method = "PubChemPy" if HAS_PUBCHEMPY else "HTTP requests"
        print(f"Searching for '{args.drug_name}' using {method}...", file=sys.stderr)
    
    # Query for SMILES
    smiles = query_drug_smiles(args.drug_name)
    
    if smiles:
        print(smiles)
        return 0
    else:
        print(f"Could not find SMILES for '{args.drug_name}'", file=sys.stderr)
        return 1

if __name__ == "__main__":
    sys.exit(main())
