import re

def digest_protein_sequence(protein_seq, cleave_pattern):
    
    """
    Splits a protein sequence into peptides using a given cleavage pattern.

    Parameters
    ----------
    protein_seq : str
        The amino-acid sequence of a single protein. This should be a continuous
        one-letter amino-acid string without whitespace.
    cleave_pattern : str
        A regular expression (regex) pattern defining the cleavage rule.
        The function uses `re.split()` with this pattern to digest the protein.

    Returns
    -------
    list of str
        A list containing all peptide fragments produced by applying the
        cleavage pattern to the input sequence.

    """
    peptides = re.split(cleave_pattern, protein_seq)
    return peptides
    

def digest_protein_collection(protein_map, cleave_pattern, min_pep_len=5, max_pep_len=30):
      
    """
    Digests an entire collection of protein sequences using a given cleavage
    pattern and optional peptide length filters.

    Parameters
    ----------
    protein_map : dict
        A dictionary mapping protein IDs (str) to amino-acid sequences (str).
        Example: {"A123": "MKWVTFISLL...", "B122": "GHSKST..."}
    cleave_pattern : str
        A regular expression (regex) defining the cleavage rule.
        Example for trypsin: ``r'(?<=[KR])'``.
    min_pep_len : int, optional
        Minimum peptide length to keep. Default is 5.
    max_pep_len : int, optional
        Maximum peptide length to keep. Default is 30.

    Returns
    -------
    dict
        A dictionary mapping protein IDs to lists of digested peptide sequences.
        Example:
            {
                "A123": ["MKWT", "FISLL", ...],
                "B122": ["GH", "SKSTL", ...]
            }

    """
    digested = {}

    for protein_id, sequence in protein_map.items():
    
        peptides = re.split(cleave_pattern, sequence)

            # Filter peptides by length
        filtered = [pep for pep in peptides if min_pep_len <= len(pep) <= max_pep_len]
        digested[protein_id] = filtered
    return digested



def compute_sequence_coverage(protein_seq, peptides):
    """
    Compute percentage of residues in the protein that are covered
    by at least one peptide.
    """
    covered = set()

    for pep in peptides:
        start = protein_seq.find(pep)
        if start != -1:
            for i in range(start, start + len(pep)):
                covered.add(i)

    return len(covered) / len(protein_seq) * 100
 
enzyme_cleavage_patterns = {
    'LysC': r'(?<=K)',
    'LysN': r'(?=K)',
    'ArgC': r'(?<=R)',
    'Trypsin': r'((?<=[KR])(?!P))'
    
}