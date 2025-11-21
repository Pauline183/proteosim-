def read_fasta(filepath):
    """
    Reads a FASTA file and returns a dictionary mapping protein IDs to sequences.

    Parameters
    ----------
    filepath : str
        Path to the FASTA file to be read. The file must contain header lines
        starting with '>' and protein sequences on the following lines.

    Returns
    -------
    dict
        A dictionary where:
            - keys   = extracted protein IDs (strings)
            - values = corresponding protein sequences as continuous strings

    Notes
    -----
    The function assumes that FASTA headers follow the UniProt-style format:
        >sp|P12345|...
    In this case, the protein ID is extracted as the middle field (e.g., "P12345").
    """
    protein_map = {}
    current_id = None
    current_sequence = []
    with open(filepath, 'r', encoding='utf-8') as fasta_handle:
        for line in fasta_handle:
            print(line)
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith('>'):
                if current_id is not None:
                    protein_map[current_id] = ''.join(current_sequence)
                    current_sequence = []
                current_id = stripped.split('|')[1]
                print(current_id)
            else:
                current_sequence.append(stripped)
    if current_id is not None:
        protein_map[current_id] = ''.join(current_sequence)
    return protein_map
