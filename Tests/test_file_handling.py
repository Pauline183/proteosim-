from Proteosim.file_handling import read_fasta
def test_read_fasta():
    tmp_fasta_path = 'data/dummy_proteins.fasta'
    protein_map = read_fasta(tmp_fasta_path)
    print(protein_map)
    
    # Replace the strings with your fasta content
    # which you expect to be now available as a dictionary
    assert protein_map["A123"] == "DUMMYSEQUENCEONE"
    assert protein_map["B123"] == "DUMMYSEQUENCETWO"
    

