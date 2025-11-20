def test_read_fasta():
    tmp_fasta_path = '/Users/paulinekipp/Documents/GitHub/advanced-python_course/data/dummy_proteins.fasta'
    protein_map = read_fasta(tmp_fasta_path)
    print(protein_map)
    
    # Replace the strings with your fasta content
    # which you expect to be now available as a dictionary
    assert protein_map["A123"] == "DUMMYSEQUENCEONE"
    assert protein_map["B123"] == "DUMMYSEQUENCETWO"

test_read_fasta()