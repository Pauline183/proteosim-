from proteosim.protein_digestion import enzyme_cleavage_patterns, digest_protein_collection, compute_sequence_coverage

def test_digest_protein_collection():
    # 1. Sehr einfache Test-Proteine mit klaren K/R-Positionen
    dummy_proteins = {
        # Schneidstellen mit Trypsin: A K R Q  -> nach K und nach R
        "P1": "AKRQ",
        # K gefolgt von P: hier DARF NICHT geschnitten werden
        "P2": "AKPRA",
    }

    # 2. Trypsin-Spaltmuster verwenden
    cleave_pattern = enzyme_cleavage_patterns["Trypsin"]

    # 3. Funktion aufrufen
    digested = digest_protein_collection(
        dummy_proteins,
        cleave_pattern=cleave_pattern,
        min_pep_len=1,
        max_pep_len=50,
    )
    return digested

    # 4. Erwartete Peptid-Listen prüfen

    # "AKRQ"   -> Trypsin schneidet hinter K und R:
    #  "AK" | "R" | "Q"
    assert digested["P1"] == ["AK", "R", "Q"]

    # "AKPRA" -> K vor P: HIER KEIN SCHNITT; aber hinter R wird geschnitten:
    #  "AKP" | "RA"
    assert digested["P2"] == ["AKP", "RA"]

def test_compute_sequence_coverage():
 # Dummy protein sequence (10 amino acids)
    dummy_prot_seq = "ABCDEFGHIJ"
    
    # Dummy peptides covering known positions
    dummy_peps = ["ABC", "FGH"]
    
    # Call the function
    coverage = compute_sequence_coverage(dummy_prot_seq, dummy_peps)
    
    # ABC covers positions 0,1,2
    # FGH covers positions 5,6,7
    # → total 6 covered out of 10 → 60%
    expected_coverage = 60.0
    
    # Check if the result matches exactly
    assert coverage == expected_coverage, f"Expected {expected_coverage}, got {coverage}"
    return coverage



