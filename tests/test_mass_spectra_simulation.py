import proteosim as ps
import matplotlib.pyplot as plt
import numpy as np

def test_calculate_mol_mass():
    aa_mass_dict = {
    'A': 71.08, 'R': 156.19, 'N': 114.10, 'D': 115.09,
    'C': 103.15, 'E': 129.12, 'Q': 128.13, 'G': 57.05,
    'H': 137.14, 'I': 113.16, 'L': 113.16, 'K': 128.17,
    'M': 131.19, 'F': 147.18, 'P': 97.12, 'S': 87.08,
    'T': 101.11, 'W': 186.21, 'Y': 163.18, 'V': 99.13,
}
    test_peptide = 'AFG'
    expected_mass = {'AFG': 275.31}
    calculated_mass = ps.calculate_mol_mass(test_peptide, aa_mass_dict)
    assert calculated_mass == expected_mass

test_calculate_mol_mass()

def test_calculate_mol_mass_collection():
    aa_mass_dict = {
    'A': 71.08, 'R': 156.19, 'N': 114.10, 'D': 115.09,
    'C': 103.15, 'E': 129.12, 'Q': 128.13, 'G': 57.05,
    'H': 137.14, 'I': 113.16, 'L': 113.16, 'K': 128.17,
    'M': 131.19, 'F': 147.18, 'P': 97.12, 'S': 87.08,
    'T': 101.11, 'W': 186.21, 'Y': 163.18, 'V': 99.13,
}
    peptides = {'AFG', 'KLM'}
    expected = {'AFG': 275.31,'KLM': 372.52}
    actual = ps.calculate_mol_mass_collection(peptides,aa_mass_dict)

    assert actual == expected

test_calculate_mol_mass_collection()

def test_calculate_mz_collection():
    peptide_mass_map = {
        "A": 71.08,
        "AF": 71.08 + 147.18,   
    }
    actual = ps.calculate_mz_collection(peptide_mass_map, charge=2)
    expected = {"A": (71.08 + 2 * 1.007) / 2,
        "AF": (218.26 + 2 * 1.007) / 2}

    assert actual == expected

test_calculate_mz_collection()

def test_fragment_peptide():
    peptide = 'PEPT'
    expected = ['P', 'PE', 'PEP','T', 'PT', 'EPT']
    actual = ps.fragment_peptide(peptide)

    assert set(actual) == set(expected)

test_fragment_peptide()