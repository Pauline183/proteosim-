amino_acid_mass_dalton = {
    'A': 71.08, 'R': 156.19, 'N': 114.10, 'D': 115.09,
    'C': 103.15, 'E': 129.12, 'Q': 128.13, 'G': 57.05,
    'H': 137.14, 'I': 113.16, 'L': 113.16, 'K': 128.17,
    'M': 131.19, 'F': 147.18, 'P': 97.12, 'S': 87.08,
    'T': 101.11, 'W': 186.21, 'Y': 163.18, 'V': 99.13,
}
def calculate_mol_mass(peptide_seq, amino_acid_mass_dict):
    """
    Calculate the molecular mass of a peptide sequence.

    This function computes the total molecular mass of a peptide by summing
    the masses of all amino acids in the sequence. The amino-acid masses must
    be provided as a dictionary mapping one-letter amino-acid codes to their
    respective monoisotopic or average masses.

    Parameters
    ----------
    peptide_seq : str
        The peptide sequence (one-letter amino-acid codes).
    amino_acid_mass_dict : dict[str, float]
        A dictionary mapping amino-acid single-letter codes to their masses.

    Returns
    -------
    dict[str, float]
        A dictionary containing one entry where the key is the peptide sequence
        and the value is its total molecular mass.

    """
    mass = 0.0
    for aa in peptide_seq:
        mass += amino_acid_mass_dict[aa]

    # Leeres dict erstellen
    result = {}

    # Eintrag hinzufügen
    result[peptide_seq] = mass

    # Zurückgeben
    return result

def calculate_mol_mass_collection(peptides, amino_acid_mass_dict):
    """
    Calculate the molecular mass of a peptide sequence.

    This function computes the total molecular mass of a peptide by summing
    the masses of all amino acids in the sequence. The amino-acid masses must
    be provided as a dictionary mapping one-letter amino-acid codes to their
    respective monoisotopic or average masses.

    Parameters
    ----------
    peptide_seq : str
        The peptide sequence (one-letter amino-acid codes).
    amino_acid_mass_dict : dict[str, float]
        A dictionary mapping amino-acid single-letter codes to their masses.

    Returns
    -------
    dict[str, float]
    A dictionary containing one entry where the key is the peptide sequence
    and the value is its total molecular mass.
    """
    result = {}

    for pep in peptides:
        mass = 0.0
        for aa in pep:
            mass += amino_acid_mass_dict[aa]
        result[pep] = mass

    return result

def calculate_mz_collection(peptide_mass_map, charge=2, proton_mass=1.007):
    """
    Convert neutral peptide masses into m/z values for a given charge state.

    This function takes a dictionary mapping peptide sequences to their
    neutral molecular masses and converts each mass into an m/z (mass-to-charge)
    value. The conversion uses the standard formula:

        m/z = (neutral_mass + charge * proton_mass) / charge

    Proton mass and charge state can be adjusted if needed.

    Parameters
    ----------
    peptide_mass_map : dict[str, float]
        A dictionary mapping peptide sequences to their neutral molecular masses.
    charge : int, optional
        The charge state used for m/z calculation. Default is +2.
    proton_mass : float, optional
        The mass of a proton added per charge. Default is 1.007 Da.

    Returns
    -------
    dict[str, float]
        A dictionary mapping each peptide sequence to its calculated m/z value.
    """

    mz_map = {}
    for pep, mass in peptide_mass_map.items():
        mz = (mass + charge * proton_mass) / charge
        mz_map[pep] = mz
    return mz_map

import numpy as np
import matplotlib.pyplot as plt

def plot_spectrum(mz_values, random_count_range=(0, 30000), seed=42):
    """
    Plot a simulated mass spectrum for a list of m/z values.

    This function creates a simple bar-plot spectrum by assigning a random
    intensity to each provided m/z value. The random intensities are drawn
    from a uniform integer range and can be reproduced by specifying a
    random seed.

    Parameters
    ----------
    mz_values : list[float]
        A list of m/z values to be plotted as peaks in the spectrum.
    random_count_range : tuple[int, int], optional
        Lower and upper bounds for the random intensity values assigned to
        each peak. Default is (0, 30000).
    seed : int, optional
        Random seed used for reproducible intensity generation. Default is 42.

    Returns
    -------
    None
        The function generates a matplotlib bar plot and does not return a value.
    """
    np.random.seed(seed)

    # Generate one random intensity for each m/z peak
    intensities = np.random.randint(
        random_count_range[0],
        random_count_range[1],
        size=len(mz_values)
    )

    # Plot the MS1 spectrum
    plt.figure(figsize=(12, 5))
    plt.bar(mz_values, intensities, width=0.5)

    plt.xlabel("m/z")
    plt.ylabel("Intensity")
    plt.title("Simulated MS1 Spectrum")

    plt.show()

def fragment_peptide(peptide):
    """
    Generate b-ion and y-ion fragment sequences for a peptide.

    This function simulates peptide backbone fragmentation as observed in
    MS/MS experiments. For a given peptide sequence, it returns all possible
    N-terminal (b-ions) and C-terminal (y-ions) fragment sequences.  
    Each fragment is represented as the corresponding peptide substring.

    Parameters
    ----------
    peptide : str
        Amino-acid sequence of the peptide to fragment.

    Returns
    -------
    list of str
        A list containing all b-ion fragments followed by all y-ion fragments.
        - b-ions: peptide[0:i]   for i in 1..len(peptide)-1
        - y-ions: peptide[-i:]   for i in 1..len(peptide)-1
    """
    fragments = []
    # b-ions
    for i in range(1, len(peptide)):
        fragments.append(peptide[:i])

    # y-ions
    for i in range(1, len(peptide)):
        fragments.append(peptide[-i:])

    return fragments

   