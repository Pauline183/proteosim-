import matplotlib.pyplot as plt
import pyteomics.achrom as achrom
def predict_lc_retention_times(peptides):
    """
    Add a short description here.

    Parameters
    ----------

    Returns
    -------
    """

     # Leeres Dictionary zur Speicherung: Peptid → Retentionszeit
    peptide_retention_times = {}

    # Schleife über jedes Peptid in der Eingabeliste
    for pep in peptides:
        # Berechne die Retentionszeit für ein Peptid mit dem Achrom-Modell (Es nimmt die Peptidsequenzund summiert die hydrophoben Beiträge aller Aminosäuren anhand der Tabelle RCs.)
        rt = achrom.calculate_RT(pep, achrom.RCs_guo_ph7_0)

        # Speichere die auf 2 Dezimalstellen gerundete RT im Dictionary
        peptide_retention_times[pep] = round(rt, 2)

    # Gib das vollständige Mapping zurück
    return peptide_retention_times

def plot_retention_time(retention_times, resolution=30):
    """
    Plot a histogram of peptide retention times.

    This function visualizes the distribution of retention times either from a list
    of numeric retention-time values or from a dictionary that maps peptides to
    their predicted retention times. If a dictionary is provided, only the numeric
    retention-time values are extracted for plotting.

    Parameters
    ----------
    retention_times : list[float] or dict[str, float]
        Either a list of retention-time values (in minutes) or a dictionary mapping
        peptide sequences to their retention-time predictions.
    resolution : int, optional
        Number of histogram bins. Higher values create a more fine-grained distribution.
        Default is 30.

    Returns
    -------
    None
        The function creates and displays a matplotlib histogram but returns no value.
    """
    # Wenn der Input ein Dictionary ist (z. B. {"PEPTIDE1": 12.3, ...}),dann wollen wir nur die **Retention-Time-Werte** extrahieren, weil nur die Zahlen ins Histogramm gehören
    # weil nur die Zahlen ins Histogramm gehören
    
    if isinstance(retention_times, dict):
        rt_values = list(retention_times.values())
    else:
        rt_values = retention_times

    

    plt.hist(rt_values, bins=resolution, alpha=0.7)
    plt.xlabel('Retention Time')
    plt.ylabel('Frequency')
    plt.title('Retention Time Distribution')
    plt.show()

def select_retention_time_window(peptide_rt_map, lower_ret_time, upper_ret_time):
    """
    Filter peptides by retention-time window.

    This function selects all peptides whose predicted retention times fall
    within a specified time interval. The input must be a dictionary mapping
    peptide sequences to their retention-time values. Only entries whose
    retention time is between the lower and upper bound (inclusive) are kept.

    Parameters
    ----------
    peptide_rt_map : dict[str, float]
        A dictionary mapping peptide sequences to their predicted retention times.
    lower_ret_time : float
        The lower bound of the retention-time window (inclusive).
    upper_ret_time : float
        The upper bound of the retention-time window (inclusive).

    Returns
    -------
    dict[str, float]
        A filtered dictionary containing only peptides with retention times
        inside the specified window.
    """
    # Dictionary comprehension:
    # Wir durchlaufen alle (peptide, retention_time)-Paare …
    # und behalten nur diejenigen, deren RT zwischen lower und upper liegt.
    selected = {
        pep: rt
        for pep, rt in peptide_rt_map.items()
        if lower_ret_time <= rt <= upper_ret_time
    }

    return selected