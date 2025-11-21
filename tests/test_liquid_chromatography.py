import proteosim as ps

def test_predict_lc_retention_times():
    peptides = ["ACDE", "FGRST", "WAAA"]
    expected = {
        "ACDE": 0.9,
        "FGRST": 9.5,
        "WAAA": 16.1
    }

    actual = ps.predict_lc_retention_times(peptides)
    assert actual == expected
    

def test_select_retention_time_window():
    peptide_rt_map = {
        "PEP1": 5.0,
        "PEP2": 12.3,
        "PEP3": 25.8
    }

    # Fenster: 0 - 20 Minuten
    selected = ps.select_retention_time_window(
        peptide_rt_map,
        lower_ret_time=0,
        upper_ret_time=20
    )

    # Erwartet: Nur PEP1 und PEP2 liegen im Bereich
    expected = ["PEP1", "PEP2"]

    # Sortieren falls Reihenfolge unterschiedlich
    assert sorted(selected) == sorted(expected)