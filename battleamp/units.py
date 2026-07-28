"""MIC unit conversion using per-peptide molecular weight.

Converting between uM and ug/ml is peptide-specific: a fixed ug/ml value maps
to a different uM value for every sequence, in proportion to its molecular
weight. Models in the benchmark report in both units, so anything that puts
them in one table must harmonise first.

This duplicates the constants and formulas in workflow/scripts/evaluate.py.
The duplication is deliberate: evaluate.py depends on scikit-learn and scipy,
and the web service should not need either just to convert a number.
tests/test_units_parity.py asserts the two stay numerically identical.
"""

VALID_UNITS = ("ug/ml", "uM")

# Monoisotopic residue weights of the 20 standard amino acids (Da).
AA_MW = {
    "A":  71.03711, "R": 156.10111, "N": 114.04293, "D": 115.02694,
    "C": 103.00919, "E": 129.04259, "Q": 128.05858, "G":  57.02146,
    "H": 137.05891, "I": 113.08406, "L": 113.08406, "K": 128.09496,
    "M": 131.04049, "F": 147.06841, "P":  97.05276, "S":  87.03203,
    "T": 101.04768, "W": 186.07931, "Y": 163.06333, "V":  99.06841,
}

WATER_MW = 18.01056  # Da, terminal H + OH


def peptide_mw(sequence):
    """Molecular weight (Da): residue weights plus one water molecule."""
    mw = WATER_MW
    for aa in sequence.upper():
        mw += AA_MW.get(aa, 0.0)
    return mw


def convert(value, sequence, src_unit, dst_unit):
    """Convert a single MIC value between uM and ug/ml.

    Returns None if value is None, so callers can map over sparse columns.
    """
    if value is None:
        return None
    if src_unit == dst_unit:
        return value
    mw = peptide_mw(sequence)
    if src_unit == "uM" and dst_unit == "ug/ml":
        return value * mw / 1000.0
    if src_unit == "ug/ml" and dst_unit == "uM":
        return value * 1000.0 / mw
    raise ValueError(f"Unknown unit conversion: {src_unit} -> {dst_unit}")


def unit_suffix(unit):
    """Filesystem- and column-safe suffix for a unit ('ug/ml' -> 'ugml')."""
    return unit.replace("/", "").replace(" ", "")
