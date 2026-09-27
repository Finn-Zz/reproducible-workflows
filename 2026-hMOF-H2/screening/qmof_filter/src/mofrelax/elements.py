"""Small element helpers kept independent of any structure library."""

from __future__ import annotations

# Metalloids are intentionally excluded. The set follows a chemistry-oriented
# definition suitable for deciding whether a dihedral rotates around an
# organic bond; it is not intended as an oxidation-state model.
METALS = frozenset(
    {
        "Li", "Be", "Na", "Mg", "Al", "K", "Ca", "Sc", "Ti", "V", "Cr",
        "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Rb", "Sr", "Y", "Zr",
        "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Cs",
        "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy",
        "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir",
        "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Fr", "Ra", "Ac", "Th", "Pa",
        "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No",
        "Lr",
    }
)


def is_metal(symbol: str) -> bool:
    return symbol in METALS


def element_class(symbol: str) -> str:
    if is_metal(symbol):
        return "M"
    if symbol == "H":
        return "H"
    return "X"

