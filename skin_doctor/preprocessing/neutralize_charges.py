from typing import List, Tuple

from nerdd_module import Problem
from nerdd_module.preprocessing import PreprocessingStep
from rdkit.Chem import AllChem, Mol, MolFromSmarts, MolFromSmiles

__all__ = ["neutralization_reactions", "NeutralizeCharges"]

patterns = [
    # Imidazoles
    ("[n+;H]", "n"),
    # Amines
    ("[N+;!H0]", "N"),
    # Carboxylic acids and alcohols
    ("[$([O-]);!$([O-][#7])]", "O"),
    # Thiols
    ("[S-;X1]", "S"),
    # Sulfonamides
    ("[$([N-;X2]S(=O)=O)]", "N"),
    # Enamines
    ("[$([N-;X2][C,N]=C)]", "N"),
    # Tetrazoles
    ("[n-]", "[nH]"),
    # Sulfoxides
    ("[$([S-]=O)]", "S"),
    # Amides
    ("[$([N-]C=O)]", "N"),
]

neutralization_reactions = [
    (MolFromSmarts(x), MolFromSmiles(y, False)) for x, y in patterns
]


class NeutralizeCharges(PreprocessingStep):
    """
    Neutralizes the molecules according to the patterns defined in
    ruleSets.neutralizationPatterns.
    """

    def __init__(self):
        super().__init__()

    def _preprocess(self, mol: Mol) -> Tuple[Mol, List[Problem]]:
        errors = []
        for reactant, product in neutralization_reactions:
            while mol and mol.HasSubstructMatch(reactant):
                rms = AllChem.ReplaceSubstructs(mol, reactant, product)
                mol = rms[0]
                try:
                    mol.UpdatePropertyCache()
                    if not mol:
                        errors.append(
                            Problem(
                                "neutralization_failed", "Neutralizing charges failed."
                            )
                        )
                except ValueError:
                    mol = None

        return mol, errors
