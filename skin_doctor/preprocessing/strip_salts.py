from typing import List, Tuple

from nerdd_module import Problem
from nerdd_module.preprocessing import PreprocessingStep
from rdkit.Chem import GetMolFrags, Mol, MolToSmiles


class StripSalts(PreprocessingStep):
    """
    Removes salts according to the rules defined in DOI: 10.1002/cmdc.201700673
    (Hit Dexter: A Machine-Learning Model for the Prediction of Frequent Hitters)
    """

    def __init__(self, remove_invalid_molecules=False):
        super().__init__()
        self.remove_invalid_molecules = remove_invalid_molecules

    def _preprocess(self, mol: Mol) -> Tuple[Mol, List[Problem]]:
        errors = []

        mol_frags = GetMolFrags(mol, asMols=True)
        smiles_split_list = []
        for frag in mol_frags:
            smiles_split_list.append((frag, frag.GetNumAtoms()))
        sorted_smiles_split_list = sorted(
            smiles_split_list, key=lambda x: x[1], reverse=True
        )
        if (
            len(sorted_smiles_split_list) > 1
            # second largest fragment is at least 70% of the largest fragment:
            and sorted_smiles_split_list[0][1] * 0.7 < sorted_smiles_split_list[1][1]
            and not (
                MolToSmiles(sorted_smiles_split_list[0][0], True)
                == MolToSmiles(sorted_smiles_split_list[1][0], True)
            )
        ):
            if self.remove_invalid_molecules:
                result_mol = None
            else:
                result_mol = mol
            errors.append(
                Problem(
                    "largest_fragment_not_unique",
                    (
                        "Could not identify the largest fragment, because there are multiple large "
                        "fragments of similar size."
                    ),
                )
            )
        else:
            if len(sorted_smiles_split_list) > 1:
                errors.append(
                    Problem(
                        "largest_fragment_not_unique",
                        "Could not identify the largest fragment, because there are multiple large "
                        "fragments of similar size.",
                    )
                )
            result_mol = sorted_smiles_split_list[0][0]

        return result_mol, errors
