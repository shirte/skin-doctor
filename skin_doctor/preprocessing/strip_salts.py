from typing import List, Tuple

from nerdd_module.preprocessing import Step
from rdkit.Chem import GetMolFrags, Mol, MolToSmiles


class StripSalts(Step):
    """
    Removes salts according to the rules defined in DOI: 10.1002/cmdc.201700673
    (Hit Dexter: A Machine-Learning Model for the Prediction of Frequent Hitters)
    """

    def __init__(self, remove_invalid_molecules=False):
        super().__init__()
        self.remove_invalid_molecules = remove_invalid_molecules

    def _run(self, mol: Mol) -> Tuple[Mol, List[str]]:
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
            and sorted_smiles_split_list[0][1] * 0.7 < sorted_smiles_split_list[1][1]
            and (
                not MolToSmiles(sorted_smiles_split_list[0][0], True)
                == MolToSmiles(sorted_smiles_split_list[1][0], True)
            )
        ):
            if self.remove_invalid_molecules:
                result_mol = None
            else:
                result_mol = mol
            errors.append("S1")
        else:
            if len(sorted_smiles_split_list) > 1:
                errors.append("S0")
            result_mol = sorted_smiles_split_list[0][0]

        return result_mol, errors
