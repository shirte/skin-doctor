from nerdd_module.preprocessing import FilterByElement, FilterByWeight, Pipeline

from .canonicalize_tautomer import CanonicalizeTautomer
from .do_smiles_roundtrip import DoSmilesRoundtrip
from .neutralize_charges import NeutralizeCharges
from .strip_salts import StripSalts


class SkinDoctorPreprocessingPipeline(Pipeline):
    def __init__(self):
        super().__init__(
            steps=[
                DoSmilesRoundtrip(remove_stereo=True),
                StripSalts(),
                NeutralizeCharges(),
                FilterByWeight(
                    min_weight=0,
                    max_weight=2500,
                    remove_invalid_molecules=True,
                ),
                FilterByElement(
                    allowed_elements=[
                        "H",
                        "B",
                        "C",
                        "N",
                        "O",
                        "F",
                        "Si",
                        "P",
                        "S",
                        "Cl",
                        "Se",
                        "Br",
                        "I",
                    ],
                    remove_invalid_molecules=True,
                ),
                CanonicalizeTautomer(
                    remove_stereo=True,
                    remove_invalid_molecules=True,
                ),
                DoSmilesRoundtrip(remove_stereo=False),
            ]
        )
