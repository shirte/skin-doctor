from nerdd_module import auto_cli

from .skin_doctor_cp_model import SkinDoctorCPModel
from .skin_doctor_ternary_model import SkinDoctorTernaryModel


@auto_cli
def main_skin_doctor_ternary():
    return SkinDoctorTernaryModel()


@auto_cli
def main_skin_doctor_cp():
    return SkinDoctorCPModel()
