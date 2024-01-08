from nerdd_module import auto_cli

from .skin_doctor_cp_model import SkinDoctorCPModel


@auto_cli
def main():
    return SkinDoctorCPModel()
