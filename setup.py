from setuptools import find_packages, setup

setup(
    name="skin-doctor",
    version="0.2.2",
    maintainer="Johannes Kirchmair",
    maintainer_email="johannes.kirchmair@univie.ac.at",
    packages=find_packages(),
    url="https://github.com/molinfo-vienna/skin-doctor",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    license="BSD 3-Clause License",
    include_package_data=True,
    install_requires=[
        "rdkit==2020.09.1",
        "scikit_learn==0.23.2",
        "pandas~=1.2.1",
        "numpy==1.19.2",
        "nerdd-module>=0.3.11",
        # avoid warnings about numpy.distutils
        "setuptools < 60.0",
        # install importlib-resources and importlib-metadata for old Python versions
        "importlib-resources>=5; python_version<'3.9'",
        "importlib-metadata>=4.6; python_version<'3.10'",
    ],
    extras_require={
        "dev": ["mypy", "ruff"],
        "test": [
            "pytest",
            "pytest-watch",
            "pytest-cov",
            "pytest-bdd",
            "hypothesis",
            "hypothesis-rdkit",
        ],
    },
    entry_points={
        "console_scripts": [
            "skindoctorcp=skin_doctor.__main__:main_skin_doctor_cp",
            "skindoctorternary=skin_doctor.__main__:main_skin_doctor_ternary",
        ],
    },
)
