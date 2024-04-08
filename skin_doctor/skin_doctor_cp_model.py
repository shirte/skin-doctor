import gzip
import os
import time
from copy import deepcopy
from multiprocessing import Pool
from typing import List

import numpy as np
import pandas as pd
from joblib import load
from nerdd_module import AbstractModel
from rdkit.Chem import MACCSkeys, Mol
from rdkit.DataStructs.cDataStructs import BitVectToText, CreateFromBitString

from .preprocessing import SkinDoctorPreprocessingPipeline

try:
    # works in python 3.9+
    from importlib.resources import files
except ImportError:
    from importlib_resources import files

import warnings

# avoid warnings in old versions of numpy
warnings.filterwarnings("ignore", message="numpy.dtype size changed")
warnings.filterwarnings("ignore", message="numpy.ufunc size changed")

__all__ = ["SkinDoctorCPModel"]

current_dir = os.path.dirname(os.path.abspath(__file__))

cores = 1

# load ml models
# first note: loading the raw pkl models takes around 10 seconds
# second note: when gzipping the models, the loading time is around 30 seconds
# third note: when multiple models are put into one file (and gzipped), the loading
#             time is around 20 seconds (compared to uncompressed)
mlm = [
    load(
        files("skin_doctor")
        .joinpath("models/skin_doctor_cp/cp_models")
        .joinpath(f"clf_{number}.pkl")
        .open("rb")
    )
    for number in range(100)
]

# load alpha values
# note: this looks slow, but it's actually quite fast
mlm_alphas = [
    [
        pd.read_csv(
            files("skin_doctor")
            .joinpath(f"models/skin_doctor_cp/cp_alpha_values")
            .joinpath(f"alphasCali_{number}_0.csv")
            .open("r")
        )
        .squeeze("columns")
        .tolist(),
        pd.read_csv(
            files("skin_doctor")
            .joinpath("models/skin_doctor_cp/cp_alpha_values")
            .joinpath(f"alphasCali_{number}_1.csv")
            .open("r")
        )
        .squeeze("columns")
        .tolist(),
    ]
    for number in range(100)
]


def get_maccs_fp(mol):
    """
    Method to calculate maccsFps

    :param mol: rdkit molecule
    :return: rdkits MACCS key
    """
    return CreateFromBitString(BitVectToText(MACCSkeys.GenMACCSKeys(mol))[1:])


def get_cp_results(maccs, CPmodelList, CPalphaList, signis, classLabels, procs):
    """
    Method that produces the results of the conformal prediction model

    :param maccs: Bitstring, MACCS fingerprint for the compound of interest
    :param CPmodelList: List of the 100 classifiers for mondrian CP
    :param CPalphaList: List of alpha values for the 100 models applied to the 100 calibration sets
    :param errorSignis: List of the four user defined error significance levels
    :param classLabels: List of the Strings for the class labels
    :return: List with the results as strings for the four error significance levels
    """
    with Pool(procs) as pool:
        gen = ((CPmodelList[x], CPalphaList[x], maccs, 2) for x in range(100))
        pValuesAllMap = pool.map_async(inner_loop_mondrian, gen)
        pValuesAll = pValuesAllMap.get()
        pValuesMedian = get_median_pValues(pValuesAll)
        resultsBoolean = [
            predict_activity_boolean(pValuesMedian, errorSigni) for errorSigni in signis
        ]
        resultsString = [
            predict_activity_string(resultBoolean, n_classes=2, classLabels=classLabels)
            for resultBoolean in resultsBoolean
        ]
    resultsString.append(np.round(pValuesMedian[0][0], 2))
    resultsString.append(np.round(pValuesMedian[0][1], 2))
    return resultsString


def inner_loop_mondrian(args):
    """
    Method that does the conformal prediction magic

    :param args: Tuple, (CPmodel, alpha values of the calibration set on this model, fingerprint of the molecule, number of classes)
    :return: p-values as results of the cp model

    """
    clf, alphasCali, X_test, n_classes = args
    # calculate the alpha values for test set
    proba_test = clf.predict_proba(X_test)
    alphasTest = []
    for counter, proba in enumerate(proba_test):
        alphas = calculate_alpha_value_multiclass(proba)
        alphasTest.append(
            alphas
        )  # alphas_test is a matrix with a line for each test compound and an alphas list in each line
    # now see at which position in calibration set the test compounds would be
    # -> calculate p values
    pValues = calculate_p_values_multiclass(alphasCali, alphasTest)
    return pd.DataFrame(pValues)


def calculate_alpha_value_multiclass(probas):
    """
    Method that calualtes the alpha values of the calibration or test set

    param probas: list with floats, that are the class probabilities returned by the underlying classifier
    return: list of alpha values
    """
    alphaList = []
    for index, proba in enumerate(probas):
        probas_copy = list(deepcopy(probas))
        del probas_copy[index]
        alpha = 0.5 - ((proba - max(probas_copy)) / 2)
        alphaList.append(alpha)
    return alphaList


def calculate_p_values_multiclass(alphas_cali, alphas_test):
    """
    Method that calculates the p-values for the test set

    param alphas_cali: list, alpha values of the calibration set
    param alphas_test: list, alpha values of the test set

    return: p-values of the test set
    """
    p_values = []
    for index, alphas in enumerate(alphas_test):
        ps = []
        for class_index, alpha in enumerate(alphas):
            # etwas umstaendlich geschrieben, weil searchsorted nur arrays in ascending order annimt
            p = (
                len(alphas_cali[class_index])
                - np.searchsorted(alphas_cali[class_index], alpha, side="left")
                + 1
            ) / (len(alphas_cali[class_index]) + 1)
            ps.append(p)
        p_values.append(ps)
    return p_values


def predict_activity_string(predictedActivities, n_classes, classLabels):
    """
    Method, that creates a class label for each compound

    param predictedActivities: list, predicted activities as boolean
    param n_classes: Number of classes the model differntiates
    param classLabels: List of strings, that are assigned to the different classes

    return: list with class labels the compounds are assigned to as strings
    """
    resultsString = []
    for pred in predictedActivities:
        if sum(pred) != 1:
            resultsString.append("Prediction not possible")
        else:
            for c in range(n_classes):
                if pred[c] == True:
                    resultsString.append(classLabels[c])
    return resultsString


def get_median_pValues(pList):
    """
    Method that calulates the median p-values for each class from a list of p-values

    param pList: list of p-values for each class
    return: median p-values
    """
    pValuesMedian = []
    for compound in range(len(pList[0])):
        perClass = []
        for state in range(len(list(pList[0]))):
            p_all = []
            for p in pList:
                p_all.append(p.iloc[compound][state])
            p_median = np.median(p_all)
            perClass.append(p_median)
        pValuesMedian.append(perClass)
    return pValuesMedian


def predict_activity_boolean(pValues, errorSignificance):
    """
    Method that decides on the predicted class membership of a compound based on the selected error significance level

    param pValues: list, p-values of the test set
    param errorSignficance: float, selected error significance level
    return: list, of predicted activities
    """
    predictedActivities = []
    for entry in pValues:
        pred_boolean = []
        for classentry in entry:
            if classentry > errorSignificance:
                pred_boolean.append(True)
            else:
                pred_boolean.append(False)
        predictedActivities.append(pred_boolean)
    return predictedActivities


def predict(
    mols,
    significance_level_1: float = 0.05,
    significance_level_2: float = 0.10,
    significance_level_3: float = 0.20,
    significance_level_4: float = 0.30,
):
    significance_levels = [
        significance_level_1,
        significance_level_2,
        significance_level_3,
        significance_level_4,
    ]

    # calculate features
    class_labels = ["Non-sensitizer", "Sensitizer"]

    # prepare descriptors
    maccs = [get_maccs_fp(mol) for mol in mols]

    # skin doctor CP
    cp_results = [
        get_cp_results(
            [i],
            mlm,
            mlm_alphas,
            significance_levels,
            class_labels,
            cores,
        )
        for i in maccs
    ]

    results = pd.DataFrame(
        dict(
            p_nonsens_class=[j[4] for j in cp_results],
            p_sens_class=[j[5] for j in cp_results],
        )
    )

    for i in range(4):
        results[f"cp_{i+1}"] = [j[i] for j in cp_results]
        results[f"significance_err_{i+1}"] = significance_levels[i]

    return results


class SkinDoctorCPModel(AbstractModel):
    def __init__(self):
        super().__init__(preprocessing_pipeline=SkinDoctorPreprocessingPipeline())

    def _predict_mols(
        self,
        mols: List[Mol],
        significance_level1: float = 0.05,
        significance_level2: float = 0.10,
        significance_level3: float = 0.20,
        significance_level4: float = 0.30,
    ) -> pd.DataFrame:
        for sl in [
            significance_level1,
            significance_level2,
            significance_level3,
            significance_level4,
        ]:
            assert 0 <= sl <= 1, "Significance levels must be between 0 and 1"

        assert (
            significance_level1
            < significance_level2
            < significance_level3
            < significance_level4
        ), "Significance levels must be in increasing order"

        return predict(
            mols,
            significance_level1,
            significance_level2,
            significance_level3,
            significance_level4,
        )
