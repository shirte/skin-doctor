from pytest_bdd import given, parsers, when


@given("the SkinDoctorCP model", target_fixture="predictor")
def skin_doctor_model():
    from skin_doctor import SkinDoctorCPModel

    return SkinDoctorCPModel()


@given(
    parsers.parse("significance level 1 is {significance_level1:f}"),
    target_fixture="significance_level1",
)
def significance_level1(significance_level1):
    return significance_level1


@given(
    parsers.parse("significance level 2 is {significance_level2:f}"),
    target_fixture="significance_level2",
)
def significance_level2(significance_level2):
    return significance_level2


@given(
    parsers.parse("significance level 3 is {significance_level3:f}"),
    target_fixture="significance_level3",
)
def significance_level3(significance_level3):
    return significance_level3


@given(
    parsers.parse("significance level 4 is {significance_level4:f}"),
    target_fixture="significance_level4",
)
def significance_level4(significance_level4):
    return significance_level4


@when(
    parsers.parse("the model generates predictions for the molecule representations"),
    target_fixture="predictions",
)
def predictions(
    representations,
    predictor,
    input_type,
    significance_level1,
    significance_level2,
    significance_level3,
    significance_level4,
):
    return predictor.predict(
        representations,
        input_type=input_type,
        significance_level1=significance_level1,
        significance_level2=significance_level2,
        significance_level3=significance_level3,
        significance_level4=significance_level4,
        output_format="record_list",
    )
