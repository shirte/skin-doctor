import pandas as pd
from pytest_bdd import given, parsers, scenario, then


@scenario(
    "features/consistent_predictions.feature",
    "Predictions stay consistent with previous versions",
)
def test_consistent_predictions():
    pass


@given(
    parsers.parse("an input molecule specified by '{input}'"),
    target_fixture="representations",
)
def representations(input):
    return [input]


@then(
    parsers.parse("the value in column '{column_name}' should contain only '{value}'")
)
def check_column_membership(predictions, column_name, value):
    if value == "(none)":
        assert all(
            pd.isnull(predictions[column_name])
        ), f"Column {column_name} must be none"
    else:
        assert all(
            value in values for values in predictions[column_name]
        ), f"Column {column_name} contains value {value}"
