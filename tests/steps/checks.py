import json

import pandas as pd
from pytest_bdd import parsers, then


@then(parsers.parse("The result should contain the columns:\n{column_names}"))
def check_result_columns(predictions, column_names):
    column_names = column_names.strip()
    for c in column_names.split("\n"):
        assert (
            c in predictions.columns
        ), f"Column {c} not in predictions {predictions.columns.tolist()}"


@then(parsers.parse("the value in column '{column_name}' should be between 0 and 1"))
def check_neighbor_column(subset, column_name):
    assert (0 <= subset[column_name]).all()
    assert (subset[column_name] <= 1).all()


@then(parsers.parse("the value in column '{column_name}' should be '{expected_value}'"))
def check_column_value(predictions, column_name, expected_value):
    value = predictions[column_name].iloc[0]

    # expected value is always provided as string
    # try to convert to float if possible
    try:
        expected_value = float(expected_value)
    except ValueError:
        pass

    if expected_value == "(none)":
        # if expected_value is the magic string "(none)", we expect None
        assert pd.isnull(value), f"Column {column_name} is assigned to {value} != None"
    else:
        # otherwise, we expect the value to be equal to the expected value
        assert (
            value == expected_value
        ), f"Column {column_name} is assigned to {value} != {expected_value}"


@then(
    parsers.parse(
        "the value in column '{column_name}' should be a subset of {superset}"
    )
)
def check_column_subset(predictions, column_name, superset):
    superset = set(json.loads(superset))

    assert all(
        set(value).issubset(superset)
        for value in predictions[column_name]
        if pd.notnull(value)
    ), f"Column {column_name} contains value not in {superset}"
