@skin_doctor
Feature: Valid predictions

  Scenario Outline: Predictions are valid
    Given a random seed set to <seed>
    And the input type is '<input_type>'
    And a list of <num_molecules> random molecules, where <num_none> entries are None
    And the representations of the molecules
    And the SkinDoctorCP model
    And significance level 1 is <significance_level1>
    And significance level 2 is <significance_level2>
    And significance level 3 is <significance_level3>
    And significance level 4 is <significance_level4>

    When the model generates predictions for the molecule representations
    And the subset of the result where the input was not None is considered

    Then The result should contain the columns:
            cp_1
            cp_2
            cp_3
            cp_4
            p_sens_class
            p_nonsens_class
    And the value in column 'cp_1' should be a subset of ["Non-sensitizer", "Sensitizer", "Prediction not possible"]
    And the value in column 'cp_2' should be a subset of ["Non-sensitizer", "Sensitizer", "Prediction not possible"]
    And the value in column 'cp_3' should be a subset of ["Non-sensitizer", "Sensitizer", "Prediction not possible"]
    And the value in column 'cp_4' should be a subset of ["Non-sensitizer", "Sensitizer", "Prediction not possible"]
    And the value in column 'p_sens_class' should be between 0 and 1
    And the value in column 'p_nonsens_class' should be between 0 and 1


  Examples:
  | seed | significance_level1 | significance_level2 | significance_level3 | significance_level4 | num_molecules | num_none | input_type |
  | 1    | 0.05                | 0.1                 | 0.2                 | 0.3                 | 10            | 0        | smiles     |
  | 2    | 0.05                | 0.1                 | 0.2                 | 0.3                 | 10            | 1        | smiles     |
  | 3    | 0.05                | 0.1                 | 0.2                 | 0.3                 | 10            | 2        | smiles     |
  | 4    | 0.05                | 0.1                 | 0.2                 | 0.3                 | 10            | 10       | smiles     |
  | 1    | 0.05                | 0.1                 | 0.2                 | 0.3                 | 10            | 0        | smiles     |
  | 2    | 0.05                | 0.1                 | 0.2                 | 0.3                 | 10            | 1        | smiles     |
  | 3    | 0.05                | 0.1                 | 0.2                 | 0.3                 | 10            | 2        | smiles     |
  | 4    | 0.05                | 0.1                 | 0.2                 | 0.3                 | 10            | 10       | smiles     |
