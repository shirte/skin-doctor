@skin_doctor
Feature: Consistent predictions

  Scenario Outline: Predictions stay consistent with previous versions
    Given an input molecule specified by '<input_smiles>'
    And the input type is 'smiles'
    And the SkinDoctorCP model
    And significance level 1 is <significance_level1>
    And significance level 2 is <significance_level2>
    And significance level 3 is <significance_level3>
    And significance level 4 is <significance_level4>
    
    When the model generates predictions for the molecule representations
    And the subset of the result where the input was not None is considered
    
    Then the value in column 'name' should be equal to '<name>'
    And the value in column 'p_sens_class' should be equal to <p_sens_class>
    And the value in column 'p_nonsens_class' should be equal to <p_nonsens_class>
    And the value in column 'cp_1' should be equal to <cp_1>
    And the value in column 'cp_2' should be equal to <cp_2>
    And the value in column 'cp_3' should be equal to <cp_3>
    And the value in column 'cp_4' should be equal to <cp_4>

    Examples:
    | significance_level1 | significance_level2 | significance_level3 | significance_level4 | name                     | input_smiles                                                                                                  | preprocessed_smiles                                                          | cp_1                        | cp_2                        | cp_3                        | cp_4                        | p_nonsens_class | p_sens_class | errors |
    | 0.05                | 0.1                 | 0.2                 | 0.3                 | Aciclovir                | C1=NC2=C(N1COCCO)N=C(NC2=O)N Aciclovir                                                                        | Nc1nc(=O)c2ncn(COCCO)c2[nH]1                                                 | ['Prediction not possible'] | ['Prediction not possible'] | ['Prediction not possible'] | ['Prediction not possible'] | 0.29            | 0.25         |        |
    | 0.05                | 0.1                 | 0.2                 | 0.3                 | Amiodarone               | CCN(CC)CCOc1c(I)cc(cc1I)C(=O)c2c3ccccc3oc2CCCC Amiodarone                                                     | CCCCc1oc2ccccc2c1C(=O)c1cc(I)c(OCCN(CC)CC)c(I)c1                             | ['Prediction not possible'] | ['Prediction not possible'] | ['Sensitizer']              | ['Sensitizer']              | 0.17            | 0.4          |        |
    | 0.05                | 0.1                 | 0.2                 | 0.3                 | Arsphenamine (Salvarsan) | C1=CC(=C(C=C1[As]=[As]C2=CC(=C(C=C2)O)N)N)O.Cl.Cl Arsphenamine (Salvarsan)                                    | None                                                                         | None                        | None                        | None                        | None                        | None            | None         | !1     |
    | 0.05                | 0.1                 | 0.2                 | 0.3                 | Cyclophosphamide         | C1CNP(=O)(OC1)N(CCCl)CCCl Cyclophosphamide                                                                    | O=P1(N(CCCl)CCCl)NCCCO1                                                      | ['Prediction not possible'] | ['Prediction not possible'] | ['Prediction not possible'] | ['Non-sensitizer']          | 0.31            | 0.23         |        |
    | 0.05                | 0.1                 | 0.2                 | 0.3                 | Doxorubicin              | C[C@H]1[C@H]([C@H](C[C@@H](O1)O[C@H]2C[C@@](Cc3c2c(c4c(c3O)C(=O)c5cccc(c5C4=O)OC)O)(C(=O)CO)O)N)O Doxorubicin | COc1cccc2c(O)c3c(O)c4c(c(O)c3c(O)c12)=C(OC1CC(N)C(O)C(C)O1)CC(O)(C(=O)CO)C=4 | ['Prediction not possible'] | ['Non-sensitizer']          | ['Non-sensitizer']          | ['Non-sensitizer']          | 0.52            | 0.09         |        |
    | 0.05                | 0.1                 | 0.2                 | 0.3                 | Hydrochlorothiazide      | O=S(=O)(N)c1c(Cl)cc2c(c1)S(=O)(=O)NCN2 Hydrochlorothiazide                                                    | NS(=O)(=O)c1cc2c(cc1Cl)NCNS2(=O)=O                                           | ['Prediction not possible'] | ['Prediction not possible'] | ['Non-sensitizer']          | ['Non-sensitizer']          | 0.47            | 0.11         |        |
    | 0.05                | 0.1                 | 0.2                 | 0.3                 | Levofloxacin             | C[C@H]1COc2c3n1cc(c(=O)c3cc(c2N4CCN(CC4)C)F)C(=O)O Levofloxacin                                               | CC1COc2c(N3CCN(C)CC3)c(F)cc3c(=O)c(C(=O)O)cn1c23                             | ['Prediction not possible'] | ['Prediction not possible'] | ['Non-sensitizer']          | ['Non-sensitizer']          | 0.42            | 0.15         |        |
    | 0.05                | 0.1                 | 0.2                 | 0.3                 | Metamizole (Sulpyrine)   | CC1=C(C(=O)N(N1C)C2=CC=CC=C2)N(C)CS(=O)(=O)[O-].O.[Na+] Metamizole (Sulpyrine)                                | Cc1c(N(C)CS(=O)(=O)O)c(=O)n(-c2ccccc2)n1C                                    | ['Prediction not possible'] | ['Non-sensitizer']          | ['Non-sensitizer']          | ['Non-sensitizer']          | 0.57            | 0.08         | S0     |
    | 0.05                | 0.1                 | 0.2                 | 0.3                 | Nifedipine               | CC1=C(C(C(=C(N1C=O)C)C(=O)OC)C2=CC=CC=C2[N+](=O)[O-])C(=O)OC Nifedipine                                       | COC(=O)C1=C(C)N(C=O)C(C)=C(C(=O)OC)C1c1ccccc1[N+](=O)[O-]                    | ['Prediction not possible'] | ['Prediction not possible'] | ['Prediction not possible'] | ['Non-sensitizer']          | 0.31            | 0.24         |        |
    | 0.05                | 0.1                 | 0.2                 | 0.3                 | Phenprocoumon            | OC=1c3ccccc3OC(=O)C=1C(CC)c2ccccc2 Phenprocoumon                                                              | CCC(c1ccccc1)c1c(O)c2ccccc2oc1=O                                             | ['Prediction not possible'] | ['Prediction not possible'] | ['Non-sensitizer']          | ['Non-sensitizer']          | 0.34            | 0.2          |        |
    | 0.05                | 0.1                 | 0.2                 | 0.3                 | Rivaroxaban              | O=C1COCCN1c2ccc(cc2)N3C[C@@H](OC3=O)CNC(=O)c4ccc(s4)Cl Rivaroxaban                                            | O=C(NCC1CN(c2ccc(N3CCOCC3=O)cc2)C(=O)O1)c1ccc(Cl)s1                          | ['Prediction not possible'] | ['Prediction not possible'] | ['Non-sensitizer']          | ['Non-sensitizer']          | 0.48            | 0.11         |        |

