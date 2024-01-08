# Skin Doctor CP

Skin Doctor CP is a machine learning model for the classification of small organic
compounds into skin sensitizers and non-sensitizers. More specifically, the core of
Skin Doctor CP is a random forest binary classifier that is enveloped in an
aggregated Mondrian conformal prediction framework. This allows users to define an
error significance level (i.e. error rate) for classification. Predictions (i.e.
sensitizer or non-sensitizer) are thus only reported for compounds for which the
expected reliability reaches or exceeds the error rate defined by the user. The model 
was trained on a curated data set of 1285 compounds measured in the local lymph node 
assay (LLNA).

Skin Doctor CP has been developed in collaboration with [Beiersdorf 
AG](https://www.beiersdorf.com).

## Usage

Skin Doctor CP can be called from the **command line**. Examples:

```bash
# input in SMILES format
skindoctorcp "CCOC(=O)N1CCN(CC1)C2=C(C(=O)C2=O)N3CCN(CC3)C4=CC=C(C=C4)OC"

# significance_level1 to significance_level4 are numbers between 0 and 1
# (default: 0.05, 0.1, 0.2, 0.3)
skindoctorcp --significance_level1 0.05 --significance_level2 0.1 \
    --significance_level3 0.2 -- significance_level4 0.3 \
    "CCN(C)C(=O)OC1=CC=CC(=C1)C(C)N(C)C"

# input can be a file
skindoctorcp molecules.sdf > result.csv

# output format can be specified
skindoctorcp --output sdf molecules.smiles > result.sdf

# more information via --help
skindoctorcp --help
```

The model can be used in **Python**. Calling the ```predict``` function of the 
```SkinDoctorCPModel``` class results in a pandas DataFrame containing the prediction 
results for each input molecule.

```python
from skindoctor import SkinDoctorCPModel

model = SkinDoctorCPModel()

# "predict" method accepts a list of SMILES representations
df_predictions = model.predict(['CCN(C)C(=O)OC1=CC=CC(=C1)C(C)N(C)C'])

# ... or a list of file paths
df_predictions = model.predict(['part1.sdf', 'part2.sdf'])
```

The result DataFrame contains the columns:
* **mol_id**: unique number identifying the input molecule
* **input**: the raw representation provided as input (e.g. OCCCCC)
* **input_type**: the representation type of the input (e.g. smiles)
* **source**: the input source (e.g. my_molecules.sdf)
* **name**: the name of the input molecule (if provided in the input)
* **input_mol**: the RDKit molecule parsed from the input representation
* **preprocessed_mol**: the RDKit molecule after preprocessing
* **errors**: a list of errors that occured during reading or preprocessing the input
* **prediction_significance_level1, prediction_significance_level2, 
prediction_significance_level3, prediction_significance_level4**: Prediction of 
Skin Doctor CP with a reliability fulfilling the selected error significance
* **p_value_sensitizing**: p-value for the sensitizing class 
* **p_value_non_sensitizing**: p-value for the non-sensitizing class


## Contribute

```bash
conda env create -f environment.yml
conda activate skin-doctor
pip install -e .[dev,test]
ptw
```


## Contributors

* Anke Wilm
* Steffen Hirte
* Axinya Tokareva