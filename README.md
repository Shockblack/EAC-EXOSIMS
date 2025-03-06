# EAC-EXOSIMS

Implementation of the Exploratory Analytical Cases (EAC) for the Habitable Worlds Observatory into an EXOSIMS `json` format.
Follows observatory `yaml` files from the [HWO GOMAP Working Group](https://github.com/HWO-GOMAP-Working-Groups/Sci-Eng-Interface/tree/main) EAC reports that are packaged into a Python via `hwo-sci-eng`.
To get this script working, I recommend building the package associated with the `hwo-sci-eng` from source due to the version on `pip` having out-of-date observatory files that will throw errors (and have out-of-date values).
To build the repository from source, clone the `Sci-Eng-Interface` repository then run the following command while in the folder to install the package into your current environment:
```
pip install .
```
You must then create an environemnt variable for the reference data, following the instructions on `Sci-Eng-Interface`, which can be stored in your `.bashrc`/`.zshrc` file (this line of code is directly pulled from their README):
```
export SCI_ENG_DIR=/Users/tumlinson/anaconda3/envs/hwotools/lib/python3.12/site-packages/hwo_sci_eng
```