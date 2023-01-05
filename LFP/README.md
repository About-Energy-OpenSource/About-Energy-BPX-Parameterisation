# About:Energy 2 Ah LFP|Graphite Cylindrical 18650 Cell Parameterisation for Battery Parameter eXchange (BPX)

About:Energy Limited, UK: 
[https://www.aboutenergy.io](https://www.aboutenergy.io)

This folder contains parameters and validation simulations for a Doyle-Fuller-Newman (DFN) model of a 2 Ah LFP|graphite cylindrical 18650 cell. The parameters are supplied by About:Energy to support the new open standard [Battery Parameter eXchange (BPX)](https://github.com/pybamm-team/BPX/), an outcome of the Faraday Institution Multi-scale Modelling project. 

About:Energy develops battery parameter sets for physics-based models by combining information from cycling data with electrochemical and physical measurements on electrodes harvested by cell teardown. For these cells, electrolyte, separator and thermal properties are informed from literature. Advanced datasheets and physics-based models for a wider range of commercially available cells will be made available as an About:Energy product in 2023.

The notebook [validation-plots](validation-plots.ipynb) shows comparisons between simulations of the parameterised DFN model and experimental data.\
[example-simulate-experiment](example-simulate-experiment.ipynb) demonstrates how to use the parameterised model to run your own simulations.\
[parameter_plot](parameter_plot.ipynb) plots the functional parameters, and converts appropriate parameters from PyBaMM to BPX format.

We note some discrepancies in voltage and capacity prediction at higher rates in the limit of low state-of-charge (SOC) - these are attributed to the use of a 1D+1D DFN model, which is not capable of accounting for inhomogeneity within the cylindrical cell. Additionally, by comparison to an NMC positive electrode (see NMC folder in this repository), the LFP positive electrode material is relatively less well described by the Fick's law diffusion and Butler-Volmer equation approximations defined by the basic BPX standard. A more refined model in conjunction with extensions to the BPX standard could widen the applicability of this parameter set.

The simulations use the package [PyBaMM (Python Battery Mathematical Modelling)](https://www.pybamm.org/), with the parameters are imported directly from the BPX JSON file using the method:
```bash
pybamm.ParameterValues.create_from_bpx("XXX")
```

## 🚀 Installation
In order to run the notebooks in this repository, you will need to install [PyBaMM](https://www.pybamm.org/) (version 22.12 or newer). We recommend installing within a [virtual environment](https://docs.python.org/3/tutorial/venv.html) in order to not alter any Python distribution files on your machine.

PyBaMM is available on GNU/Linux, MacOS and Windows. For more detailed instructions on how to install PyBaMM, see [the PyBaMM documentation](https://pybamm.readthedocs.io/en/latest/install/GNU-linux.html#user-install).

### Linux/Mac OS
To install the requirements on Linux/Mac OS use the following terminal commands, replacing "XXX" with the repository name:

1. Clone the repository
```bash
https://github.com/About-Energy-OpenSource/XXX.git
```
2. Change into the `XXX` directory 
```bash
cd XXX
```
3. Create a virtual environment
```bash
virtualenv env
```
4. Activate the virtual environment 
```bash
source env/bin/activate
```
5. Install the required packages
```bash 
pip install -r requirements.txt
```

### Windows
To install the requirements on Windows use the following terminal commands, replacing "XXX" with the repository name:

1. Clone the repository
```bash
https://github.com/About-Energy-OpenSource/XXX.git
```
2. Change into the `XXX` directory 
```bash
cd XXX
```
3. Create a virtual environment
```bash
virtualenv env
```
4. Activate the virtual environment 
```bash
\path\to\env\Scripts\activate
```
where `\path\to\env` is the path to the environment created in step 3 (e.g. `C:\Users\'Username'\env\Scripts\activate.bat`).

5. Install the required packages
```bash 
pip install -r requirements.txt
```

As an alternative, you can set up [Windows Subsystem for Linux](https://docs.microsoft.com/en-us/windows/wsl/about). This allows you to run a full Linux distribution within Windows.

### Troubleshooting
**Problem**: `ModuleNotFoundError: No module named 'wheel'`.

**Solution**: Try `pip install wheel` before `pip install -r requirements.txt`.

## 📫 Get in touch
If you have any questions or would like more information on battery parameterisation services, please get in touch via email: <contact@aboutenergy.co.uk>.

