# XRB Analysis

## Overview

This analysis framework is designed to load lightcurves generated and saved in the `XRB_Data_Prep` directory. It utilizes these lightcurves to perform two distinct analyses, which are executed through the `XRB_Sens.py` script. Configuration settings are managed through `config.py`.

## Analyses

### 1. Individual Source Searches

The first analysis focuses on individual source searches. This process involves several key functions:

- **`do_lc_trials`**: Runs lightcurve trials for individual sources.
- **`collect_lc_trials`**: Gathers results from the lightcurve trials.
- **`find_lc_nsig`**: Calculated the sens/discovery potential of the lightcurve results.
- **`plot_lc_bias`**: Plots the bias in the lightcurve analysis.

### 2. Stacking Analysis

The second analysis employs a stacking approach, utilizing all available lightcurves. The relevant functions for this analysis include:

- **`do_stacking_trials`**: Runs trials on the stacked lightcurves.
- **`collect_stacking_trials`**: Gathers results from the stacking trials.
- **`find_stacking_nsig`**: Calculates the sens/discovery potential of the stacking results.
- **`plot_stacking_bias`**: Plots the bias in the stacking analysis.

## Resource Management

Given that the analyses are resource-intensive, two additional functions are provided to facilitate the submission of trials to NPX:

- **`submit_do_lc_trials`**: Submits individual source lightcurve trials to NPX for processing.
- **`submit_do_stacking_trials`**: Submits stacking trials to NPX for efficient resource management.

## Requirements

To successfully run the analysis, the following Python packages are required:

- **csky**: Version to be finalized once branch is merged.
- **numpy**: For numerical operations and data handling.
- **click**: For creating command-line interfaces.
- **pandas**: For data manipulation and analysis.
- **Submitter**: For generating DAGs to manage trial submissions. Available in ssclafani949/Submitter

## Usage


## Setup
The first step will be to install all requirements listed above.

Then you will need to modify `config.py` to point to your directories
**condor** is a cluster located at UMD.  So we will focus on the NPX implementation.

ana_dir and base_dir should be created in your data/user/ based on your username.  But feel free to modify if you want another temporary location
The rest of the keys are local to the XRB repo or are just read only.

submitter_config_npx will need to be modified to point to a virtual environment with csky and and the few other packages.  Here we are using 
cvmfs py3-v4.2.1
but it can be done with a more recent version if necessary

## Running Trials

Trials can no be run on cobalt by running the command:
eg:
```python XRB_Sens.py do-lc-trials --n-trials 100 --n-sig 0 --src_gamma 2.0 --thresh 0.0 --lag 0.0 --name 'Cyg_X_dash_1' --seed 0 --inject_gp```

```python XRB_Sens.py do-lc-trials --n-trials 100 --n-sig 10 --src_gamma 2.0 --thresh 0.0 --lag 0.0 --name 'Cyg_X_dash_1' --seed 0 --noinject_gp```

will run 100 background trials (injecting the GP) and 100 signal trials at ns = 10 (not injecting any GP)

the trials take a lot of time to run.  I have presaved the trials for all the parameters that can be used to calculate sensitvity, discovery potential and plot biases.

These are available here: 


`/data/user/ssclafani/data/analyses/XRB_baseline_v1.0/combo/lc/stacking/TSD_chi2.dict`

`/data/user/ssclafani/data/analyses/XRB_baseline_v1.0/combo/lc/TSD_chi2.dict`

