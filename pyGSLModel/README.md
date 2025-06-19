# pyGSLModel

pyGSLModel is a python package designed to facilitate performing metabolic simulations of human GlycoSphingoLipid (GSL) metabolism. It also provides an easy to use strategy for performing transcriptomic integration with iMAT on TCGA and GTEX expression data (accessed from the Xena database).

## Quick Start Guide

pyGSLModel includes functions to download the Human-GEM () genome scale metabolic model of human metabolism and to subsequently convert gene names (with mygene API) and prune the model with fastcore (with pyfastcore) to create a smaller GSL specific model if desired. This has notable performance benefits.

```python
# Importing functions
from pyGSLModel import download_model, convert_genes, prune_model

# Downloading the HUMAN_GEM
model = download_model()

# Converting genes from ensembl IDs to gene names
model = convert_genes(model)

# Pruning the model with fastcore to  focus on GSL metabolism
model = prune_model(model)
```

The pruning step takes longer as the function runs a number of test simulations on the full-scale model with different obejctives to ensure it doesn't miss any core reactions, it then subsequently performs fastcore pruning with the pyfastcore package. If you know you want to use the pruned GSL focused model, we made a dedicated function that downloads a pre-pruned and tidied model.

```python
# Importing function
from pyGSLModel import download_GSL_model

# Downloading the pruned GSL model
model = download_GSL_model()
```

Once your model is prepared, the simplest thing is to run a simulation.
To do this you can make use of the `run_metabolic_model` function. This function allows you to simulate your model given a particular method and objective function.
The methods include standard FBA which would be enforced via `method = "LFBA"`, a multiple-FBA approach in which multiple solutions are found for different lienar paths within GSL metabolism based on the selected objective and recombined into a distribution representative ensemble model: `method = "mFBA"`. To specify the GSL objective function for simulation you select one of `D14_Neuron` for a Day 14 I3Neuron, `D28_Neuron` for a Day 28 I3Neuron, `AC` for an IAstrocyte or `MG` for an IMicroglia in the `objective_choice=` argument. These define the set of reactions in the objective function that the solution will optimise for.

```python
# Importing function
from pyGSLModel import run_metabolic_model

# Running a simulation
solution = run_metabolic_model(model, method="mFBA", objective_choice="D14_Neuron")
```

The `run_metabolic_model` function also allows for easy analysis of gene knockouts with the `knockout` argument. This is by default set to `"WT"` which won't change the model, however if you input a gene symbol, this will be knocked out of the model before the simulation is performed e.g., `knockout = "B4GALNT1"`.

```python
# Running a simulation with a gene knocked out.
solution_ko = run_metabolic_model(model, method="mFBA", objective_choice="D14_Neuron", knockout="B4GALNT1")
```

As pyGSLModels focus is to simplify analysis of GSL metabolism, there are also included functions to analyse the metabolic simulation outputs.

```python
# Importing functions
from pyGSLModel import tabulate_model_results, plot_model_results

# First you can tabule the model results for key GSLs
results_df = tabulate_model_results(model, solution)

# You can then use that data frame as input to the default plotting function
results_fig = plot_model_results(results_df)
```
