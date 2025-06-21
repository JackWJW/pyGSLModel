# pyGSLModel

pyGSLModel is a python package designed to facilitate performing metabolic simulations of human GlycoSphingoLipid (GSL) metabolism. It also provides an easy to use strategy for performing transcriptomic integration with iMAT on TCGA and GTEX expression data (accessed from the Xena database).

## Quick Start Guide

### Model Preparation
pyGSLModel includes functions to download the Human-GEM () genome scale metabolic model of human metabolism and to subsequently convert gene names (with mygene API) and prune the model with fastcore (with pyfastcore) to create a smaller GSL specific model if desired. This has notable performance benefits.

```python
# Imports
from pyGSLModel import download_model, convert_genes, prune_model

# Downloading the HUMAN_GEM
model = download_model()

# Converting genes from ensembl IDs to gene names
model = convert_genes(model)

# Pruning the model with fastcore to  focus on GSL metabolism
model = prune_model(model)

# Removing unwanted transport reactions
model = remove_GSL_transport(model)
```

The pruning step takes longer as the function runs a number of test simulations on the full-scale model with different obejctives to ensure it doesn't miss any core reactions, it then subsequently performs fastcore pruning with the pyfastcore package. If you know you want to use the pruned GSL focused model, we made a dedicated function that downloads a pre-pruned and tidied model.

```python
# Imports
from pyGSLModel import download_GSL_model

# Downloading the pruned GSL model
model = download_GSL_model()
```

### Standard Simulations of Metabolism
With a prepared model, you are ready to run some metabolic simulations!
To do this you can make use of the `run_metabolic_model` function. This function allows you to simulate your model given a particular method and objective function.
The methods include standard FBA which would be enforced via `method = "FBA"`, a multiple-FBA approach in which multiple solutions are found for different lienar paths within GSL metabolism based on the selected objective and recombined into a distribution representative ensemble model: `method = "mFBA"`. To specify the GSL objective function for simulation you select one of `D14_Neuron` for a Day 14 I3Neuron, `D28_Neuron` for a Day 28 I3Neuron, `AC` for an IAstrocyte or `MG` for an IMicroglia in the `objective_choice=` argument. These define the set of reactions in the objective function that the solution will optimise for.

```python
# Imports
from pyGSLModel import run_metabolic_model

# Running a simulation
solution = run_metabolic_model(model, method="mFBA", objective_choice="D14_Neuron")
```

The `run_metabolic_model` function also allows for easy analysis of gene knockouts with the `knockout` argument. This is by default set to `"WT"` which won't change the model, however if you input a gene symbol, this will be knocked out of the model before the simulation is performed e.g., `knockout = "B4GALNT1"`.

```python
# Running a simulation with a gene knocked out.
solution_ko = run_metabolic_model(model, method="mFBA", objective_choice="D14_Neuron", knockout="B4GALNT1")
```

### Analysing and Plotting Results
As pyGSLModels focus is to simplify analysis of GSL metabolism, there are also included functions to analyse the metabolic simulation outputs.

```python
# Imports
from pyGSLModel import tabulate_model_results, plot_model_results

# First you can tabule the model results for key GSLs
results_df = tabulate_model_results(model, solution)

# You can then use that data frame as input to the default plotting function
results_fig = plot_model_results(results_df)
```

### Transcriptomic Integration
pyGSLModel includes a number of functions to integrate transcriptomic data with your metabolic model to generate fluxes.

There are three functions to integrate transcriptomic data with the iMAT method (implemented here with imatpy). This method doesn't require a conventional objective function, instead it uses the gene expression data to drive the simulation. This is useful for complex multicellular organisms and tissues, where a clear biological objective is difficult to discern.

Two of the functions utilise TCGA and GTEX data from the TCGA_TARGET_GTEX dataset (Xena).

The first of these functions: `TCGA_iMAT_Integrate()` takes TCGA data averaged over multiple samples for each cancer/normal tissue and performs iMAT based on genes of GSL metabolism. You can tune the function by adjusting the upper/lower quantiles for how the data is processed. An `upper_quantile = 0.25` and a `lower_quantile = 0.75` means the top 25% of expression values will be assigned a +1, the bottom 25% will be assigned -1 and the remaining values assigned a 0. This is the format needed for iMAT to optimise the model to prioritise high expression gene associated reactions.

The second of these functions: `TCGA_iMAT_sample_integrate()` performs simarly to `TCGA_iMAT_Integrate()`, however it doesn't average the data for a particular cancer/tissue, instead it calculates fluxes for each individual sample in the data set for the given `tissue` argument. This produces a lot more data and needs to run a lot of simulations so can take some time depending on the tissue you choose.

```python
# Imports
from pyGSLModel import TCGA_iMAT_integrate, TCGA_iMAT_sample_integrate

# Performing iMAT on TCGA and GTEX data where expression has been averaged across samples of the same disease/tissue
results_avg_df = TCGA_iMAT_integrate(model, upper_quantile = 0.25, lower_quantile = 0.75, epsilon=1, threshold=0.01)

# Performing iMAT on TCGA and GTEX data for each individual sample for a given tissue
results_sample_df = TCGA_iMAT_sample_integrate(model, tissue="Brain", upper_quantile = 0.25, lower_quantile=0.75, epsilon=1, threshold=0.01)
```

The third iMAT integration method allows the user to input their own data. The required format is a dataframe where each column represents a different sample and the index contains the gene symbols.

```python
# Imports
import pandas as pd
from pyGSLModel import iMAT_integrate

# Making a mock dataframe
d = {
    "Gene" : ["B4GALNT1", "ST3GAL5", "ST8SIA1"],
    "Sample_1" : [10,5,1],
    "Sample_2" : [1,5,10],
    "Sample_3" : [5,1,10],
}

mock_df = pd.DataFrame(d)
mock_df = mock_df.set_index("Gene")

# Perfoming iMAT
results_custom_df = iMAT_integrate(model, mock_df, upper_quantile = 0.25, lower_quantile = 0.75, epsilon=1, threshold=0.01)
```

## Dependencies and License
pyGSLModel is under a GNU license.
Dependecies include:
- requests (Apache)
- cobra (GNU)
- pyfastcore (MIT)
- mygene (BSD)
- pandas (BSD)
- matplotlib (PSF)
- seaborn (BSD)
- imatpy (MIT)
- numpy (BSD)
- networkx (BSD)
- pyvis (BSD)