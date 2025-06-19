# Imports
from .model_preparation import *
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

### Function 8: Generating a results table ###
def tabulate_model_results(model, sol):
    """
    Takes the model object and solution object from a simulation and produces a pandas dataframe containing the results with respect to GSL related reactions.

    Inputs:
    - model : A cobra.Model object
    - sol : A solution object output from run_metabolic model
    """
    # Creating a list of metabolites that we want to Record
    keep_metabolites = [
    "MAM01904g", # GA1
    "MAM01905g", # GA2
    "MAM01941g", # GD1a
    "MAM01942g", # GD1alpha
    "MAM01943g", # GD1b
    "MAM01945g", # GD1c
    "MAM01946g", # GD2
    "MAM01947g", # GD3
    "MAM01912g", # GB5
    "MAM01959g", # GB4
    "MAM01960g", # GB3
    "MAM02008g", # GM1
    "MAM02010g", # GM1b
    "MAM02011g", # GM2
    "MAM02015g", # GM3
    "MAM02023g", # GQ1b
    "MAM02025g", # GQ1c
    "MAM02028g", # GT1a
    "MAM02030g", # GT1b
    "MAM02031g", # GT1c
    "MAM02032g", # GT2
    "MAM02033g", # GT3
    "MAM02328g", # LacCer_Pool
    "MAM02346g", # LC3
    "MAM02347g", # LC4
    "MAM02330g", # paragloboside
    "MAM02904g" # sialylparagloboside
    ]

    # Collect every reaction involving any of the metabolites in our list
    keep_reactions   = [
    rxn.id
    for rxn in model.reactions
    if any(m.id in keep_metabolites for m in rxn.metabolites)
    ]

    #Building a dataframe to store the results
    rows = []
    for rxn_id in keep_reactions:
        rxn = model.reactions.get_by_id(rxn_id)
        reactants = [f"{met.id} ({met.name})"
                        for met, c in rxn.metabolites.items() if c < 0]
        products  = [f"{met.id} ({met.name})"
                        for met, c in rxn.metabolites.items() if c > 0]
        product_keep = [met.name
                        for met, c in rxn.metabolites.items()
                        if c > 0 and met.id in keep_metabolites]
        genes = [g.id for g in rxn.genes]
        rows.append({
            "Reaction ID":   rxn.id,
            "Reactants":     ", ".join(reactants),
            "Products":      ", ".join(products),
            "Key Product":  ", ".join(product_keep),
            "Genes":         ", ".join(genes)
        })
    temp_df = pd.DataFrame(rows)
    temp_df = temp_df[temp_df["Key Product"] != ""].copy()
    temp_df.set_index("Reaction ID", inplace=True)

    #Assigning the flux values
    flux_map = dict(sol.fluxes)
    temp_df["Flux (mmol/gDW/hr)"] = temp_df.index.map(flux_map)
    temp_df["Relative GSL Flux (%)"] = temp_df["Flux (mmol/gDW/hr)"] / temp_df["Flux (mmol/gDW/hr)"].sum() * 100
    results_data = temp_df.reset_index().sort_values("Flux (mmol/gDW/hr)",ascending=False)
    return results_data

### Function 9: Plotting simulation results ###
def plot_model_results(data):
    """
    Takes the results data from tabulate_model_results and plots two barcharts showing Flux against Key Product and Flux against Genes.

    Inputs:
    - data : the tabulted pandas data frame from tabulate_model_results
    """
    sns.set_style("ticks")
    df_vis = data[["Key Product", "Genes", "Flux (mmol/gDW/hr)","Relative GSL Flux (%)"]]

    fig, axs = plt.subplots(ncols=2,figsize=(12,6))
    sns.barplot(data=data,x="Key Product", y="Relative GSL Flux (%)",errorbar=None,ax=axs[0],color="skyblue",edgecolor="black")
    sns.barplot(data=data,x="Genes",y="Relative GSL Flux (%)",errorbar=None,ax=axs[1],color="skyblue",edgecolor="black")

    for ax in fig.axes:
        ax.tick_params(rotation = 90)
        ax.set(xlabel=None)

    return fig