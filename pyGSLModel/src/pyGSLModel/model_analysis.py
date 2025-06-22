# Imports
from .model_preparation import *
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from cobra import Model
from cobra.core.solution import Solution
import networkx as nx
from pyvis.network import Network

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
    plt.close(fig)
    return fig

def visualise_flux_network(model, solution, file_path="./flux_network_graph.html",height:str="600px",width:str="800px",met_col="#F6A6B2",rxn_col="#AED6F1"):
    """
    Generates a diagramatic visualisation of your model solution, with edges and nodes weighted by flux.
    Writes an html object for the visualisation

    Inputs:
    - model : cobrapy model object
    - solution : model solution object
    - file_path : file path and file name for where to save the output e.g., './flux_network_graph.html'
    - height : height as string e.g., "100px"
    - width : width as string e.g., "100px" or "100%"
    - met_col : defines colour of metaboltie nodes
    - rxn_col : defines colour of reaction/gene nodes
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
    keep_rxns = [
        rxn for rxn in model.reactions
        if any(m.id in keep_metabolites for m in rxn.metabolites)
    ]

    # 2) Build bipartite graph: metabolites ↔ reactions
    G = nx.DiGraph()
    # Add metabolite nodes
    for mid in keep_metabolites:
        if mid in model.metabolites:
            met = model.metabolites.get_by_id(mid)
            G.add_node(mid,
                       label=met.name,
                       type="met",
                       color=met_col,
                       size=20)

    # Add reaction nodes and edges
    max_flux = max(abs(solution.fluxes.get(r.id, 0.0)) for r in keep_rxns) or 1.0
    for rxn in keep_rxns:
        flux = float(solution.fluxes.get(rxn.id, 0.0))
        abs_flux = abs(flux)
        # size reaction by flux
        G.add_node(
            rxn.id,
            label=",".join(g.id for g in rxn.genes) or rxn.id,
            type="rxn",
            color=rxn_col,
            size= (abs_flux / max_flux) * 40 + 10,
        )
        # connect reactants → reaction → products
        for met, coeff in rxn.metabolites.items():
            if met.id not in keep_metabolites:
                continue
            if coeff < 0:
                src, tgt = met.id, rxn.id
            else:
                src, tgt = rxn.id, met.id
            G.add_edge(
                src,
                tgt,
                value=(abs_flux / max_flux) * 40 + 1,
                title=f"flux={flux:.3g}"
            )

    # 3) Build pyvis Network
    net = Network(
        directed=True,
        height=height,
        width=width,
        cdn_resources="remote"
    )
    # Load from networkx
    net.from_nx(G)

    # 4) Tweak physics/layout
    net.force_atlas_2based(gravity=-50, central_gravity=0.01, spring_length=200)
    net.write_html(file_path)
    print(f"html file wrtten to: {file_path}")