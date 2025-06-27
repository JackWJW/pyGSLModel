# Importing Libraries
import requests, io

import cobra
from cobra.io import read_sbml_model
from cobra.manipulation import rename_genes
from cobra.core.solution import Solution

from pyfastcore import Fastcore

import mygene

### Function 1: Model download ###
def download_model():
    """
    Downloads the Human-GEM metabolic model from: https://raw.githubusercontent.com/SysBioChalmers/Human-GEM/main/model/Human-GEM.xml and reads in the model with cobrapy. Returns the model
    """
    # Downloading Genome Scale Model (HUMAN-GEM)
    print("Downloading  and Reading in Model")
    url = "https://raw.githubusercontent.com/SysBioChalmers/Human-GEM/main/model/Human-GEM.xml"
    input_model = requests.get(url)
    input_model.raise_for_status()

    #Reading in the model
    model = read_sbml_model(io.StringIO(input_model.text))
    print("Model succesfully downloaded and read in.")
    return model

### Function 2: Pruned Model Download ###
def download_GSL_model():
    """
    Downloads a pruned and tidied version of the HUMAN-GEM to preserve core GSL reactions and reads in the model with cobrapy. Returns the model
    """
    # Downloading Genome Scale Model (HUMAN-GEM)
    print("Downloading  and Reading in Model")
    url = "https://raw.githubusercontent.com/JackWJW/pyGSLModel/main/tidied_GSL_model/tidied_GSL_model.xml"
    input_model = requests.get(url)
    input_model.raise_for_status()

    #Reading in the model
    model = read_sbml_model(io.StringIO(input_model.text))
    print("Model succesfully downloaded and read in.")
    return model

### Function 2: Converting gene IDs to gene symbols ###
def convert_genes(model):
    """
    Takes the current model and converts Ensembl IDs to gene symbols (names) via the mygene api.

    Inputs:
    - model : a cobra.Model object
    """
    # Getting a list of ensembl ids from the model
    ensembl_ids = [g.id for g in model.genes]
    # Prepare mygene
    mg = mygene.MyGeneInfo()
    # Set up the query 
    res = mg.querymany(
        ensembl_ids,
        scopes='ensembl.gene',
        fields='symbol',
        species='human',
        as_dataframe=False
    )
    # Generate the mapping
    id2symbol = {
        entry['query']: entry.get('symbol',entry['query'])
        for entry in res
    }
    # Apply the mapping
    rename_genes(model,id2symbol)
    return model

### Function 3: Preparing the objective function dictionaries ###
def prepare_objective(model,objective_choice="D14_Neuron"):
    """
    Takes input model and defines Glycosphinglipid focused objective functions for Day 14 I3Neurons, Day 28 I3 Neurons, IMicroglia and IAstrocytes

    Inputs:
    - model : a cobra.Model object
    - objective_choice : "D14_Neuron", "D28_Neuron", "MG", "AC"
    """
    if objective_choice == "D14_Neuron":
        # Day 14 I3Neuron GSL composition
        obj_composition = {
            'MAR08167' : 0.08, #GM1b - A
            'MAR08186' : 0.19, #GD1b - D
            'MAR08189' : 0.29, #GT1b - D
            'MAR08181' : 0.08, #GQ1b - D
            'MAR08184' : 0.03, #GM3 - M
            'MAR08190' : 0.03, #GM2 - M
            'MAR08185' : 0.08, #GM1 - M
            'MAR08188' : 0.23, #GD1a - M
        }
        rxn_list = {}
        for rid, w in obj_composition.items():
            rxn = model.reactions.get_by_id(rid)
            rxn_list[rxn] = w

    elif objective_choice == "D28_Neuron":
        # Day 28 I3Neuron GSL composition
        obj_composition = {
            'MAR08167' : 0.17, #GM1b - A
            'MAR08181' : 0.06, #GQ1b - D
            'MAR08186' : 0.16, #GD1b - D
            'MAR08189' : 0.24, #GT1b - D
            'MAR08184' : 0.01, #GM3 - M
            'MAR08185' : 0.08, #GM1 - M
            'MAR08188' : 0.18, #GD1a - M
            'MAR08190' : 0.1, #GM2 - M
        }
        rxn_list = {}
        for rid, w in obj_composition.items():
            rxn = model.reactions.get_by_id(rid)
            rxn_list[rxn] = w

    elif objective_choice == "AC":
        # IAstrocyte GSL composition
        obj_composition = {
            'MAR08165' : 0.23, #GA2 - A
            'MAR08148' : 0.05, #Gb3 - B
            'MAR08149' : 0.03, #Gb4 - B
            'MAR08150' : 0.07, #Gb5/a2-3SpGB - B
            'MAR08179' : 0.03, #GD3 - D
            'MAR08277' : 0.03, #pGB - L
            'MAR08184' : 0.5, #GM3 - M
            'MAR08190' : 0.06, #GM2 - M
        }
        rxn_list = {}
        for rid, w in obj_composition.items():
            rxn = model.reactions.get_by_id(rid)
            rxn_list[rxn] = w
    
    elif objective_choice == "MG":
        # IMicroglia GSL composition
        obj_composition = {
            'MAR08165' : 0.02, #GA2 - A
            'MAR08168' : 0.07, #GD1alpha - A
            'MAR08148' : 0.11, #Gb3 - B
            'MAR08149' : 0.1, #Gb4 - B
            'MAR08150' : 0.07, #Gb5/a2-3SpGB - B
            'MAR08179' : 0.02, #GD3 - D
            'MAR08277' : 0.16, #pGB - L
            'MAR08190' : 0.02, #GM2 - M
            'MAR08184' : 0.43, #GM3 - M
        }
        rxn_list = {}
        for rid, w in obj_composition.items():
            rxn = model.reactions.get_by_id(rid)
            rxn_list[rxn] = w
    return rxn_list

### Function 4: Pruning model
def prune_model(model, objective_choice="D14_Neuron",remove_transport="Yes"):
    """
    Prunes the model using fastcore methods. Core-reactions will be collected via a series of trial simulations and metabolic subsytems involved in GSl metabolism. These core reactions will be used to prune the model.

    
    Inputs:
    - model : 
    - objective_choice : 
    - remove_transport : 
    """
    rxn_list = prepare_objective(model=model,objective_choice=objective_choice)
    model.objective = rxn_list
    
    if remove_transport == "Yes":
        target_subsystems = ["Sphingolipid metabolism", "Blood group biosynthesis"]
    elif remove_transport == "No":
        target_subsystems = ["Sphingolipid metabolism", "Blood group biosynthesis", "Transport reactions"]
    else:
        raise ValueError("Invalid remove_transport input, should be one of 'Yes' or 'No'")
    
    # Creating a list of reactions
    print("Collecting Core Reactions")
    subsystem_core_reactions_list = []
    for rxn in model.reactions:
        if rxn.subsystem in target_subsystems:
            subsystem_core_reactions_list.append(rxn.id)
                
    
    sim_core_reactions_list = []

    for obj_c in ["AC","MG","D14_Neuron","D28_Neuron"]:
        sol = run_metabolic_model(model,method="FBA",objective_choice=obj_c)
        for rid, flux in sol.fluxes.items():
            if flux > 0:
                sim_core_reactions_list.append(rid)
    
    core_reactions_multiples = subsystem_core_reactions_list + sim_core_reactions_list

    core_reaction_names = list(set(core_reactions_multiples))
    core_reaction_names.append("MAR00920")

    core_reactions = []
    for r in core_reaction_names:
        rxn = model.reactions.get_by_id(r)
        core_reactions.append(rxn)
    
    print("Core reactions collected, now pruning with FastCore")
    # Creating a fastcore method instance
    fc_builder = Fastcore(model=model, core_reactions=core_reactions)

    # Running fastcore
    fc_builder.fast_core()

    # Building model
    model = fc_builder.build_context_specific_model()
    print("Pruning with FastCore complete")
    return model

### Function 5: runnign multiple linear FBAs and then weights the solutions and combines into an ensemble network ###
def multi_fba(model, objective_choice):
    """
    Takes model and objective input and performs multiple linear FBA, one optimising for each lipid synthetic path before then weighting the flux outputs and combining into an ensemble network to best model the target reaction distribution.

    Inputs:
    - model : a cobra.Model object
    - objective_choice : One of "D14_Neuron", "D28_Neuron", "AC" or "MG". Defines the cell type lipid profile target.
    """

    # Defining Multiple Objectives for Each Pathway Given the objective choice
    if objective_choice == "D14_Neuron":
        ### Path A ###
        obj_composition_A = {
            'MAR08167' : 0.08, #GM1b - A
        }
        
        rxn_list_A = {}
        for rid, w in obj_composition_A.items():
            rxn = model.reactions.get_by_id(rid)
            rxn_list_A[rxn] = w
        
        ### Path D ###
        obj_composition_D = {
            'MAR08186' : 0.19, #GD1b - D
            'MAR08189' : 0.29, #GT1b - D
            'MAR08181' : 0.08, #GQ1b - D
        }
        
        rxn_list_D = {}
        for rid, w in obj_composition_D.items():
            rxn = model.reactions.get_by_id(rid)
            rxn_list_D[rxn] = w

        ### Path M ###
        obj_composition_M = {
            'MAR08184' : 0.03, #GM3 - M
            'MAR08190' : 0.03, #GM2 - M
            'MAR08185' : 0.08, #GM1 - M
            'MAR08188' : 0.23, #GD1a - M
        }
        
        rxn_list_M = {}
        for rid, w in obj_composition_M.items():
            rxn = model.reactions.get_by_id(rid)
            rxn_list_M[rxn] = w       

        ### Running Path A simulation ###
        model_A = model.copy()
        model_A.objective = rxn_list_A
        sol_A = model_A.optimize()

        ### Running Path D simulation ###
        model_D = model.copy()
        model_D.objective = rxn_list_D
        sol_D = model_D.optimize()   

        ### Running Path M simulation ###
        model_M = model.copy()
        model_M.objective = rxn_list_M
        sol_M = model_M.optimize()

        ### Computing Weights ###
        W_A = sum(obj_composition_A.values())
        W_D = sum(obj_composition_D.values())
        W_M = sum(obj_composition_M.values())

        total_W = W_A + W_D + W_M

        p_A, p_D, p_M = W_A/total_W, W_D/total_W, W_M/total_W

        ### Calculating the combined flux ###
        v_A, v_D, v_M = sol_A.fluxes, sol_D.fluxes, sol_M.fluxes

        combined_flux = (p_A*v_A) + (p_D*v_D) + (p_M*v_M)
        combined_flux = combined_flux.clip(upper=1000)

        ### Creating the solution object ###
        sol_combined = Solution(
            status = "optimal",
            objective_value=0.0,
            fluxes=combined_flux.to_dict()
        )

    if objective_choice == "D28_Neuron":
        ### Path A ###
        obj_composition_A = {
            'MAR08167' : 0.17, #GM1b - A
        }
        
        rxn_list_A = {}
        for rid, w in obj_composition_A.items():
            rxn = model.reactions.get_by_id(rid)
            rxn_list_A[rxn] = w
        
        ### Path D ###
        obj_composition_D = {
            'MAR08181' : 0.06, #GQ1b - D
            'MAR08186' : 0.16, #GD1b - D
            'MAR08189' : 0.24, #GT1b - D
        }
        
        rxn_list_D = {}
        for rid, w in obj_composition_D.items():
            rxn = model.reactions.get_by_id(rid)
            rxn_list_D[rxn] = w

        ### Path M ###
        obj_composition_M = {
            'MAR08184' : 0.01, #GM3 - M
            'MAR08185' : 0.08, #GM1 - M
            'MAR08188' : 0.18, #GD1a - M
            'MAR08190' : 0.1, #GM2 - M
        }
        
        rxn_list_M = {}
        for rid, w in obj_composition_M.items():
            rxn = model.reactions.get_by_id(rid)
            rxn_list_M[rxn] = w       

        ### Running Path A simulation ###
        model_A = model.copy()
        model_A.objective = rxn_list_A
        sol_A = model_A.optimize()

        ### Running Path D simulation ###
        model_D = model.copy()
        model_D.objective = rxn_list_D
        sol_D = model_D.optimize()   

        ### Running Path M simulation ###
        model_M = model.copy()
        model_M.objective = rxn_list_M
        sol_M = model_M.optimize()

        ### Computing Weights ###
        W_A = sum(obj_composition_A.values())
        W_D = sum(obj_composition_D.values())
        W_M = sum(obj_composition_M.values())

        total_W = W_A + W_D + W_M

        p_A, p_D, p_M = W_A/total_W, W_D/total_W, W_M/total_W

        ### Calculating the combined flux ###
        v_A, v_D, v_M = sol_A.fluxes, sol_D.fluxes, sol_M.fluxes

        combined_flux = (p_A*v_A) + (p_D*v_D) + (p_M*v_M)
        combined_flux = combined_flux.clip(upper=1000)

        ### Creating the solution object ###
        sol_combined = Solution(
            status = "optimal",
            objective_value=0.0,
            fluxes=combined_flux.to_dict()
        )

    if objective_choice == "AC":
        ### Path A ###
        obj_composition_A = {
            'MAR08165' : 0.23, #GA2 - A
        }
        
        rxn_list_A = {}
        for rid, w in obj_composition_A.items():
            rxn = model.reactions.get_by_id(rid)
            rxn_list_A[rxn] = w

        ### Path B ###
        obj_composition_B = {
            'MAR08148' : 0.05, #Gb3 - B
            'MAR08149' : 0.03, #Gb4 - B
            'MAR08150' : 0.07, #Gb5/a2-3SpGB - B
        }
        
        rxn_list_B = {}
        for rid, w in obj_composition_B.items():
            rxn = model.reactions.get_by_id(rid)
            rxn_list_B[rxn] = w
        
        ### Path D ###
        obj_composition_D = {
            'MAR08179' : 0.03, #GD3 - D
        }
        
        rxn_list_D = {}
        for rid, w in obj_composition_D.items():
            rxn = model.reactions.get_by_id(rid)
            rxn_list_D[rxn] = w

        ### Path L ###
        obj_composition_L = {
            'MAR08277' : 0.03, #pGB - L
        }
        
        rxn_list_L = {}
        for rid, w in obj_composition_L.items():
            rxn = model.reactions.get_by_id(rid)
            rxn_list_L[rxn] = w

        ### Path M ###
        obj_composition_M = {
            'MAR08184' : 0.5, #GM3 - M
            'MAR08190' : 0.06, #GM2 - M
        }
        
        rxn_list_M = {}
        for rid, w in obj_composition_M.items():
            rxn = model.reactions.get_by_id(rid)
            rxn_list_M[rxn] = w       

        ### Running Path A simulation ###
        model_A = model.copy()
        model_A.objective = rxn_list_A
        sol_A = model_A.optimize()

        ### Running Path B simulation ###
        model_B = model.copy()
        model_B.objective = rxn_list_B
        sol_B = model_B.optimize()

        ### Running Path D simulation ###
        model_D = model.copy()
        model_D.objective = rxn_list_D
        sol_D = model_D.optimize()   

        ### Running Path L simulation ###
        model_L = model.copy()
        model_L.objective = rxn_list_L
        sol_L = model_L.optimize()   

        ### Running Path M simulation ###
        model_M = model.copy()
        model_M.objective = rxn_list_M
        sol_M = model_M.optimize()

        ### Computing Weights ###
        W_A = sum(obj_composition_A.values())
        W_B = sum(obj_composition_B.values())
        W_D = sum(obj_composition_D.values())
        W_L = sum(obj_composition_L.values())
        W_M = sum(obj_composition_M.values())

        total_W = W_A + W_B + W_D + W_L + W_M

        p_A, p_B, p_D, p_L, p_M = W_A/total_W, W_B/total_W, W_D/total_W, W_L/total_W, W_M/total_W

        ### Calculating the combined flux ###
        v_A, v_B, v_D, v_L, v_M = sol_A.fluxes, sol_B.fluxes, sol_D.fluxes, sol_L.fluxes, sol_M.fluxes

        combined_flux = (p_A*v_A) + (p_B*v_B) + (p_D*v_D) + (p_L*v_L) + (p_M*v_M)
        combined_flux = combined_flux.clip(upper=1000)

        ### Creating the solution object ###
        sol_combined = Solution(
            status = "optimal",
            objective_value=0.0,
            fluxes=combined_flux.to_dict()
        )

    if objective_choice == "MG":
        ### Path A ###
        obj_composition_A = {
            'MAR08165' : 0.02, #GA2 - A
            'MAR08168' : 0.07, #GD1alpha - A
        }
        
        rxn_list_A = {}
        for rid, w in obj_composition_A.items():
            rxn = model.reactions.get_by_id(rid)
            rxn_list_A[rxn] = w

        ### Path B ###
        obj_composition_B = {
            'MAR08148' : 0.11, #Gb3 - B
            'MAR08149' : 0.1, #Gb4 - B
            'MAR08150' : 0.07, #Gb5/a2-3SpGB - B
        }
        
        rxn_list_B = {}
        for rid, w in obj_composition_B.items():
            rxn = model.reactions.get_by_id(rid)
            rxn_list_B[rxn] = w
        
        ### Path D ###
        obj_composition_D = {
            'MAR08179' : 0.02, #GD3 - D
        }
        
        rxn_list_D = {}
        for rid, w in obj_composition_D.items():
            rxn = model.reactions.get_by_id(rid)
            rxn_list_D[rxn] = w

        ### Path L ###
        obj_composition_L = {
            'MAR08277' : 0.16, #pGB - L
        }
        
        rxn_list_L = {}
        for rid, w in obj_composition_L.items():
            rxn = model.reactions.get_by_id(rid)
            rxn_list_L[rxn] = w

        ### Path M ###
        obj_composition_M = {
            'MAR08190' : 0.02, #GM2 - M
            'MAR08184' : 0.43, #GM3 - M
        }
        
        rxn_list_M = {}
        for rid, w in obj_composition_M.items():
            rxn = model.reactions.get_by_id(rid)
            rxn_list_M[rxn] = w       

        ### Running Path A simulation ###
        model_A = model.copy()
        model_A.objective = rxn_list_A
        sol_A = model_A.optimize()

        ### Running Path B simulation ###
        model_B = model.copy()
        model_B.objective = rxn_list_B
        sol_B = model_B.optimize()

        ### Running Path D simulation ###
        model_D = model.copy()
        model_D.objective = rxn_list_D
        sol_D = model_D.optimize()   

        ### Running Path L simulation ###
        model_L = model.copy()
        model_L.objective = rxn_list_L
        sol_L = model_L.optimize()   

        ### Running Path M simulation ###
        model_M = model.copy()
        model_M.objective = rxn_list_M
        sol_M = model_M.optimize()

        ### Computing Weights ###
        W_A = sum(obj_composition_A.values())
        W_B = sum(obj_composition_B.values())
        W_D = sum(obj_composition_D.values())
        W_L = sum(obj_composition_L.values())
        W_M = sum(obj_composition_M.values())

        total_W = W_A + W_B + W_D + W_L + W_M

        p_A, p_B, p_D, p_L, p_M = W_A/total_W, W_B/total_W, W_D/total_W, W_L/total_W, W_M/total_W

        ### Calculating the combined flux ###
        v_A, v_B, v_D, v_L, v_M = sol_A.fluxes, sol_B.fluxes, sol_D.fluxes, sol_L.fluxes, sol_M.fluxes

        combined_flux = (p_A*v_A) + (p_B*v_B) + (p_D*v_D) + (p_L*v_L) + (p_M*v_M)
        combined_flux = combined_flux.clip(upper=1000)

        ### Creating the solution object ###
        sol_combined = Solution(
            status = "optimal",
            objective_value=0.0,
            fluxes=combined_flux.to_dict()
        )

    return sol_combined

### Function 6: Removing unwanted reactions associated with GSL trafficking ###
def remove_GSL_transport(model):
    """
    Takes the model as input and removes a number of reactions associated with GSL trafficking reactions that can impact model solutions.

    Inputs:
      - model : a cobra.Model object
    """
    # Removing Unwated Reactions
    deg_rxns_remove = []
    deg_rxns_list = [
        "MAR12059", # Extracellular to Golgi transport of GM2
        "MAR12051", # Extracellular to Golgi trasport of GD3
        "MAR01351", # Synthesis of GA2 from galactosyl Glucosyl Ceramide by B4GALNT1
        "MAR08904", # GT1a transport from cytosol to golgi
        "MAR08882", # GD1c transport from cytosol to golgi
        "MAR08902", # GQ1b transport from cytosol to golgi
        "MAR08145", # Laccer transport from cytosol to golgi
        # "MAR08230", # 9-O-acetylated-GD3 transport from cytosol to golgi
        # "MAR08232", # 9-O-acetylated-GT3 transport from cytosol to golgi
    ]

    for drid in deg_rxns_list:
        deg_rxns_remove.append(model.reactions.get_by_id(drid))
    
    model.remove_reactions(deg_rxns_remove)
    return model

### Function 7: Performing a simulation ###
def run_metabolic_model(model,method="FBA",objective_choice="D14_Neuron",knockout="WT"):
    """
    Performs constraint-based metabolic simulation utilising the desired solution method.
    You can also select the objective function for the simulation and choose to knockout a gene.
    
    Inputs:
    - model : a cobra.Model object that should be the model of interest
    - method : should be one of "FBA" or "mFBA" for Linear FBA, or multiple linear FBA (mFBA) respctively.
    - objective_choice : defines which cell styled objective to use for the simulation and should be one of "D14_Neuron", "D28_Neuron", "AC" or "MG".
    - knockout : allows you to input a gene in the format "ABCDE" to be knocked out of the model before performing the simulation. This is set to "WT" for Wild-type by default.
    """
    # Checking method
    if method not in ["FBA", "mFBA"]:
        raise ValueError(f"Invalid method ({method}), please select one of 'FBA','mFBA'")
    
    # Checking objective_choice
    if objective_choice not in ["D14_Neuron", "D28_Neuron", "AC", "MG"]:
        raise ValueError(f"Invalid objective_choice ({objective_choice}), please select one of 'D14_Neuron','D28_Neuron','AC','MG'")

    # Running the simulation
    model = model.copy()

    rxn_list = prepare_objective(model,objective_choice)
    model.objective = rxn_list

    if knockout != "WT":
        try:
            model.genes.get_by_id(knockout).knock_out()
        except:
            print(f"Gene {knockout} not in model")

    if method == "FBA":
        sol = model.optimize()
    elif method == "mFBA":
        sol = multi_fba(model,objective_choice=objective_choice)
    else:
        raise ValueError("Selected solution method should be one of 'FBA' or 'mFBA'")
        
    return sol