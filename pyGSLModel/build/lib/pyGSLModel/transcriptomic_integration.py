from .model_preparation import *
from .model_analysis import *

from imatpy.parse_gpr import gene_to_rxn_weights
from imatpy.imat import imat

import numpy as np

### Function 10 Integrating transcriptomic data for Glioblastoma/Glioma into the model ###
def TCGA_DEFlux_integrate(model, objective_choice="AC",tissue="Brain",comparison="TCGA",fc_threshold=0.5):
    """
    This function integrates TCGA transcriptomic data for Cancer data accessed and tidied from the Xena database into the metabolic model.
    This is acheived by scaling the reaction bounds  with a differential E-flux style approach based on the relative expression in a cancer sample when compared to normal tissue samples for a predefined set of GSl-synthetic reactions.
    The function runs the transcripomic integrated metabolic simulations and returns a pandas DataFrame with one row per sample (indeced by sample ID), wich one column per reaction, with the flux encoded as values.

    Inputs:
    - model : a cobrapy model object.
    - objective_choice : defines which cell styled objective to use for the simulation and should be one of "D14_Neuron", "D28_Neuron", "AC" or "MG".
    - tissue : defines which tissue to select, this will be used to select the cancer and normal tissues. For example, is "Brain" is selected, Glioblastoma and Glioma data will be accessed alongside TCGA normal tissue data.
               Options include: 'White blood cell', 'Adrenal gland', 'Bladder', 'Brain', 'Breast', 'Cervix', 'Bile duct', 
                                'Colon', 'Lymphatic tissue', 'Esophagus', 'Head and Neck region', 'Kidney', 'Liver', 'Lung', 
                                'Lining of body cavities', 'Ovary', 'Pancreas', 'Paraganglia', 'Prostate', 'Rectum', 'Skin', 
                                Stomach', 'Testis', 'Thymus', 'Thyroid Gland', 'Uterus', 'Endometrium', 'Eye'.
                Please make sure your objective_choice is reasonable given the tissue of choice.
    - comparison : This defines which 'normal' data should be used to calculate relative expression. "TCGA" is the preferred options as this minimises challenges associated to batch effects.
                 Option 1: "TCGA" which uses Solid Tissue Normal samples taken adjacent to cancer tissue.
                 Option 2: "GTEX" which uses GTEX normal tissue data.
    - fc_threshold : Differential expression foldchange threshold where below which flux will be weighted based on FC (for reactions which are active in normal samples) or reaction upper bound set to 0 (for reactions which are not normally active) (default 0.5).
    """
    # Checking objective_choice
    if objective_choice not in ["D14_Neuron", "D28_Neuron", "AC", "MG"]:
        raise ValueError(f"Invalid objective_choice ({objective_choice}), please select one of 'D14_Neuron','D28_Neuron','AC','MG'")
    
    # Downloading Transcriptomic Data
    df_input = pd.read_csv("https://raw.githubusercontent.com/JackWJW/pyGSLModel/main/Xena-TCGA_TARGET_GTEX_Data/Xena_Data.tsv",sep='\t')

    # Checking tissue input is in the dataframe
    all_sites = df_input["_primary_site"].unique().tolist()
    if tissue not in all_sites:
        raise ValueError(f"Error: tissue '{tissue}' not found in dataset."
                         f"Please use one of: {all_sites}")

    # Creating datasets for the relative expression calculations based onuse input
    if comparison == "TCGA":
        df_input = df_input[(df_input["_primary_site"]==tissue)&(df_input["_study"]=="TCGA")].copy()
        df_cancer = df_input[df_input["_sample_type"]=="Primary Tumor"].copy()
        df_normal = df_input[df_input["_sample_type"]=="Solid Tissue Normal"].copy()
        if df_normal.shape[0] < 3:
            raise ValueError("Error: Insufficient 'Solid Tissue Normal' samples in TCGA data, try a different tissue or GTEX comparison")
        
    elif comparison == "GTEX":
        df_input = df_input[(df_input["_primary_site"]==tissue)&(df_input["_study"].isin(["TCGA","GTEX"]))].copy()
        df_cancer = df_input[df_input["_sample_type"]=="Primary Tumor"].copy()
        df_normal = df_input[df_input["_study"]=="GTEX"].copy()
        if df_normal.shape[0] < 3:
            raise ValueError("Error: Insufficient 'GTEX' samples in dataset, try a different tissue or TCGA comparison")

    else:
        print("invalid option for comparison selected, exiting function")
        raise ValueError("Error: Invalid comparison input, should be one of 'TCGA' or 'GTEX'")

    # Calculating and storing the average expression for each gene in the normal samples

    GENE_LIST = ['A4GALT', 'B3GALNT1', 'B3GALT4', 'B3GALT5', 'B3GNT3', 'B3GNT5', 
                'B4GALNT1', 'B4GALT5', 'B4GALT6', 'GAL3ST1', 'ST3GAL1', 'ST3GAL2', 
                'ST3GAL4', 'ST3GAL5', 'ST3GAL6', 'ST6GALNAC2', 'ST6GALNAC3', 'ST6GALNAC4', 
                'ST6GALNAC5', 'ST6GALNAC6', 'ST8SIA1', 'ST8SIA5', 'UGCG', 'UGT8']
    
    df_cancer[GENE_LIST] = 2**df_cancer[GENE_LIST] - 1
    df_normal[GENE_LIST] = 2**df_normal[GENE_LIST] - 1
    
    for gene_col in GENE_LIST:
        normal_exp = df_normal[gene_col].mean()
        ew = df_cancer[gene_col] / normal_exp
        df_cancer[f"{gene_col}_EW"] = np.clip(ew, 0.1, 10)

    # Running a base simulation to get starting flux values
    sol_base = run_metabolic_model(model,method="mFBA",objective_choice=objective_choice)
    print("Base model Solved")

    # Preparing results from base simulation
    base_df = tabulate_model_results(model, sol_base)
    rxn_ids = base_df["Reaction ID"].tolist()
    base_flux_series = base_df.set_index("Reaction ID")["Flux (mmol/gDW/hr)"]

    #Performing simulation for each sample
    sim_counter = 0
    sim_total = df_cancer.shape[0]

    all_rows = {}
    for _, row in df_cancer.iterrows():
        sim_counter += 1
        print(f"Simulations Performed:{sim_counter}/{sim_total}")

        sample_id = row["sample"]
        model_copy = model.copy()

        # Calculating Flux Bounds
        for rid in rxn_ids:
            rxn = model_copy.reactions.get_by_id(rid)
            base_flux = float(base_flux_series.loc[rid])

            gene_ews = []
            for g in rxn.genes:
                col_name = f"{g.id}_EW"
                if col_name in row:
                    gene_ews.append(row[col_name])
            
            if len(gene_ews) > 0:
                EW = max(gene_ews)
            else:
                EW = 1.0
            
            if base_flux > 0:
                if EW < fc_threshold:
                    rxn.upper_bound = base_flux * EW
            else:
                if EW < fc_threshold:
                    rxn.upper_bound = 0
        
        # Running the simulation 
        sol_sample = run_metabolic_model(model_copy,method="mFBA",objective_choice=objective_choice)

        # Storing flux results for the reactions of interest with each reaction as a column and the flux as the value in that column for that sample.
        sample_df = tabulate_model_results(model_copy, sol_sample)
        sample_df = sample_df[["Key Product", "Relative GSL Flux (%)"]].copy().set_index("Key Product")
        sample_df = sample_df.T.copy()
        sample_df["sample"] = sample_id
        all_rows[f"df_{sim_counter}"] = sample_df

    # Building the dataframe
    flux_by_sample_df = pd.concat(all_rows.values(), axis=0,ignore_index=True).set_index("sample")

    return flux_by_sample_df

def TCGA_iMAT_integrate(model, upper_quantile = 0.25, lower_quantile = 0.75, epsilon=1, threshold=0.01):
    """
    Performs iMAT transcriptomic integration (using the imatpy package) to generate fluxes for different cancers from the TCGA and normal tissue from GTEX (accessed via Xena).
    Uses the average gene expression for different TCGA cancers and GTEX tissue data, using only gene expresssion for GSL metabolism related genes.

    Inputs:
    - model : a cobrapy model object
    - upper_quantile : Defines the upper bound percentage of gene expression values for a sample to be assigned 1 for iMAT. If 0.25, the top 25% would be assigned 1
    - lower_quantile : Defines the lower bound percentage of gene expression values for a sample to be assigned -1 for iMAT. If 0.75, the bottom 25% would be assigned -1
    - epsilon : iMAT maximises the sum of high expressing reactions with flux > epsilon (default 1)
    - threshold : Alongside epsilon, iMAT maximises the sum of low expressing reactions with flux < threshold (default 0.01)
    """
    # Downloading and selecting cancer transcriptomic data
    df_input = pd.read_csv("https://raw.githubusercontent.com/JackWJW/pyGSLModel/main/Xena-TCGA_TARGET_GTEX_Data/Xena_Data.tsv",sep='\t')
    df_cancer = df_input[(df_input["_sample_type"].isin(["Primary Tumor", "Normal Tissue"]))&(df_input["_study"].isin(["TCGA","GTEX"]))].copy()

    # Defining list of genes for analysis
    GENE_LIST = ['A4GALT', 'B3GALNT1', 'B3GALT4', 'B3GALT5', 'B3GNT3', 'B3GNT5', 
                'B4GALNT1', 'B4GALT5', 'B4GALT6', 'GAL3ST1', 'ST3GAL1', 'ST3GAL2', 
                'ST3GAL4', 'ST3GAL5', 'ST3GAL6', 'ST6GALNAC2', 'ST6GALNAC3', 'ST6GALNAC4', 
                'ST6GALNAC5', 'ST6GALNAC6', 'ST8SIA1', 'ST8SIA5', 'UGCG', 'UGT8']

    # Preparing averaged data for each cancer
    df_cancer = df_cancer[GENE_LIST+["primary disease or tissue"]].copy()
    df_avgs = df_cancer.groupby(["primary disease or tissue"]).mean().T.copy()

    # Defining a helper function to convert expression into 1, 0 or -1 for high, neutral and low expressed genes
    def convert_col(col):
        u_q = col.quantile(upper_quantile)
        l_q = col.quantile(lower_quantile)

        converted = col.copy()
        converted[col > l_q] = 1
        converted[col < u_q] = -1
        converted[(col >= u_q) & (col <= l_q)] = 0

        return converted

    # Converting Values
    df_converted = df_avgs.apply(convert_col).copy()

    # Perform iMAT simulations for each cancer
    all_rows = {}
    all_genes = [g.id for g in model.genes]
    colnames = df_converted.columns.to_list()
    imat_counter = 0
    imat_total = len(colnames)
    for col in colnames:
        imat_counter += 1
        print(f"Simulations Performed:{imat_counter}/{imat_total}")
        model_copy = model.copy()
        model_weights = pd.Series(df_converted[col])
        model_weights = model_weights.reindex(all_genes, fill_value=0)
        imat_weights = gene_to_rxn_weights(model=model_copy,gene_weights=model_weights)
        imat_results = imat(model=model_copy,rxn_weights=imat_weights,epsilon=epsilon,threshold=threshold)

        #Tabulating results
        sample_df = tabulate_model_results(model_copy, imat_results)
        sample_df = sample_df[["Key Product", "Flux (mmol/gDW/hr)"]].copy().set_index("Key Product")
        sample_df = sample_df.T.copy()
        sample_df["Cancer"] = col
        all_rows[f"{col}_iMAT"] = sample_df

    # Building the dataframe
    imat_data = pd.concat(all_rows.values(), axis=0,ignore_index=True).set_index("Cancer")

    return imat_data

def TCGA_iMAT_sample_integrate(model, tissue, upper_quantile = 0.25, lower_quantile=0.75, epsilon=1, threshold=0.01):
    """
    Performs iMAT transcriptomic integration (using the imatpy package) to generate fluxes for different cancer and normal tissue samples from TCGA and GTEX (accessed via Xena).
    Generates a flux preduction based on gene expression via iMAT for every sample corresponding to the specified tissue of interest.
    Returns a dataframe containing metadata information for each sample such as tissue, disease, DFI, etc as well as flux data through each key lipid.

    Inputs:
    - model : a cobrapy model object
    - tissue : defines which tissue to select, this will be used to select the cancer and normal tissues.
               Options include: 'White blood cell', 'Adrenal gland', 'Bladder', 'Brain', 'Breast', 'Cervix', 'Bile duct', 
                                'Colon', 'Lymphatic tissue', 'Esophagus', 'Head and Neck region', 'Kidney', 'Liver', 'Lung', 
                                'Lining of body cavities', 'Ovary', 'Pancreas', 'Paraganglia', 'Prostate', 'Rectum', 'Skin', 
                                Stomach', 'Testis', 'Thymus', 'Thyroid Gland', 'Uterus', 'Endometrium', 'Eye'.
    - upper_quantile : Defines the upper bound percentage of gene expression values for a sample to be assigned 1 for iMAT. If 0.25, the top 25% would be assigned 1
    - lower_quantile : Defines the lower bound percentage of gene expression values for a sample to be assigned -1 for iMAT. If 0.75, the bottom 25% would be assigned -1
    - epsilon : iMAT maximises the sum of high expressing reactions with flux > epsilon (default 1)
    - threshold : Alongside epsilon, iMAT maximises the sum of low expressing reactions with flux < threshold (default 0.01)
    """
    # Downloading transcriptomic data
    df_input = pd.read_csv("https://raw.githubusercontent.com/JackWJW/pyGSLModel/main/Xena-TCGA_TARGET_GTEX_Data/Xena_Data.tsv",sep='\t')

    # Checking tissue input is in the dataframe
    all_sites = df_input["_primary_site"].unique().tolist()
    if tissue not in all_sites:
        raise ValueError(f"Error: tissue '{tissue}' not found in dataset."
                         f"Please use one of: {all_sites}")

    

    # Storing list of key genes:
    GENE_LIST = ['A4GALT', 'B3GALNT1', 'B3GALT4', 'B3GALT5', 'B3GNT3', 'B3GNT5', 
            'B4GALNT1', 'B4GALT5', 'B4GALT6', 'GAL3ST1', 'ST3GAL1', 'ST3GAL2', 
            'ST3GAL4', 'ST3GAL5', 'ST3GAL6', 'ST6GALNAC2', 'ST6GALNAC3', 'ST6GALNAC4', 
            'ST6GALNAC5', 'ST6GALNAC6', 'ST8SIA1', 'ST8SIA5', 'UGCG', 'UGT8']
    
    # Data preparation
    df_input = df_input[(df_input["_primary_site"]==tissue)&(df_input["_study"].isin(["TCGA", "GTEX"]))].copy()
    df_mini = df_input[["sample","_sample_type","primary disease or tissue"]].copy()
    df_format = df_input[GENE_LIST+["sample"]].set_index("sample").T.copy()

    # Defining a helper function to convert expression into 1, 0 or -1 for high, neutral and low expressed genes
    def convert_col(col):
        u_q = col.quantile(upper_quantile)
        l_q = col.quantile(lower_quantile)

        converted = col.copy()
        converted[col > l_q] = 1
        converted[col < u_q] = -1
        converted[(col >= u_q) & (col <= l_q)] = 0

        return converted

    # Converting Values
    df_converted = df_format.apply(convert_col).copy()

    # Perform iMAT simulations for each cancer
    all_rows = {}
    all_genes = [g.id for g in model.genes]
    colnames = df_converted.columns.to_list()
    imat_counter = 0
    imat_total = len(colnames)
    for col in colnames:
        imat_counter += 1
        print(f"Simulations Performed:{imat_counter}/{imat_total}")
        model_copy = model.copy()
        model_weights = pd.Series(df_converted[col])
        model_weights = model_weights.reindex(all_genes, fill_value=0)
        imat_weights = gene_to_rxn_weights(model=model_copy,gene_weights=model_weights)
        imat_results = imat(model=model_copy,rxn_weights=imat_weights,epsilon=epsilon,threshold=threshold)

        #Tabulating results
        sample_df = tabulate_model_results(model_copy, imat_results)
        sample_df = sample_df[["Key Product", "Flux (mmol/gDW/hr)"]].copy().set_index("Key Product")
        sample_df = sample_df.T.copy()
        sample_df["sample"] = col
        all_rows[f"{col}_iMAT"] = sample_df

    # Building the dataframe
    imat_data = pd.concat(all_rows.values(), axis=0,ignore_index=True)
    imat_data_merged = pd.merge(df_mini, imat_data, on = "sample").copy()

    return imat_data_merged

def iMAT_multi_integrate(model, data, upper_quantile = 0.25, lower_quantile = 0.75, epsilon=1, threshold=0.01):
    """
    Performs iMAT transcriptomic integration (using the imatpy package) to generate fluxes for a user supplied dataframe. Columns should be samples, Index should be Gene symbols, with expression as values.

    Inputs:
    - model : a cobrapy model object
    - data : a pandas dataframe with Gene Symbols as the index, Columns as the samples and values as normalised and Log2(x+1) transformed expression data.
    - upper_quantile : Defines the upper bound percentage of gene expression values for a sample to be assigned 1 for iMAT. If 0.25, the top 25% would be assigned 1
    - lower_quantile : Defines the lower bound percentage of gene expression values for a sample to be assigned -1 for iMAT. If 0.75, the bottom 25% would be assigned -1
    - epsilon : iMAT maximises the sum of high expressing reactions with flux > epsilon (default 1)
    - threshold : Alongside epsilon, iMAT maximises the sum of low expressing reactions with flux < threshold (default 0.01)
    """

    # Defining a helper function to convert expression into 1, 0 or -1 for high, neutral and low expressed genes
    def convert_col(col):
        u_q = col.quantile(upper_quantile)
        l_q = col.quantile(lower_quantile)

        converted = col.copy()
        converted[col > l_q] = 1
        converted[col < u_q] = -1
        converted[(col >= u_q) & (col <= l_q)] = 0

        return converted

    # Converting Values
    df_converted = data.apply(convert_col).copy()

    # Perform iMAT simulations for each cancer
    all_rows = {}
    all_genes = [g.id for g in model.genes]
    colnames = df_converted.columns.to_list()
    imat_counter = 0
    imat_total = len(colnames)
    for col in colnames:
        imat_counter += 1
        print(f"Simulations Performed:{imat_counter}/{imat_total}")
        model_copy = model.copy()
        model_weights = pd.Series(df_converted[col])
        model_weights = model_weights.reindex(all_genes, fill_value=0)
        imat_weights = gene_to_rxn_weights(model=model_copy,gene_weights=model_weights)
        imat_results = imat(model=model_copy,rxn_weights=imat_weights,epsilon=epsilon,threshold=threshold)

        #Tabulating results
        sample_df = tabulate_model_results(model_copy, imat_results)
        sample_df = sample_df[["Key Product", "Flux (mmol/gDW/hr)"]].copy().set_index("Key Product")
        sample_df = sample_df.T.copy()
        sample_df["Sample"] = col
        all_rows[f"{col}_iMAT"] = sample_df

    # Building the dataframe
    imat_data = pd.concat(all_rows.values(), axis=0,ignore_index=True).set_index("Sample")

    return imat_data

def iMAT_integrate(model, data, upper_quantile = 0.25, lower_quantile = 0.75, epsilon=1, threshold=0.01):
    """
    Performs iMAT transcriptomic integration (using the imatpy package) to generate fluxes for a user supplied dataframe. Columns should be samples, Index should be Gene symbols, with expression as values.

    Inputs:
    - model : a cobrapy model object
    - data : a pandas series with Gene Symbols as the index, and values as normalised and Log2(x+1) transformed expression data.
    - upper_quantile : Defines the upper bound percentage of gene expression values for a sample to be assigned 1 for iMAT. If 0.25, the top 25% would be assigned 1
    - lower_quantile : Defines the lower bound percentage of gene expression values for a sample to be assigned -1 for iMAT. If 0.75, the bottom 25% would be assigned -1
    - epsilon : iMAT maximises the sum of high expressing reactions with flux > epsilon (default 1)
    - threshold : Alongside epsilon, iMAT maximises the sum of low expressing reactions with flux < threshold (default 0.01)
    """

    # Defining a helper function to convert expression into 1, 0 or -1 for high, neutral and low expressed genes
    def convert_col(col):
        u_q = col.quantile(upper_quantile)
        l_q = col.quantile(lower_quantile)

        converted = col.copy()
        converted[col > l_q] = 1
        converted[col < u_q] = -1
        converted[(col >= u_q) & (col <= l_q)] = 0

        return converted

    # Converting Values
    df_converted = data.apply(convert_col).copy()

    # Perform iMAT simulation
    all_genes = [g.id for g in model.genes]
    model_copy = model.copy()
    model_weights = df_converted.iloc[:, 0]
    model_weights = model_weights.reindex(all_genes, fill_value=0)
    imat_weights = gene_to_rxn_weights(model=model_copy,gene_weights=model_weights)
    imat_results = imat(model=model_copy,rxn_weights=imat_weights,epsilon=epsilon,threshold=threshold)

    return imat_results