from .model_preparation import *
from .model_analysis import *
from .transcriptomic_integration import *

import sys
import json

from huggingface_hub import hf_hub_download
import joblib
import numpy as np
import pandas as pd

from skorch.classifier import NeuralNetClassifier
import torch
import torch.nn as nn
from scipy.special import expit

# Definingh FocalLoss for ANN
class FocalLoss(nn.Module):
    def __init__(self, alpha=0.25, gamma=5, reduction='mean'):
        super().__init__()
        self.alpha = alpha
        self.gamma = gamma
        self.reduction = reduction
    
    def forward(self, logits, targets):
        targets = targets.view(-1,1).type_as(logits)
        bce_loss = nn.functional.binary_cross_entropy_with_logits(logits, targets, reduction='none')
        pt = torch.exp(-bce_loss)
        focal_loss = self.alpha * (1 - pt) ** self.gamma * bce_loss
        if self.reduction == 'mean':
            return focal_loss.mean()
        elif self.reduction == 'sum':
            return focal_loss.sum()
        else:
            return focal_loss
        
# Class for the ANN model (binary classification)
class DeepBinary(nn.Module):
    def __init__(self, hidden_dim=64, num_layers=4, dropout_rate=0.25):
        super().__init__()
        layers = []
        layers.append(nn.LazyLinear(hidden_dim))
        layers.append(nn.LayerNorm(hidden_dim))
        layers.append(nn.ReLU())
        layers.append(nn.Dropout(dropout_rate))

        for _ in range(num_layers - 1):
            layers.append(nn.Linear(hidden_dim, hidden_dim))
            layers.append(nn.LayerNorm(hidden_dim))
            layers.append(nn.ReLU())
            layers.append(nn.Dropout(dropout_rate))

        layers.append(nn.Linear(hidden_dim, 1))  # final logit
        self.model = nn.Sequential(*layers)

    def forward(self, x):
        return self.model(x)

# Defining NeuralNet class which will be necessary for use with skorch and skopt
class NeuralNetBinaryClassifier(NeuralNetClassifier):
    def predict_proba(self, X):
        logits = self.forward(X).detach().cpu().numpy()
        probs = expit(logits)
        return np.hstack((1 - probs, probs))
    
def predict_proba_from_module(module: nn.Module, X: np.ndarray, device: str = "cpu") -> np.ndarray:
    """Return Nx2 array like sklearn predict_proba (cols [p0,p1])"""
    module = module.to(device)
    module.eval()
    with torch.no_grad():
        X_t = torch.from_numpy(np.asarray(X, dtype=np.float32)).to(device)
        logits = module(X_t).cpu().numpy().ravel()  # shape (N,)
        probs1 = expit(logits)  # sigmoid
        probs0 = 1.0 - probs1
        return np.vstack([probs0, probs1]).T

#Collecting models from hugging face

repo_id = "JackWJW/LGG_Prognosis_Ensemble"

svm_file = hf_hub_download(repo_id=repo_id, filename="LGG_Prediction_Models/models/SVM_pipeline.joblib")
rf_file  = hf_hub_download(repo_id=repo_id, filename="LGG_Prediction_Models/models/RandomForest_pipeline.joblib")
lr_file  = hf_hub_download(repo_id=repo_id, filename="LGG_Prediction_Models/models/LogisticRegression_pipeline.joblib")
xgb_file = hf_hub_download(repo_id=repo_id, filename="LGG_Prediction_Models/models/XGBoost.joblib")

ann_state_file   = hf_hub_download(repo_id=repo_id, filename="LGG_Prediction_Models/models/ANN_state_dict.pt")
ann_config_file  = hf_hub_download(repo_id=repo_id, filename="LGG_Prediction_Models/models/ANN_config.json")
ann_scaler_file  = hf_hub_download(repo_id=repo_id, filename="LGG_Prediction_Models/models/ANN_preprocessor_scaler.joblib")  # optional
ann_vt_file      = hf_hub_download(repo_id=repo_id, filename="LGG_Prediction_Models/models/ANN_preprocessor_varthresh.joblib")  # optional

scaler_file = hf_hub_download(repo_id=repo_id, filename="LGG_Prediction_Models/ensemble/scaler.joblib")
beta_file   = hf_hub_download(repo_id=repo_id, filename="LGG_Prediction_Models/ensemble/beta_vec.npy")
ens_thresh_file = hf_hub_download(repo_id=repo_id, filename="LGG_Prediction_Models/ensemble/ens_threshold.json")

# Load them
svm = joblib.load(svm_file)
rf  = joblib.load(rf_file)
lr  = joblib.load(lr_file)
xgb = joblib.load(xgb_file)

# ANN: load config + state_dict, reconstruct module
ann_config = json.load(open(ann_config_file))
ann_module = DeepBinary(
    hidden_dim=ann_config["module_args"].get("hidden_dim"),
    num_layers=ann_config["module_args"].get("num_layers"),
    dropout_rate=ann_config["module_args"].get("dropout_rate"),
)
# load state dict (CPU safe)
state = torch.load(ann_state_file, map_location="cpu")
ann_module.load_state_dict(state)

# optional: load ANN preprocessor(s)
try:
    ann_scaler = joblib.load(ann_scaler_file)
except Exception:
    ann_scaler = None
try:
    ann_vt = joblib.load(ann_vt_file)
except Exception:
    ann_vt = None

# ensemble components
scaler = joblib.load(scaler_file)
beta = np.load(beta_file)
ens_thresh = json.load(open(ens_thresh_file)).get("ensemble_threshold")

def calculate_GSL_score(data):
    """
    Calculates the GSL risk scores based on RNA-seq data (TPM input format).

    Inputs:
    - data : a pandas dataframe with Gene Symbols as the index, Columns as the samples and values as TPM expression data.
    """

    #Selecting genes of interest
    GENE_LIST = ['A4GALT', 'ABO', 'B3GALNT1', 'B3GALT1', 'B3GALT4', 'B3GALT5', 
    'B3GNT2', 'B3GNT5', 'B4GALNT1', 'B4GALT5', 'B4GALT6', 
    'FUT1', 'FUT2', 'FUT3', 'FUT5', 'FUT6', 'FUT9', 
    'GAL3ST1', 'GCNT2', 'ST3GAL1', 'ST3GAL2', 'ST3GAL3', 'ST3GAL4', 
    'ST3GAL5', 'ST3GAL6', 'ST6GALNAC2', 'ST6GALNAC3', 'ST6GALNAC4', 'ST6GALNAC5', 
    'ST6GALNAC6', 'ST8SIA1', 'ST8SIA5', 'UGCG', 'UGT8']

    filter_df = data.loc[GENE_LIST]

    #Log2 normalising the data
    df_input = np.log2(filter_df+0.001)

    #downloading model
    model = download_GSL_model()

    #Performing simulation
    sim_results = iMAT_multi_integrate(model=model, data=df_input,upper_quantile=0.3,lower_quantile=0.7,epsilon=100,threshold=10)
    sim_df = sim_results[0]

    #Creating combined df
    combined_df = pd.merge(df_input.T,sim_df,left_index=True,right_index=True)

    #Generating predictions
    probs = {}
    probs["SVM"] = svm.predict_proba(combined_df)[:,1]
    probs["RandomForest"] = rf.predict_proba(combined_df)[:,1]
    probs["XGBoost"] = xgb.predict_proba(combined_df)[:,1]
    probs["LogisticRegression"] = lr.predict_proba(combined_df)[:,1]

    ann_X = combined_df.copy()
    if ann_vt is not None:
        ann_X = ann_vt.transform(ann_X)
    if ann_scaler is not None:
        ann_X = ann_scaler.transform(ann_X)
    # get ANN probs using the pure PyTorch module
    ann_probs2col = predict_proba_from_module(ann_module, np.asarray(ann_X, dtype=np.float32), device="cpu")  # Nx2
    probs["ANN"] = ann_probs2col[:, 1]

    #Ensemble calculation
    model_list = ["SVM","RandomForest","XGBoost","LogisticRegression","ANN"]
    L_res = np.vstack([probs[m] for m in model_list]).T
    Z_res = scaler.transform(L_res)
    eta_res = Z_res @ beta

    results_dict = {}
    results_dict["GSL_Score"] = eta_res

    ensemble_preds = (eta_res >= ens_thresh).astype(int)
    results_dict["GSL_Class"] = ensemble_preds

    #Generating results dataframe
    results_df = combined_df.copy()
    results_df["GSL_Score"] = results_dict["GSL_Score"]
    results_df["GSL_Class"] = results_dict["GSL_Class"]
    class_dict = {1:"High", 0:"Low"}
    results_df["GSL_Risk"] = results_df["GSL_Class"].map(class_dict)

    return results_df


