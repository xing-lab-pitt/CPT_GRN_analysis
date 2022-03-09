import scvelo as scv
import pandas as pd 
import numpy as np
from anndata import AnnData
import loompy
from matplotlib import pyplot as plt
from sklearn.preprocessing import StandardScaler,MinMaxScaler
import matplotlib.patches as mpatches
import pickle
import os
import mnnpy
import leidenalg
from sknetwork.clustering import Louvain, BiLouvain, modularity, bimodularity
from sknetwork.visualization import svg_graph, svg_digraph, svg_bigraph
from scipy import sparse
from IPython.display import SVG
import community
import gseapy as gp
import autograd
from autograd import grad, jacobian
from sklearn.linear_model import LinearRegression
from scipy import stats,signal
from community import community_louvain
from sklearn.preprocessing import StandardScaler,MinMaxScaler
from sklearn.cross_decomposition import PLSRegression

from sklearn.feature_selection import f_regression, mutual_info_regression
import scanpy

from sklearn.mixture import GaussianMixture

from sklearn.metrics import silhouette_score
from scipy.sparse import issparse
from scvelo.preprocessing.utils import get_mean_var,materialize_as_ndarray





def filter_dispersion(adata,n_bins=20,min_disp = 0.5,max_disp = np.inf,min_mean = 0.1,max_mean = 3):
    
    mean, var = materialize_as_ndarray(get_mean_var(adata.X))

    mean[mean == 0] = 1e-12  # set entries equal to zero to small value
    dispersion = var / mean
    # logarithmized mean as in Seurat
    dispersion[dispersion == 0] = np.nan
    dispersion = np.log(dispersion)
    mean = np.log1p(mean)
    
    
    
    df = pd.DataFrame()  
    
    df["mean"]=mean
    df["dispersion"] = dispersion

    df["mean_bin"] = pd.cut(df["mean"], bins=n_bins)
    disp_grouped = df.groupby("mean_bin")["dispersion"]

    disp_mean_bin = disp_grouped.mean()
    disp_std_bin = disp_grouped.std(ddof=1)


    # retrieve genes that have nan std (i.e. single gene fell in one bin)
    # and implicitly set them to have a normalized disperion of 1
    one_gene_per_bin = disp_std_bin.isnull()

    disp_std_bin[one_gene_per_bin] = disp_mean_bin[one_gene_per_bin].values
    disp_mean_bin[one_gene_per_bin] = 0

    # normalized dispersion
    mu = disp_mean_bin[df["mean_bin"].values].values
    std = disp_std_bin[df["mean_bin"].values].values
    df["dispersion_norm"] = ((df["dispersion"] - mu) / std).fillna(0)
    
    dispersion_norm = df["dispersion_norm"].values
    gene_subset = np.logical_and.reduce(
        (
            mean > min_mean,
            mean < max_mean,
            dispersion_norm > min_disp,
            dispersion_norm < max_disp,
        )
    )

    adata.var["means"] = df["mean"].values
    adata.var["dispersions"] = df["dispersion"].values
    adata.var["dispersions_norm"] = df["dispersion_norm"].values
    adata.var["highly_variable"] = gene_subset
    
    return adata,df,gene_subset