# CPT-GRN Documentation

## 1. Gene Binarization (`1_gene_binarization.ipynb`)

This notebook is designed for selecting genes that exhibit bimodal distributions in their corresponding processes. A gene is considered binarizable if:
1. K-means clustering yields a high silhouette score.
2. Any of the split values lie between the centers of the two clusters (sv > c1 and sv < c2).

Silhouette scores cannot identify unimodal distributions, and Silverman's bandwidth method is sensitive to data fluctuations, potentially leading to spurious modes.

**Parameters:**
- `cl_score`: Silhouette score of K-means clustering for filtering genes based on the bimodality of each gene’s distribution.

**Return:**
- Genes with high variance that show significant bimodal distribution and can be binarized.

## 2. RNA Velocity with Dynamo (`2_rna_velocity_Dynamo.ipynb`)

This script calculates single-cell RNA velocity using Dynamo. The function `dynamo.tl.cell_velocities` projects high-dimensional velocity vectors onto given low-dimensional embeddings and computes the cell transition (fp_transition_matrix).

**Return:**
- Processed data with Dynamo, including RNA velocity and Fokker-Planck transition matrix between different cells.

## 3. Trajectory Simulation (`3_trajectory_simulation.ipynb`)

Based on the fp_transition_matrix from the previous step, single-cell trajectories are simulated using the transition matrix. The mean value of the first sample population (cell type) is calculated, and cells are selected around this center. In each simulation, the trajectory starts from one of these cells and transitions to the next cell based on its transition probability in the transition graph. Trajectories ending in the final sample population are used to calculate reaction coordinates (RCs).

**Parameters:**
- `dwell_thres`: The threshold duration for trajectories in the target cell state. Trajectories exceeding this threshold are defined as ‘reactive trajectories’ and used for calculating reaction coordinates.

**Return:**
- Simulated ‘reactive’ single-cell trajectories.

## 4. Reaction Coordinate Calculation (`4_reaction_coordinate.ipynb`)

Reaction coordinates are calculated using both the simulated trajectories and all single-cell data.

**Parameters:**
- `nrc`: Number of reaction coordinates.
- `traj_w`: Weights of simulated trajectories in computing reaction coordinates.
- `sf`: Smoothing factor for B-spline interpolation.

**Return:**
- Reaction coordinates (array: nrc × n_PCs) of the CPT process.

## 5. Frustration Calculation with PLSR (`5_frustration_with_plsr.ipynb` and `6_locdfr.ipynb`)

These notebooks calculate the frustration score of gene regulation networks (GRNs) in each single cell. The GRN is inferred using PLSR, and statistically significant regulatory relationships are selected using locfdr. The single-cell RNA vector is binarized to calculate frustration.

**Parameters:**
- `fdrate`: False discovery rate for identifying significant regulations.

**Return:**
- Sparse GRN saved in `csr_matrix` format.
- Frustration score of the GRN.

## 6. Community Analysis (`7_community_analysis.ipynb`)

This notebook uses the Leiden algorithm to divide the GRN into different communities and calculate inter-community regulations during CPT.

**Parameters:**
- `res_v`: Resolution of community detection. Higher values lead to more communities. Default value is 1.0.

**Return:**
- Inter-community interactions of GRN in each single cell and a scatter plot.

## 7. Gene Correlation Binarization (`binarize_gene_corr.ipynb`)

This notebook sets the threshold for genes that cannot be binarized in the first step. Genes binarized using the method in the first step serve as references. By correlating with these reference genes, genes that could not be binarized initially can now be selected and binarized.

**Parameters:**
- `high_thres`: Threshold value array for the proportion of binarized genes correlated with the current gene in the ON state. Higher values improve binarization of the current gene.
- `low_thres`: Threshold value array for the proportion of binarized genes correlated with the current gene in the OFF state. Lower values improve binarization of the current gene.

**Return:**
- Genes that can be binarized using the correlation method and the centers of both clusters for binarization.

## 8. Frustration Calculation with GRISLI (`frustration_with_GRISLI.ipynb`)

This notebook uses GRISLI as an alternative method for inferring gene regulation networks and calculating frustration.

**Parameters:**
- `tig_score_thres`: Threshold of Tigress score for identifying significant regulations.

**Return:**
- Frustration score of the GRN inferred with GRISLI.

## 9. Jacobian Matrix Analysis (`jacobian_stat.ipynb`)

This notebook studies the variation of the Jacobian matrix during CPT. The Jacobian matrix is calculated using Dynamo.

**Return:**
- Jacobian matrix of each single cell and statistics of single-cell Jacobian matrices.