# CPT_GRN_analysis

This repository contains scripts for analyzing Cell Phenotypic Transition (CPT) using Gene Regulation Network (GRN) inference. The workflow integrates RNA sequencing data and dynamic systems theory to study gene expression programs in cell transitions. Below are the steps included in the analysis:

## Workflow

1. **Gene Filtering**: Select genes that exhibit binary-like responses during cell transitions.
2. **RNA Velocity Analysis**:
   - Utilize Dynamo to calculate RNA velocity.
   - Construct transition graphs representing cell state dynamics.
3. **Trajectory Simulation**:
   - Simulate cell trajectory on the transition graph using Markov models and Fokker-Planck equations.
4. **Reaction Coordinate Calculation**:
   - Estimate reaction coordinates (RCs) that describe the transition progression across states.
5. **GRN Inference Using PLSR**:
   - Infer gene regulatory networks using Partial Least Squares Regression (PLSR).
   - Apply local false discovery rate (LFDR) to refine the network.
6. **Frustration Score Calculation**:
   - Compute the frustration score along the RCs to identify conflicts between gene regulations and their expression states.
7. **Advanced GRN Inference and Analysis**:
   - Use GRISLI for detailed GRN inference.
   - Analyze community structures within the GRN using methods like community detection.
   - Assess the network for intercommunity interactions and their implications on cell phenotype transitions.

## Additional Details

- **Data Sources**: Analysis performed on datasets from GEO, including studies on epithelial-mesenchymal transition (EMT) and development of pancreatic endocrine and dentate gyrus granule cells.
- **Statistical Analysis**: Reaction coordinates were derived from data points using methods like the finite temperature string method. Network dynamics were characterized by calculating changes in gene-gene interactions and network heterogeneity.
- **Community Detection**: Implemented to observe the dynamics within and between different gene communities, revealing how intercommunity interactions contribute to the overall gene expression reprogramming during CPTs.

## Citation

Please refer to the associated publication for detailed methodology and findings:
[paper](https://www.biorxiv.org/content/10.1101/2021.09.21.461257v2)
```
@article {Wang2021.09.21.461257,
	author = {Weikang Wang and Ke Ni and Dante Poe and Jianhua Xing},
	title = {Transiently increased intercommunity regulation characterizes concerted cell phenotypic transition},
	elocation-id = {2021.09.21.461257},
	year = {2023},
	doi = {10.1101/2021.09.21.461257},
	publisher = {Cold Spring Harbor Laboratory},
	URL = {https://www.biorxiv.org/content/early/2023/12/19/2021.09.21.461257},
	eprint = {https://www.biorxiv.org/content/early/2023/12/19/2021.09.21.461257.full.pdf},
	journal = {bioRxiv}
}
```
## License

This project is available under the [License Name]. See the LICENSE file for more info.
