# Subsampling

This repository has R and RCPP codes that perform the subsampling procedure described in https://arxiv.org/abs/2012.02122. 

The "**estimating function**" file contains the function that implements the procedure for a user-provided dataset. The required input arguments are listed within the code and are given here as well for convenience.
*Please note that this function requires that the entire dataset can be loaded into the R session.
Soon there will be a version that implements the reservoir sampling for a user-provided dataset divided to batches that are stored on the hard drive.*
These are the arguments that the function should receive:

V = observed times

D = status (T/F)

X = covariate matrix

R = recruitment (left truncation) times

q0 = number of subsampled censored observations for the uniform estimator (being the pilot estimator for L and A). Defaults to be the same as q.

q = number of subsampled censored observations in the second round for the L and A 

method: one of "U"/"L"/"A" (standing for uniform, L-optimal, A-optimal)

**As for the other files**:
The "**Simulations - Data generation and analysis**" file contains the data generation process and the analysis workflow for the simulations described in the manuscript.

The "**UKB analysis**" file contains the real data analysis workflow - both the main time-dependent coefficients analysis, as well as the cross-validated AUC(t) derivations.
The analysis is demonstrated using a synthetic dataset (**synthUKB.csv**) similar in nature to the UKB, so results obtained by analyzing it should not be similar to those in the paper.

The RCPP files are required for performing the subsampling procedure. For the reservoir sampling there is a slight difference in one of the functions, hence the two different
files. These codes will soon be unified and integrated into a package.
