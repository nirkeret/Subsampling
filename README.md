# Subsampling

This repository has R and RCPP codes that perform the subsampling procedure described in https://arxiv.org/abs/2012.02122. 

The "**estimating function**" file contains the function that implements the procedure for a user-provided dataset. The required input arguments are listed within the code and are given here as well for convenience. **Please note that this function requires that the entire dataset can be loaded into the R session.**
Soon there will be a version that implements the reservoir sampling for a user-provided dataset divided to batches that are stored on the hard disk.
These are the arguments that the function should receive:

V = vector of observed times.

D = event status (T/F).

X = covariate matrix.

R = vector of delayed-entry (left truncation) times.

q0 = number of subsampled censored observations for the uniform pilot estimator (relevant only for the “L” and “A” methods). Defaults to be the same as q.

q = number of subsampled censored observations.

method: one of "U"/"L"/"A" (standing for uniform, L-optimal, A-optimal).

**How to use the code**: The Rmarkdown file: "**Simulations - Data generation and analysis**" contains the data generation process and the analysis workflow for the simulations described in the manuscript. After the data are generated, the estimating function is used, so it provides a clear demonstration for its usage.

The "**UKB analysis**" file contains the real data analysis workflow for the time-dependent coefficients analysis provided in the manuscript.
The analysis is demonstrated using a synthetic dataset (**synthUKB.csv**) similar in nature to the UKB (but with only 100k observations), so results obtained by analyzing it should not be similar to those in the paper. The CSV file is compressed in a zip file named "synthUKB" - please unzip the file to obtain the synthetic dataset. The time-dependent coefficient analysis uses the reservoir-sampling algorithm described in the manuscript. **Please note that this analysis does not use the "estimating function" described above**, because the latter requires that the entire dataset be loaded in the RAM. In general, the reservoir sampling algorithm requires that the data be stored in batches on the hard disk, and then the batches can be sequentially loaded. With the time-dependent coefficients analysis, since we expand the data into time-dependent form "on the fly", we don't have to store the batches on the hard disk before hand.

Currently the reservoir-sampling procedure is not implemented in a user-friendly function, but the process can be understood by the "UKB analysis" Rmarkdown file. Please note that the RCPP file required for the reservoir sampling is different than the one required for the "estimating function", which is used throughout the simulations. These RCPP codes will soon be unified and integrated into an R package, together with a friendly implemenation for the reservoir sampling algorithm.
