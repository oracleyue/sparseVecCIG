# Sparse Vector CIGs

This is a faster algorithm to estimate sparse vector conditional
independence graphs, by optimizing the following $l_0$ penalized
log-likelihood:
\[
  F(\Omega) = -\log \mathrm{det}(\Omega) + \mathrm{tr}(\mathbf{S} \Omega) +
  \lambda \sum_{i \neq j} I(\Omega_{ij} \neq 0),
\]
where $\mathbf{S}$ is the sample covariance matrix.

## Algorithm description

A detailed paper is in progress.

## List of functions and scripts

Core functions:

- `bcdSpML.m`: the main funciton to solve the optimization problem.

- `calcLambda.m`: chooses the "optimal" lambda (the regularization
  parameter) via AIC or BIC criterion.
  
- `sprandOm.m`: randomly generate sparse inverse covariance matrix
  $\Omega$.
  
- `imshowOm.m`: an assistant function to visualize matrix $\Omega$ (any
  sparse matrices) to check its sparsity structure.
  
Demo scripts:

- `main.m`: a simple script to use `bcdSpML()` to estimate sparse $\Omega$.

- `main_lambda.m`: a simple demo to use AIC/BIC to choose an "optimal"
  lambda.
  
- `benchmark_single.m` and `benchmark_multiple.m`: scripts to benchmark
  computation time for our algorithm and the Goran Marjanovic's algorithm
  ([link](https://ieeexplore.ieee.org/abstract/document/8461742)). Here we
  didn't include Goran's codes due to copy right issues. One may contact
  the original authors to get their codes, and place in the folder
  `./Goran/` to run these two scripts.
  
  
  
Last modified on 12 Jul 2019