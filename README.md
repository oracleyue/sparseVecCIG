# Sparse Vector CIGs

This is a fast algorithm to estimate sparse vector conditional independence
graphs, by optimizing the following l0 penalized log-likelihood: <p
align="center"> <img
src="https://raw.githubusercontent.com/oracleyue/sparseVecCIG/master/img/log-likelihood.png"
width="400"> </p> where S is the sample covariance matrix. The implemented
algorithm is developed based on an early work [(Marjanovic & Solo, ICASSP
2018)](https://ieeexplore.ieee.org/abstract/document/8461742).

## Algorithm description

A detailed article is in preparation.

## List of functions and scripts

Core functions:

- `bcdSpML.m`: the main funciton to solve the optimization problem.

- `calcLambda.m`: chooses the "optimal" lambda (the regularization
  parameter) via AIC or BIC criterion.
  
- `sprandOm.m`: randomly generate sparse inverse covariance matrix Omega.
  
- `imshowOm.m`: an assistant function to visualize matrix Omega to check
  its sparsity structure.
  
Demo scripts:

- `main.m`: a simple script to use `bcdSpML()` to estimate sparse Omega.

- `main_lambda.m`: a simple demo to use AIC/BIC to choose an "optimal"
  lambda.
  
- `benchmark_single.m` and `benchmark_multiple.m`: scripts to benchmark
  computation time for our algorithm and the Goran Marjanovic's algorithm
  ([link](https://ieeexplore.ieee.org/abstract/document/8461742)). Here we
  didn't include Goran's codes due to copy right issues. One may contact
  the original authors to get their codes, and place in the folder
  `./Goran/` to run these two scripts.
  
  
  
*Last modified on 12 Jul 2019*