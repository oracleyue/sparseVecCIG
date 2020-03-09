# Sparse Vector CIGs

This is a fast algorithm to estimate sparse vector conditional independence
graphs, by optimizing the following l0 penalized log-likelihood:

<p align="center"> <img src="/img/log-likelihood.png" width="400"> </p> 

where S is the sample covariance matrix. The implemented algorithm is developed based on an early work
[(Marjanovic & Solo, ICASSP 2018)](https://ieeexplore.ieee.org/abstract/document/8461742).

## Algorithm description

Refer to our ICASSP 2020 for more details:

*Zuogong Yue, Padmavathi Sundaram and Victor Solo (2020). A Fast Algorithm
for Block-Wise Sparse Inverse Covariance Estimation. In 2020 IEEE
International Conference on Acoustics, Speech and Signal Processing
(ICASSP). IEEE, pp. (to be released online)*.

## List of functions and scripts

Functions: 

+ model estimation and selection

  - `spMLE.m`: the main funciton to solve the optimization problem.

  - `calcLambda.m`: chooses the "optimal" lambda (the regularization
    parameter) via AIC or BIC criterion.
    
    - `setOptions.m`: a quick/easy interface to generate the argument
      `algOpt` for `calcLambda()` calls.
    
+ random model generation and visualization
  
  - `sprandOm.m`: randomly generate sparse inverse covariance matrix Omega.
  
  - `imshowOm.m`: an assistant function to visualize matrix Omega to check
    its sparsity structure.
  
Demo scripts:

- `main.m`: a simple script to use `spMLE()` to estimate sparse Omega.

- `main_lambda.m`: a simple demo to use AIC/BIC to choose an "optimal"
  lambda by calling `calcLambda()`.
  
- `main_benchmark.m`: a script to benchmark computation time for our
  algorithm and the Goran Marjanovic's algorithm
  ([link](https://ieeexplore.ieee.org/abstract/document/8461742)). Here we
  didn't include Goran's codes due to copy right issues. One may contact
  the authors to get their codes, and place them in the folder `./goran/`
  in order to run this script.
  
  
  
*Last modified on 12 Jul 2019*
