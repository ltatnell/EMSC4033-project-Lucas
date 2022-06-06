# EMSC4033-project
By Lucas Tatnell - u6079899

This repository contains code to be used as a mini-package for the following things;

1. Fitting an arbitrary number of orthogonal polynomials to a given set of data
2. Calculating probabilities for that fit, given the error in the data and the goodness of fit

Specifically, this code is catered towards the analysis of rare earth element (REE) patterns, and further contains code which can determine the most probable amount of orthogonal polynomials to describe data in a REE dataset.

Instructions on how to use the code are in Report and Instructions.md. A demonstration of the code working on a REE dataset can be found in Code/Code demo.ipynb

## What does the code do?

Orthogonal polynomial fit to data allows quantification of the constribution of increasing order independant polynomials to data. Polynomials are of the form;

$$ data(x) = \lambda_0 + \lambda_1 f_1^{orth} + \lambda_2 f_2^{orth} + \lambda_3 f_3^{orth} + ... $$

and

  $$  f_0^{orth} = 1 $$
  
  $$    f_N^{orth} = \prod_{j=1}^{N\ge 1} x - \alpha_j  $$

for any $N \ge 0$ where $\alpha_j$ are chosen such that the inner product of $f_k^{orth}$ and $f_p^{orth}$ is 0 if $k \neq p$ for the given data. 
$\lambda$ coefficients are chosen to minimise $X^2$ for a given set of data.

This mini-package allows the calculation of the realtive probability of fit for N polynomails given the error in the data. This can be used to decide how many polynomials to fit. Probabilities are calculated using 

$$ p(d|m)  \propto  {{1} \over {[(2\pi)^N*|C_{\lambda s}|]^{1/2}}  } exp(-{1\over2} * \phi(m))$$ 

where $d$ is the data, $m$ is the model, $N$ is the number of $\lambda$ fitted to data, $C_{\lambda s}$ is the covariance matrix of the lambdas, and $\phi(m)$ is the error in the data ($X^2$ statistic).
Here, the relative probability calculated is proportional to the true probability, so can be used to find the most likely model, or N.


Source code is found in Code/src/python/Make_lambdas.py  
There are three key functions provided here for probability analysis of REE patterns:

```
best_fit_anomaly(ree, N, std_dev=2)

probability_of_N_lambdas(ree, min_N, max_N, std_dev = 2):

dataset_probability(data, min_N, max_N, std_dev = 2):
```

best_fit_anomaly() calculates the best-fit anomaly for N lambdas out of Eu anomaly, Ce anomaly, both, or none for a REE pattern 'ree'.

probability_of_N_lambdas() calculates the highest probability of fit given the error in the data for lambdas between min_N and max_N. 
This uses best_fit_anomaly() to calculate lambdas

dataset_probability() runs probability_of_all_lambdas() over a pandas.DataFrame with 14 columns, one for each REE, and returns the number of lambdas most likely to describe each REE pattern. 

A number of other functions are provided, which are used to code these three functions. Using these functions, one can decide quantitatively hoe many lambda coefficients should be used to analyse data with.

