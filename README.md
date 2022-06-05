# EMSC4033-project
By Lucas Tatnell - u6079899

This repository contains code to be used as a library for the following things;

1. Fitting an arbitrary number of orthogonal polynomials to a given set of data
2. Calculate probabilities for that fit, given the error in the data and the goodness of fit

Specifically, this code is catered towards the analysis of rare earth element (REE) patterns, and further contains code which can determine the most probable amount of orthogonal polynomials to describe data in a REE dataset.

Instructions on how to use the code are in Report and Instructions.md. A demonstration of the code working on a REE dataset can be found in Code/Code demo.ipynb

## What does the code do?

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

