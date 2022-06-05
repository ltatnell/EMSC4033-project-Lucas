import pytest
from python.Make_lambdas import *


def test_lambdas_to_data(tolerance = 0.001):
    
    test_lambdas = [1,10,100,1000,10000]
    
    data = lambdas_to_data(test_lambdas,radii)
    
    #these are from the program AlambdaR (https://lambdar.rses.anu.edu.au/alambdar/)
    correct_data = [2.903,2.007,1.451,1.128,0.867,0.807,0.754,0.699,0.640,0.586,0.546,0.525,0.529,0.557]
    
    diff = [a-b for a,b, in zip(data,correct_data)]
    
    #check if the difference between my modelled data is significantly different to that of AlambdaR
    #AlambdaR only goes to 3 decimal places, so tolerance here is 0.001
    assert sum(diff) < tolerance, "the data differs from the data given by AlambdaR by {}".format(sum(diff))
    
    
def test_fit_lambdas_tolerance(tolerance = 1.0):
    #tolerance in %
    
    #from alambdar
    test_data = [2.903,2.007,1.451,1.128,0.867,0.807,0.754,0.699,0.640,0.586,0.546,0.525,0.529,0.557]
    
    lambdas = fit_lambdas(test_data, N = 5)
    
    #from alambdar
    correct_lambdas = [1,10,100,1000,10000]
    
    test = [a/b for a,b in zip(lambdas[0],correct_lambdas)]

    test_difference = all((1-tolerance) <= x <= (1+tolerance) for x in test)
    
    #make sure all data doesnt differ from correct data too much
    assert test_difference, "Fitted lambdas differ too much from correct lambdas. Differences in each lambda = {}".format(test)
   

def test_fit_lambdas_different_N(tolerance = 0.00001):
    
    test_data = [2.903,2.007,1.451,1.128,0.867,0.807,0.754,0.699,0.640,0.586,0.546,0.525,0.529,0.557]
    
    N4 = fit_lambdas(test_data, N = 4)[0]
    N5 = fit_lambdas(test_data, N = 5)[0]
    
    #all of N 4 should be the same as the first 4 of N5
    test_difference = [a-b for a,b in zip(N4,N5)]
    
    test = all(a<tolerance for a in test_difference)
    
    assert test, "lambdas differ too much. Difference between first 4 lambdas when calculating 4 lambdas vs 5 lambdas is {}".format(test_difference)
    
    
def test_probability():
    
    #these data were constructed from 5 lambdas (lambda 0 to lambda 4)
    test_data = [2.903,2.007,1.451,1.128,0.867,0.807,0.754,0.699,0.640,0.586,0.546,0.525,0.529,0.557]
    
    
    test_probabilities = probability_of_N_lambdas(test_data, 1, 7, std_dev = 1)[0]
    
    max_index = test_probabilities.index(max(test_probabilities))
    
    assert max_index == 4, "Function is saying {} lambdas are mose probably fitting data constructed from 5 lambdas".format(max_index+1)