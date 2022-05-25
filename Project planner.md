## EMSC4033 project plan

# A beysian approach to rare earth elements lambda coefficients

## Executive summary

Rare earth element (REE) pattern lambda parameterisation involves deconvolution of a REE pattern into a sum of orthogonal polynomials up to the fourth degree. I aim to use a beysian approach to this to decide the appropriate degree of polynomials which is most likely to describe a given REE pattern, given the error in that pattern.

## Goals

- Write an algorithm which obtains a probability of the model given the data for all polynomials up to degree N.
- Use this to show graphically whether high lambda 4 in data is reflective of real earth processes, or a result of high error.
- Investigate whether lambda 5 or higher could be observed to be a result of earth processes in high quality data.

## Background and Innovation  

The shape of a rare earth element (REE) pattern is a result of the petrogenetic processes which formed the rock which contains that REE pattern. In some cases, when minerals in which REEs are compatible are involved in petrogenesis, information about the petrogenic history can be obtained by analysing the shape of a REE pattern. Lambda coefficients for REE patterns are the coefficients to orthogonal polynomials which best describe a given REE pattern (O'neill 2016). O'neill (2016) uses conventional 95% confidence intervals to determine statistical significance of lambda coefficients. This is convenient and simple to calculate, however relies on an assumed standard deviation of uncertainty for REE measurements equal to the mean standard deviation for each measurement. However, in many cases instrumental counting uncertainty is not fully representative of the true uncertainty, which can depend on several factors. This "mean" standard deviation is inferred by setting the X<sup>2</sup> function equal to the degrees of freedom for each point in a dataset and taking the mean. A beysian approach to this proble could provide a better way to treat uncertainty for large datasets. If both the uncertainty on measurements and the total number of polynomials fitted to the pattern is taken as variable, a probability density plot for the most likely number of polynomials can be made. This will allow one to decide how many polynomials to fit to their data which have statistical significance given the scatter in a dataset.

To do this, orthogonal polynomial constants for polynomials up to degree N must be calculated. A function to do this is available in the python package "Pyrolite" (Morgans et. al., 2021). I will then randomly sample a "mean" standard deviation given the distribution of standard deviations in a dataset for polynomials up to degree N. A probability of fit can then be calculated for polynomial of degree N with randomly sampled standard deviation. If this is less than the probability of the previous randomly sampled data, then it may be randomly discarded. Running this algorithm an arbitrary amount of times allows us to create a probability density function with respect to the number of polynomials fitted to the data.

## Resources & Timeline

The python package "Pyrolite" (Morgans et. al., 2021) has a function to calculate orthogonal polynomial constants, and my program "ClambdaR" has functions for fitting data to a set of orthogonal polynomials. I will have to create a function to randomly sample a "mean" standard deviation from the dataset, given polynomials of degree N. I will then create a function to calculate the probability of fit for polynomails of degree N with a randomply sampled standard deviation. I will then create a function which randomly samples N and standard deviations an arbitrary number of times to obtain a probability density plot for N, the number of polynomials most likely to describe the data given the error in the data (this is similar to transdimensional inferrence, eg. Sambridge et. al. 2013). 
  
I estimate it will take a week to create these functions, and another week to debug and validate the results.

## Testing, validation, documentation

I'm a little unsure at this stage how I am going to test this.
I can test the generalised functions which create lambda coefficients against the conventional method of calculating lambdas to 4th degree.
Will have to think about how to test whether my probablilities are real.



  

## References

Morgan, W., Tom, ChetanNathwani, NicolasPietteLauziere, Martin, B., Alessandro, G., Louise, S., and RichardScottOz, 2021, morganjwilliams/pyrolite: 0.3.1.  
O’Neill, H. S. C., 2016, The Smoothness and Shapes of Chondrite-normalized Rare Earth Element Patterns in Basalts: Journal of Petrology, v. 57, no. 8, p. 1463-1508.  
Sambridge, M., Bodin, T., Gallagher, K., and Tkalcic, H., 2013, Transdimensional inference in the geosciences: Philos Trans A Math Phys Eng Sci, v. 371, no. 1984, p. 20110547.  

