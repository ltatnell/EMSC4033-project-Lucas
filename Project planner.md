## EMSC4033 project plan

# A beysian approach to rare earth elements lambda coefficients

## Executive summary

Rare earth element (REE) pattern lambda parameterisation involves deconvolution of a REE pattern into a sum of orthogonal polynomials up to the fourth degree. I aim to use a beysian approach to this to decide the appropriate degree of polynomials which is most likely to describe a given REE pattern, given the error in that pattern.

## Goals

- Write an algorithm which obtains a probability of the model given the data for all polynomials up to degree N.
- Use this to show graphically whether high lambda 4 in data is reflective of real earth processes, or a result of high error.
- Investigate whether lambda 5 or higher could be observed to be a result of earth processes in high quality data.

## Background and Innovation  

The shape of a rare earth element (REE) pattern is a result of the petrogenetic processes which formed the rock which contains that REE pattern. In some cases, when minerals in which REEs are compatible are involved in petrogenesis, information about the petrogenic history can be obtained by analysing the shape of a REE pattern. Lambda coefficients for REE patterns are the coefficients to orthogonal polynomials which best describe a given REE pattern (O'neill 2016). O'neill (2016) uses conventional 95% confidence intervals to determine statistical significance of lambda coefficients. This is convenient and simple to calculate, however relies on an assumed standard deviation of uncertainty for REE measurements equal to the mean standard deviation for each measurement. However, in many cases instrumental counting uncertainty is not fully representative of the true uncertainty, which can depend on several factors. A beysian approach to this proble could provide a better way to treat uncertainty for large datasets. If both the uncertainty on measurements and the total number of polynomials fitted to the pattern is taken as variable, a probability density plot for the most likely number of polynomials can be made. This will allow one to decide how many polynomials to fit to their data which have statistical significance given the scatter in a dataset.

## Resources & Timeline



_What do you have at your disposal already that will help the project along. Did you convince somebody else to help you ? Are there already some packages you can build upon. What makes it possible to do this project in the time available. Do you intend to continue this project in the future ?_

(For example:
  - I’ll be using data of X from satellite and then also data from baby blue seals…
  - I’ll step on existing package Y and build extra functionality on top of class W.
  - I’ll use textbook Z that describes algorithms a, b, c
  - …
)

## Testing, validation, documentation

**Note:** You need to think about how you will know your code is correct and achieves the goals that are set out above (specific tests that can be implemented automatically using, for example, the `assert` statement in python.)  It can be really helpful if those tests are also part of the documentation so that when you tell people how to do something with the code, the example you give is specifically targetted by some test code.

_Provide some specific tests with values that you can imagine `assert`ing_

