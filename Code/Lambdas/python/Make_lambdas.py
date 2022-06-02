# This code provides functions to construct orthogonal polynomials in terms of the ionic radii of the 14 rare earth elements (REE)

import numpy as np
import math
from scipy.optimize import curve_fit

import sympy

test_data = (1.447,5.241,1.014,5.722,2.277,1.005,3.285,0.646,4.336,0.906,2.715,0.412,2.501,0.351)

test_data_2=(-3.086,-1.517,-0.532,0.180,1.428,2.107,2.894,3.769,4.686,5.506,6.166,6.626,6.868,6.901)

reenames  =  ("La","Ce","Pr",
              "Nd","Sm","Eu",
              "Gd","Tb","Dy",
              "Ho","Er","Tm",
              "Yb","Lu")

# Radii from O'Neill 2016
radii  =  (1.160, 1.143, 1.126,  # La Ce Pr
           1.109, 1.079, 1.066, # Nd Sm Eu
           1.053, 1.040, 1.027, # Gd Tb Dy
           1.015, 1.004, 0.994, # Ho Er Tm
           0.985, 0.977) # Yb Lu

#CI chondrite values from O'Neill 2016
ON16 = (0.2472, 0.6308, 0.0950, # La Ce Pr
          0.4793, 0.1542, 0.0592, # Nd Sm Eu
          0.2059, 0.0375, 0.2540, # Gd Tb Dy
          0.0554 ,0.1645, 0.0258, # Ho Er Tm
          0.1684, 0.0251) # Yb Lu

test_data_norm = [math.log(a/b) for a,b in zip(test_data,ON16)]

#this function is from the package pyrolite
#dont know how to install packages on here
def orthogonal_polynomial_constants(xs, degree=3, rounding=None, tol=10 ** -15):
    r"""
    Finds the parameters
    :math:`(\beta_0), (\gamma_0, \gamma_1), (\delta_0, \delta_1, \delta_2)` etc.
    for constructing orthogonal polynomial functions `f(x)` over a fixed set of values
    of independent variable `x`.
    Used for obtaining lambda values for dimensional reduction of REE data [#ref_1]_.

    Parameters
    ----------
    xs : :class:`numpy.ndarray`
        Indexes over which to generate the orthogonal polynomials.
    degree : :class:`int`
        Maximum polynomial degree. E.g. 3 will generate constant, linear, and quadratic
        polynomial components.
    tol : :class:`float`
        Convergence tolerance for solver.
    rounding : :class:`int`
        Precision for the orthogonal polynomial coefficents.

    Returns
    -------
    :class:`list`
        List of tuples corresponding to coefficients for each of the polynomial
        components. I.e the first tuple will be empty, the second will contain a single
        coefficient etc.

    Notes
    -----
        Parameters are used to construct orthogonal polymomials of the general form:

        .. math::

            f(x) &= a_0 \\
            &+ a_1 * (x - \beta) \\
            &+ a_2 * (x - \gamma_0) * (x - \gamma_1) \\
            &+ a_3 * (x - \delta_0) * (x - \delta_1) * (x - \delta_2) \\

    See Also
    --------
    :func:`~pyrolite.util.lambdas.calc_lambdas`
    :func:`~pyrolite.geochem.transform.lambda_lnREE`

    References
    ----------
    .. [#ref_1] O’Neill HSC (2016) The Smoothness and Shapes of Chondrite-normalized
           Rare Earth Element Patterns in Basalts. J Petrology 57:1463–1508.
           doi: `10.1093/petrology/egw047 <https://dx.doi.org/10.1093/petrology/egw047>`__
    """
    xs = np.array(xs)
    x = sympy.var("x")
    params = []
    for d in range(degree):
        ps = sympy.symbols("{}0:{}".format(chr(945 + d), d))

        if d:
            eqs = []
            for _deg in range(d):
                q = 1
                if _deg:
                    q = x ** _deg
                for p in ps:
                    q *= x - p
                eqs.append(q)
            sums = []
            for q in eqs:
                sumq = 0.0
                for xi in xs:
                    sumq += q.subs(dict(x=xi))
                sums.append(sumq)

            guess = np.linspace(np.nanmin(xs), np.nanmax(xs), d + 2)[1:-1]
            result = sympy.solvers.solvers.nsolve(sums, ps, list(guess), tol=tol)
            if rounding is not None:
                result = np.around(np.array(result, dtype=float), decimals=rounding)
            params.append(tuple(result))
        else:
            params.append(())  # first parameter
    return params


saved_orthogonal_polynomial_constants = {}

def get_orthogonal_polynomial_constants(xs, degree=3, rounding=None, tol=10 ** -15):
    r"""
    This retrieves saved orthogonal polynomial constants if they have already been calculated, or calculates them if they have not been calculated.
    
    from orthogonal_polynomial_constants:
    ----------
    Finds the parameters=
    :math:`(\beta_0), (\gamma_0, \gamma_1), (\delta_0, \delta_1, \delta_2)` etc.
    for constructing orthogonal polynomial functions `f(x)` over a fixed set of values
    of independent variable `x`.
    Used for obtaining lambda values for dimensional reduction of REE data [#ref_1]_.

    Parameters
    ----------
    xs : :class:`numpy.ndarray`
        Indexes over which to generate the orthogonal polynomials.
    degree : :class:`int`
        Maximum polynomial degree. E.g. 3 will generate constant, linear, and quadratic
        polynomial components.
    tol : :class:`float`
        Convergence tolerance for solver.
    rounding : :class:`int`
        Precision for the orthogonal polynomial coefficents.

    Returns
    -------
    :class:`list`
        List of tuples corresponding to coefficients for each of the polynomial
        components. I.e the first tuple will be empty, the second will contain a single
        coefficient etc.

    Notes
    -----
        Parameters are used to construct orthogonal polymomials of the general form:

        .. math::

            f(x) &= a_0 \\
            &+ a_1 * (x - \beta) \\
            &+ a_2 * (x - \gamma_0) * (x - \gamma_1) \\
            &+ a_3 * (x - \delta_0) * (x - \delta_1) * (x - \delta_2) \\

    See Also
    --------
    :func:`~pyrolite.util.lambdas.calc_lambdas`
    :func:`~pyrolite.geochem.transform.lambda_lnREE`

    References
    ----------
    .. [#ref_1] O’Neill HSC (2016) The Smoothness and Shapes of Chondrite-normalized
           Rare Earth Element Patterns in Basalts. J Petrology 57:1463–1508.
           doi: `10.1093/petrology/egw047 <https://dx.doi.org/10.1093/petrology/egw047>`__
    
    """
    
    #save a unique string for the given inputs
    save_str_name = str(xs)+str(degree)+str(rounding)+str(tol)
    
    #retrieve those constants or make new ones
    if save_str_name in saved_orthogonal_polynomial_constants:
        return saved_orthogonal_polynomial_constants[save_str_name]
    
    else:
        
        new_constants = orthogonal_polynomial_constants(xs, degree=degree, rounding=rounding, tol=tol)
        saved_orthogonal_polynomial_constants[save_str_name] = new_constants
        return(new_constants)
    
    
def lambdas_to_data(lambdas,x_data,x_fit = radii):
    r"""
    Fits N coefficients ("lambdas") of polynomials which are orthogonal at points x_fit to the data y_fit.
    
    Parameters
    ----------
    lambdas <- this is a list or tuple of orthogonal polynomial coefficients
    x_data <- this is a list of the points at which to return fitted data 
    x_fit <- this is a list of the points at which to construct orthogonal polynomials
    
    Returns
    ---------- 
    A list of modelled data at points "x_data" fitted to orthogonal polynomials where "lambdas" are the coefficients
    and the polynomials are orthogonal at "x_fit".
    """
    
    N = len(lambdas)
    
    parameters = get_orthogonal_polynomial_constants(x_fit, degree = N)
    
    x_sum = [0]*len(x_data)
    
    for n, parameter in enumerate(parameters):
        
        #same lambda for all x_data
        product = [lambdas[n]]*len(x_data)
        
        for param in parameter:
            #calculate the contribution of each lambda to each data
            mult = [data - param for data in x_data]
            
            product = [a * b for a,b in zip(product,mult)]
        
        #sum contributions of all lambdas
        x_sum = [a+b for a,b in zip(product,x_sum)]
    
    return x_sum


def fit_lambdas(y_fit,x_data = radii, x_fit = radii, N = 3, std_dev = 2):
    r"""
    Fits N coefficients ("lambdas") of polynomials which are orthogonal at points x_fit to the data y_fit.
    
    Parameters
    ----------
    y_fit <- a list of data to fit orthogonal polynomials to.
    x_data <- a list of data for the x locations of y_fit
    x_fit <- a list of data to construct orthogonal polynomials with respect to
    N <- maximum degree of orthogonal coefficients to calculate
    std_dev <- mean instrumental error for y_fit in %
    
    Returns
    ----------  
    (lambdas, covariance matrix, reduced chi squared)
    lambdas are the fitted orthogonal polynomial coefficients
    covariance matrix is the covariance of the fitted lambdas (ideally a diagonal NxN matrix)
    chi squared is the goodness of fit
    """
    
    #get parameters for N lambdas
    parameters = get_orthogonal_polynomial_constants(x_fit, degree = N)

    #create a function to calculate REE at radii x depending on N lambdas
    lambdas = []
    for n in range(0,N):
        lambdas.append('lambda' + str(n))
    
    #the number of lambdas is variable, so we define a nested function
    def fit_function(x, *lambdas):
        #this is the modelled REE at radii x
        fit_sum = 0
        
        #for each lambda we are fitting...
        for i,lam in enumerate(lambdas):
            #calculate the contribution of that lambda
    
            lam_prod = lam
            
            for parameter in parameters[i]:
            
                lam_prod = lam_prod * (x - parameter)

            #sum the contribution of each lambdas
            fit_sum += lam_prod

        #for some reason, curve_fit doesnt like the output of this function unless it is coerced to this thing
        return(np.array(fit_sum,dtype='float64'))
    
    #initial guess at the lambdas. This is just here to force curve_fit to fit the correct amount of parameters.
    p0 = [1.0]+[0]*(N-1)

    #least squares fit to N polynomials
    fit = curve_fit(fit_function,p0=p0,xdata = x_data, ydata = y_fit)
    
    lambdas = fit[0]
    
    #calculate goodness of fit chi squared statistic
    fit_data = lambdas_to_data(lambdas,x_data,x_fit)
    
    residuals = [a - b for a,b in zip(fit_data,y_fit)]
    
    residuals_sq = sum([i**2 for i in residuals])

    degrees_of_freedom = len(y_fit) - N - 1
    
    red_chi_sq = (residuals_sq/(std_dev*0.01))/degrees_of_freedom
    
    return (*fit, red_chi_sq)
    
    
def probability_of_lambdas(cov_matrix, chi_sq):
    r"""
    Calculates the probability of the model given the data from a covariance matrix and a reduced chi_squared statistic
    
    Parameters
    ----------
    cov_matrix <- a square array of covariances
    chi_sq <- chi squared statistic for the data that the covariance matrix was calculated from
    
    Returns
    ----------  
    A float of relative probability. This is self-relative and should be normalised to itself in deciding which model is the most likely.
    
    """
    
    
    N = len(cov_matrix)
    
    determinant = np.linalg.det(cov_matrix)
    
    probability = math.exp(-chi_sq/2)/(((2*math.pi)**N * determinant)**(1/2))
    
    return(probability)

    
#------ From here functions are only applicable to rare earth element patterns


def best_fit_anomaly(ree, N, std_dev=2):
    r"""
    This function decides what anomaly (Eu,Ce, both, or none) fits the data best for a rare earth element (REE) pattern based on minimising the reduced chi squared statistic.
    
    Parameters
    ----------
    ree <- this should be a length 14 iterable in order of the REE, with missing data set to None. REE should be of the form ln(REE/REE-CI) where REE-CI is the REE contents in CI-chondrite.
    N <- number of orthogonal polynomials to fit to.
    std_dev <- the mean instrumental error of the data in %.
    
    Returns
    ----------
    lambdas, covariance matrix, reduced chi squared, anomaly
    
    lambdas are coefficients to N orthogonal polynomials to fit the data
    covariance matrix is the covariance of those lambdas
    reduced chi squared is the reduced chi squared statistic for lambdas 
    anomaly is one of 0,1,2,3. 0 means no anomaly was fitted. 1 means Eu anomaly fitted. 2 means Ce anomaly fitted. 3 means Eu and Ce anomaly fitted.
    """
    
    #find locations of Nones
    #make filters for possible anomalies
    none_indexes = {i for i,v in enumerate(ree) if v == None}
    
    filt_Eu = none_indexes | {5}
    filt_Ce = none_indexes | {1}
    filt_Eu_Ce = none_indexes | {5,1}
    
    radii_none = [v for i,v in enumerate(radii) if i not in none_indexes]
    ree_none = [v for i,v in enumerate(ree) if i not in none_indexes]
    
    radii_Eu = [v for i,v in enumerate(radii) if i not in filt_Eu]
    ree_Eu = [v for i,v in enumerate(ree) if i not in filt_Eu]
    
    radii_Ce = [v for i,v in enumerate(radii) if i not in filt_Ce]
    ree_Ce = [v for i,v in enumerate(ree) if i not in filt_Ce]
    
    radii_Eu_Ce = [v for i,v in enumerate(radii) if i not in filt_Eu_Ce]
    ree_Eu_Ce = [v for i,v in enumerate(ree) if i not in filt_Eu_Ce]
 
    #calculate lambdas and chi squared for all anomalies
    lambdas_none = fit_lambdas(ree_none,x_data = radii_none, N = N, std_dev = std_dev)
    lambdas_Eu = fit_lambdas(ree_Eu,x_data = radii_Eu, N = N, std_dev = std_dev)
    lambdas_Ce = fit_lambdas(ree_Ce,x_data = radii_Ce, N = N, std_dev = std_dev)
    lambdas_Eu_Ce = fit_lambdas(ree_Eu_Ce,x_data = radii_Eu_Ce, N = N, std_dev = std_dev)
    
    
    #take the lambdas with the lowest chi squared
    chi_squareds = [lambdas_none[2],lambdas_Eu[2],lambdas_Ce[2],lambdas_Eu_Ce[2]]
    
    lowest_chi_sq_index = chi_squareds.index(min(chi_squareds))
    
    lambdas_out = [lambdas_none,lambdas_Eu,lambdas_Ce,lambdas_Eu_Ce][lowest_chi_sq_index]
    
    return(*lambdas_out, lowest_chi_sq_index)

    
def probability_of_N_lambdas(ree, min_N, max_N, std_dev = 2):
    r"""
    Calculates the self-relative probabilities for fitting N lambdas to rare earth element (REE) data. N iterates from min_N to max_N.
    
    Parameters
    ----------
    ree <- iterable of 14 REE with missing data set to None
    min_N <- minimum orthogonal polynomial degree to calculate
    max_N <- maximum orthogonal polynomial degree to calculate
    std_dev <- instrumental error standard deviation in %.
    
    Returns
    ----------
    (probabilities),(N)
    
    probabilities is a tuple of self-relative probabilities relating to N polynomials fitted to data. N is a tuple of which degree polynomial each probability relates to.
    """
    
    probabilities = [None] * (max_N-min_N+1)
    
    for N in range(min_N, max_N+1):
        
        lambdas = best_fit_anomaly(ree, N = N, std_dev = std_dev)
        
        probabilities[N-min_N] = probability_of_lambdas(lambdas[1],lambdas[2])
    
    #normalise to self
    probabilities_sum = sum(probabilities)
    
    relative_probabilities = [i/probabilities_sum for i in probabilitites]
    
    return(tuple(relative_probabilities), tuple([i for i in range(min_N, max_N+1)]))
    
    
    
    


    