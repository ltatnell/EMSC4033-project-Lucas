# This code provides functions to construct orthogonal polynomials in terms of the ionic radii of the 14 rare earth elements (REE)

reenames <- c("La",
              "Ce",
              "Pr",
              "Nd",
              "Sm",
              "Eu",
              "Gd",
              "Tb",
              "Dy",
              "Ho",
              "Er",
              "Tm",
              "Yb",
              "Lu")

# Radii from O'Neill 2016
radii <- c(1.160, # La
           1.143, # Ce
           1.126, # Pr
           1.109, # Nd
           1.079, # Sm
           1.066, # Eu
           1.053, # Gd
           1.040, # Tb
           1.027, # Dy
           1.015, # Ho
           1.004, # Er
           0.994, # Tm
           0.985, # Yb
           0.977) # Lu

names(radii) <- reenames

def orthogonal_polynomial_constants(xs, degree=3, rounding=None, tol=10 ** -14){

    xs = np.array(xs)
    x = var("x")
    params = []
    for d in range(degree){
        ps = symbols("{}0:{}".format(chr(945 + d), d))
        logger.debug("Generating {} DIM {} equations for {}.".format(d, d, ps))
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
        }
    return params
}
