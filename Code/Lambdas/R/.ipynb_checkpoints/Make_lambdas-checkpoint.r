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





