###########################################################
###########################################################
###    Created on Wed May 24 11:27:54 2017              ###
###    Updated on Thu May 25 13:36:15 2017              ###
###    By Samuel Low                                    ###
###    Atmospheric Density Model                        ###
###    U.S. Standard Atmosphere Table 1976              ###
###    Valid only for altitudes 86km to 1000km          ###
###########################################################
###########################################################

import math

# Given some mean anomaly, M, find the eccentric anomaly E from the relation
# M = E - e*sin(E), where M is input in radians.

def SolveKepEqn(M,e):
    E1 = M # Initialise eccentric anomaly
    residual = 1.0 # Initialise convergence residual
    while residual >= 0.00001:
        fn = E1 - (e*math.sin(E1)) - M
        fd = 1 - (e*math.cos(E1))
        E2 = E1 - (fn/fd)
        residual = abs(E2-E1) # Compute residual
        E1 = E2 # Update the eccentric anomaly
    
    return E2