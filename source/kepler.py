###############################################################################
###############################################################################
##                                                                           ##
##     _____ ___  ___  ___  _____      __  __                                ##
##    |  _  | _ \| _ \|_ _||_   _|    |  \/  |                               ##
##    | |_| |   <| _ < | |   | |   _  | \  / |                               ##
##    |_____|_|\_|___/|___|  |_|  |_| |_|\/|_|                               ##
##                                                     v 1.1                 ##
##                                                                           ##
##    FILE DESCRIPTION:                                                      ##
##                                                                           ##
##    A solver for Kepler's equation; taking in the mean anomaly and the     ##
##    orbit eccentricity as input, the program solves for the eccentric      ##
##    anomaly of the spacecraft via a Newton-Raphson approach.               ##
##                                                                           ##
##    Written by Samuel Y. W. Low.                                           ##
##    First created 13-Jan-2021 10:34 AM (+8 GMT)                            ##
##    Last modified 13-Jan-2021 10:34 AM (+8 GMT)                            ##
##                                                                           ##
###############################################################################
###############################################################################

import math

def solve_kepler(M,e):
    '''Given some mean anomaly, M (rad), and eccentricity, e, this function 
    finds the eccentric anomaly E from the relation M = E - e*sin(E).
    
    Parameters
    ----------
    M : float
        Mean anomaly (radians) between +/- pi
    e : float
        Osculating eccentricity between (not inclusive) 0 and 1

    Returns
    -------
    E : float
        Eccentric anomaly (radians)

    '''
    
    E1 = M # Initialise eccentric anomaly
    residual = 1.0 # Initialise convergence residual
    while residual >= 0.00001:
        fn = E1 - (e*math.sin(E1)) - M
        fd = 1 - (e*math.cos(E1))
        E2 = E1 - (fn/fd)
        residual = abs(E2-E1) # Compute residual
        E1 = E2 # Update the eccentric anomaly
        
    return E2