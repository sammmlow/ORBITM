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
##    Atmospheric density model based on the U.S. Standard Atmosphere 1976.  ##
##    Returns the atmospheric density (km/m^3) given an altitude input (km). ##
##    Valid only for altitudes between 86km to 1000km.                       ##
##                                                                           ##
##    Written by Samuel Y. W. Low.                                           ##
##    First created 24-May-2017 11:27 AM (+8 GMT)                            ##
##    Last modified 19-Sep-2021 22:27 PM (-7 GMT)                            ##
##                                                                           ##
###############################################################################
###############################################################################

import math

# Input R is the altitude from the surface of the Earth (km)
def density(R):
    '''Atmospheric density model based on the U.S. Standard Atmosphere 1976.
    Returns the atmospheric density (km/m^3) given an altitude input (km).
    Valid only for altitudes between 86km to 1000km.
    
    Parameters
    ----------
    R : float
        Radial altitude from the surface of the Earth (km)

    Returns
    -------
    density : float
        Atmospheric density (kg/m^3)

    '''
    
    co = [0,0,0,0,0] # Coefficient
    
    if R <= 86:
        co = [ 0.0000000000, 0.0000000000, 0.0000000000,-0.1523325,  0.202941]
    elif R <= 91:
        co = [ 0.0000000000,-3.322622e-06, 9.111460e-04,-0.2609971,  5.944694]
    elif R <= 100:
        co = [ 0.0000000000, 2.873405e-05,-0.008492037,  0.6541179, -23.62010]
    elif R <= 110:
        co = [-1.240774e-05, 0.005162063, -0.8048342,    55.55996,  -1443.338]
    elif R <= 120:
        co = [ 0.0000000000,-8.854164e-05, 0.03373254,  -4.390837,   176.5294]
    elif R <= 150:
        co = [ 3.661771e-07,-2.154344e-04, 0.04809214,  -4.884744,   172.3597]
    elif R <= 200:
        co = [ 1.906032e-08,-1.527799e-05, 0.004724294, -0.6992340,  20.50921]
    elif R <= 300:
        co = [ 1.199282e-09,-1.451051e-06, 6.910474e-04,-0.1736220, -5.321644]
    elif R <= 500:
        co = [ 1.140564e-10,-2.130756e-07, 1.570762e-04,-0.07029296,-12.89844]
    elif R <= 750:
        co = [ 8.105631e-12,-2.358417e-09,-2.635110E-06,-0.01562608,-20.02246]
    elif R <= 1000:
        co = [-3.701195e-12,-8.608611e-09, 5.118829e-05,-0.06600998,-6.137674]
    else:
        co = [0,0,0,0,0]
        
    # Compute the exponent term.
    AtmosExp  = co[0]*(R**4)
    AtmosExp += co[1]*(R**3)
    AtmosExp += co[2]*(R**2)
    AtmosExp += co[3]*(R)
    AtmosExp += co[4]*(1)
    
    # Return the atmospheric density (kg/m^3)
    return math.exp(AtmosExp)