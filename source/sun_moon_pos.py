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
##    Functions to calculate sun and moon position for a given UTC time.     ##
##                                                                           ##
##    Written by Ozan Kilic                                                  ##
##    Created 14-Jun-2022 10:11 AM (+1 GMT)                                  ##
##                                                                           ##
###############################################################################
###############################################################################

import datetime
import numpy as np
import math

def get_julian_datetime(date):
    """
    Convert a datetime object into julian float.
    Args:
        date: datetime-object of date in question

    Returns: float - Julian calculated datetime.
    Raises: 
        TypeError : Incorrect parameter type
        ValueError: Date out of range of equation
    """

    # Ensure correct format
    if not isinstance(date, datetime.datetime):
        raise TypeError('Invalid type for parameter "date" - expecting datetime')
    elif date.year < 1801 or date.year > 2099:
        raise ValueError('Datetime must be between year 1801 and 2099')

    # Perform the calculation
    julian_datetime = 367 * date.year - int((7 * (date.year + int((date.month + 9) / 12.0))) / 4.0) + int(
        (275 * date.month) / 9.0) + date.day + 1721013.5 + (
                          date.hour + date.minute / 60.0 + date.second / math.pow(60,
                                                                                  2)) / 24.0 - 0.5 * math.copysign(
        1, 100 * date.year + date.month - 190002.5) + 0.5

    return julian_datetime

def approxECISunPosition(UTCTime):
    
    """
    % Given Gregorian year(s), day(s), and time(s), calculate Earth-centered
    % inertial (ECI) coordinates (in meters) of the vector(s) pointing from
    % Earth's center to the center of the Sun.
    % Encapsulated by Seth B. Wagenman, from stargazing.net/kepler/sun.html as
    % corrected by https://astronomy.stackexchange.com/a/37199/34646 along with
    % other minor corrections based on recent revisions to the ecliptic plane.
    % Set time(s) relative to the J2000.0 epoch
    %
    Seth Wagenman (2022). approxECISunPosition 
    (https://www.mathworks.com/matlabcentral/fileexchange/78766-approxecisunposition), 
    MATLAB Central File Exchange. Retrieved June 14, 2022.
    """
    
    fractionalDay = get_julian_datetime(UTCTime) - 2451545.0;
    
    # Mean longitude of the Sun
    meanLongitudeSunDegrees = 280.4606184 + \
    	((36000.77005361 / 36525) * fractionalDay) # (degrees)

    # Mean anomaly of the Sun
    meanAnomalySunDegrees = 357.5277233 + \
    	((35999.05034 / 36525) * fractionalDay) # (degrees)


    meanAnomalySunRadians = meanAnomalySunDegrees * np.pi / 180

    # Ecliptic longitude of the Sun
    eclipticLongitudeSunDegrees = meanLongitudeSunDegrees + \
    	(1.914666471 * np.sin(meanAnomalySunRadians)) + \
    	(0.918994643 * np.sin(2 * meanAnomalySunRadians)) #(degrees)

    eclipticLongitudeSunRadians = eclipticLongitudeSunDegrees * np.pi / 180
    
    # Obliquity of the ecliptic plane formula mostly derived from:
    # https://en.wikipedia.org/wiki/Ecliptic#Obliquity_of_the_ecliptic
    # Formula deals with time denominated in centuries
    epsilonDegrees = 23.43929 - ((46.8093 / 3600) * fractionalDay / 36525) #(degrees)

    epsilonRadians = epsilonDegrees * np.pi / 180
    
    # **** Calculations for unit vectors direction to the Sun are in ECI ****
    # Now calculate the Earth to Sun unit vector(s) USING RADIANS, NOT DEGREES

    unitEarthSunVectorTransposed = np.array([[np.cos(eclipticLongitudeSunRadians).T],
                                             [(np.cos(epsilonRadians)*np.sin(eclipticLongitudeSunRadians)).T],
                                             [(np.sin(epsilonRadians)*np.sin(eclipticLongitudeSunRadians)).T]])

    unitEarthSunVector = unitEarthSunVectorTransposed.T

    # Scale up from unit vector(s); distances vary--Earth's orbit is elliptical
    distEarthToSunInAstronomicalUnits = \
     	1.000140612 - \
     	(0.016708617 * np.cos(meanAnomalySunRadians)) - \
     	(0.000139589 * np.cos(2 * meanAnomalySunRadians)) # in Astronomical Units
    # (new international definition for an A. U.)
    distEarthToSunInMeters = 149597870700 * distEarthToSunInAstronomicalUnits
    
    # Multiply unit vector(s) by distance to Sun for scaling to full magnitude
    # Use meters for generic code compatibility
    vector = distEarthToSunInMeters*unitEarthSunVector

    return vector

def approxECIMoonPosition(UTCTime):
    """
    Reference:
    "An alternative lunar ephemeris model for on-board flight software use",
    D. G. Simpson, Proceedings of the 1999 NASA/GSFC Flight Mechanics Symposium,
    p. 175-184).                                                      
                                                                                          
    Meysam Mahooti (2022). Moon Position 
    (https://www.mathworks.com/matlabcentral/fileexchange/56041-moon-position), 
    MATLAB Central File Exchange. Retrieved June 14, 2022.
    """
    mjd = get_julian_datetime(UTCTime) - 2400000.5;
    
    century2day = 36525
    xcoeffs = np.array([383.0e3, 31.5e3, 10.6e3, 6.2e3, 3.2e3, 2.3e3, 0.8e3])
    ycoeffs = np.array([351.0e3, 28.9e3, 13.7e3, 9.7e3, 5.7e3, 2.9e3, 2.1e3])
    zcoeffs = np.array([153.2e3, 31.5e3, 12.5e3, 4.2e3, 2.5e3, 3.0e3, 1.8e3])
    xa = np.array([8399.685, 70.990, 16728.377, 1185.622, 7143.070, 15613.745, 8467.263])
    xp = np.array([5.381, 6.169, 1.453, 0.481, 5.017, 0.857, 1.010])
    ya = np.array([8399.687, 70.997, 8433.466, 16728.380, 1185.667, 7143.058, 15613.755])
    yp = np.array([3.811, 4.596, 4.766, 6.165, 5.164, 0.300, 5.565])
    za = np.array([8399.672, 8433.464, 70.996, 16728.364, 1185.645, 104.881, 8399.116])
    zp = np.array([3.807, 1.629, 4.595, 6.162, 5.167, 2.555, 6.248])

    t = (mjd - 51544.5)/century2day # time in Julian centuries from J2000
    
    xterms = xa * t + xp
    yterms = ya * t + yp
    zterms = za * t + zp
    
    # Moon position vector in ECI, km unit
    r_moon = [np.dot(xcoeffs, np.sin(xterms)),
              np.dot(ycoeffs, np.sin(yterms)),
              np.dot(zcoeffs, np.sin(zterms))]
    
    return r_moon
