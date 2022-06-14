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
##    This routine runs the orbit maintenance simulation following the       ##
##    decay model that I had simplified in 2018. The link to my paper is     ##
##    below. This routine does not require STK or any external libraries     ##
##    besides the basic imports you see below. It also runs at a fraction    ##
##    of a second, as compared to a full STK simulation with AstroGator,     ##
##    because it does not actually compute the orbit propagation, but it     ##
##    computes the nominal decay value based on the model, and solves only   ##
##    a two-body Kepler's equation, whereas STK does the propagation with    ##
##    remarkably high fidelity and thus is expected to be more accurate.     ##
##    In general, use this offline routine if you simply need a quick ball-  ##
##    park figure for your mission's Delta-V needs. Otherwise, use STK's.    ##
##                                                                           ##
##    Link: https://digitalcommons.usu.edu/smallsat/2018/all2018/364/        ##
##                                                                           ##
##    Low, S.Y.W., &; Chia, Y.X. (2018). “Assessment of Orbit Maintenance    ##
##    Strategies for Small Satellites”, 32nd Annual AIAA/USU Conference      ##
##    on Small Satellites, Logan, Utah, Utah State University, USA.          ##
##                                                                           ##
##    Written by Samuel Y. W. Low.                                           ##
##    First created 12-Oct-2020 10:18 AM (+8 GMT)                            ##
##    Last modified 19-Sep-2021 22:27 PM (-7 GMT)                            ##
##                                                                           ##
###############################################################################
###############################################################################

# Import basic utilities
import math
import datetime
import numpy as np
import datetime

# Import local libraries
from source import atmos
from source import kepler
from source import sun_moon_pos

# To the User - for changes to program inputs, it is recommended to change 
# the parameters in the ORBITM GUI, rather than through the `config.txt` 
# file or even through the code below.

def orbmrun(tstart, tfinal, sc_Cd, sc_area_d,
            orb_a, orb_e, orb_i, orb_R, orb_w, orb_m,
            maintenance_tolerance, maintenance_fro,
            sc_mass, isp_min, isp_max):
    '''This function runs a rapid orbit maintenance scenario using a 
    simplified orbit decay model from the paper "Assessment of Orbit
    Maintenance Strategies for Small Satellites (Small Sat 2018)". The
    atmospheric density model used is the US Standard Atmosphere 1976,
    from the file "atmos.py". No full orbit propagation is performed here.
    
    Parameters
    ----------
    tstart : string
        Epoch start string (in format 1-Jan-1900-12:00:00.000) 
    tfinal : string
        Epoch final string (in format 1-Jan-1900-12:00:00.000) 
    sc_Cd : float
        Drag coefficient of spacecraft (Cd)
    sc_area_d : float
        Average drag surface area of spacecraft (m^2)
    orb_a : float
        Initial osculating semi-major axis (km)
    orb_e : float
        Initial osculating eccentricity (between 0 and 1)
    orb_i : float
        Initial osculating inclination (degrees)
    orb_R : float
        Initial osculating right ascension of ascending node (degrees)
    orb_w : float
        Initial osculating argument of periapsis (degrees)
    orb_m : float
        Initial mean anomaly (degrees)
    maintenance_tolerance : float
        The tolerance band in which a thruster fire would be triggers. For
        example, if the nominal altitude was 400km, with a 10km tolerance,
        then thrusters will be triggered when the altitude crosses 390km.
    maintenance_fro : int
        Flag for the user to perform orbit maintenance only, or orbit
        maintenance with the inclusion of a repeating ground track.
    sc_mass : float
        Wet mass of the spacecraft (including fuel, in kg)
    isp_min : float
        For propulsion sizing plot; the axis minimum for specific impulse (s).
    isp_max : float
        For propulsion sizing plot; the axis maximum for specific impulse (s).

    Returns
    -------
    sat_epochs : list
        An N-length list comprising a "linspace" vector of time in seconds
        for the entire simulation, starting at t = 0s.
    sat_altitude : list
        An N-length list comprising the geodetic altitudes of the satellite.
    sat_mean_sma : list
        An N-length list comprising the mean semi-major axis of the satellite.
    total_DV : float
        Total Delta-V used in the orbit maintenance scenario (m/s).
    total_impulse : float
        Total impulse used in the orbit maintenance scenario (kg m/s).

    '''
    
    print("You are now running ORBITM's offline orbit maintenance. \n")
    
    ##########################################################################
    ##########################################################################
    
    # Begin by preparing sub-functions for date-time conversion.
    
    months_dict = {'Jan':1, 'Feb':2, 'Mar':3, 'Apr':4,
                   'May':5, 'Jun':6, 'Jul':7, 'Aug':8,
                   'Sep':9, 'Oct':10,'Nov':11,'Dec':12}
    
    months_dict_inv = {'01':'Jan','02':'Feb','03':'Mar','04':'Apr',
                       '05':'May','06':'Jun','07':'Jul','08':'Aug',
                       '09':'Sep','10':'Oct','11':'Nov','12':'Dec'}
    
    def timestr2datetime(string):
        strSplit = string.split()
        dt = datetime.datetime(int(strSplit[2]),
                               int(months_dict[strSplit[1]]),
                               int(strSplit[0]),
                               int(strSplit[3][0:2]),
                               int(strSplit[3][3:5]),
                               int(strSplit[3][6:8]))
        return dt
    
    # Read the start epoch string as a datetime object
    tstart_dt = timestr2datetime(tstart)
    
    # Read the final epoch string as a datetime object
    tfinal_dt = timestr2datetime(tfinal)
    
    # Set the time step (default is 60 minutes)
    tstep = 3600
    tstep_dt = datetime.timedelta(seconds=3600)
    
    ##########################################################################
    ##########################################################################
    
    # Define all useful constants and sub-functions needed to compute decay.
    UNIGCONST = 6.67430e-11
    EARTHMASS = 5.97237e24      # in kg
    EARTHRADEQT = 6378.137      # in km
    MUU = UNIGCONST * EARTHMASS # GM
    
    def kep2cart(a, e, i, RAAN, omega, nu, mu):
        p = a*(1-e**2) # semi-latus rectum, [km]
        r = p/(1+e*np.cos(nu)) # orbit radius, [km]    
        # h = np.sqrt(mu*a*(1-e^2)) # angular momentum    
        x = r*(np.cos(RAAN)*np.cos(omega+nu) - np.sin(RAAN)*np.sin(omega+nu)*np.cos(i)) # x-position, [km]
        y = r*(np.sin(RAAN)*np.cos(omega+nu) + np.cos(RAAN)*np.sin(omega+nu)*np.cos(i)) # y-position, [km]
        z = r*(np.sin(i)*np.sin(omega+nu)) # z-position, [km]    
        cart = [x,y,z] # cartesian state vector
        return cart
    
    cart_state = kep2cart(orb_a, orb_e, orb_i, orb_R, orb_w, orb_m, MUU)
    
    # Keplerian orbit period, with SMA as the semi-major axis (m)
    def OrbitPeriod( SMA ):
        return 2 * math.pi * math.sqrt( (SMA**3) / MUU )
    
    # Orbit linear velocity, with R as radial distance from Earth Center (m)
    def Velocity( R, SMA ):
        return math.sqrt( MUU * ( (2/R) - (1/SMA) ) )
    
    # Deceleration from drag, with R as radial distance from Earth Center (m)
    def DragAccel( R, SMA ):
        ALT = R - ( EARTHRADEQT * 1000 )  # Altitude (m)
        A2M = sc_area_d / sc_mass         # Area to mass ratio
        VEL = Velocity( R, SMA )          # Velocity magnitude (m/s)
        ATM = atmos.density( ALT / 1000 ) # Density (kg/m^3)                
        return 0.5 * ATM * sc_Cd * A2M * (VEL**2)
    
    r_factor = 0.5                        # Reflection factor    
    def SRPAccel( r_factor ):        
        A2M = sc_area_d / sc_mass         # Area to mass ratio
        return 4.65e-6 * A2M * (1+r_factor) 
    
    def sun_moonAccel(date,cart_state):       
        cart_state = np.array(cart_state)*1e3 # Satellite state
        sun_pos = (sun_moon_pos.approxECISunPosition(date))[0]            # Sun position
        moon_pos = np.array(sun_moon_pos.approxECIMoonPosition(date))*1e3 # Moon position  
        mu_sun = 1.32712440042e20  # Solar gravitational paramater
        mu_moon = 4.9048695e12     # Lunar gravitational paramater
        # Relative position vector of satellite w.r.t. point mass 
        d_sun = cart_state - sun_pos 
        d_moon = cart_state - moon_pos 
        # Acceleration 
        a_sun = -mu_sun*(d_sun/(np.linalg.norm(d_sun)**3)+sun_pos/(np.linalg.norm(sun_pos)**3))
        a_moon = -mu_moon*(d_moon/(np.linalg.norm(d_moon)**3)+moon_pos/(np.linalg.norm(moon_pos)**3))        
        return np.linalg.norm(a_sun)+np.linalg.norm(a_moon)
        
    # The decay rate of the altitude (dR/dt) in meters (see reference).
    def DecayRate(R,SMA,datetime):
        DragAcceleration = DragAccel( R, SMA )
        SRP_Accel = SRPAccel(r_factor)
        SM_Accel = sun_moonAccel(datetime,cart_state)
        return -1 * (DragAcceleration+SRP_Accel+SM_Accel) * OrbitPeriod(SMA) / math.pi
    
    # Define a function that computes the Delta-V needed for Hohmann transfer.
    def HohmannTransferDV(rd,maintenance_tolerance,maintenance_fro):
        R1 = rd
        if maintenance_fro == True:
            R2 = rd + (2*maintenance_tolerance*1000)
        else:
            R2 = rd + (maintenance_tolerance*1000)
        DV1 = math.sqrt(MUU/R1) * (math.sqrt((2*R2)/(R1+R2)) - 1)
        DV2 = math.sqrt(MUU/R2) * (1 - math.sqrt((2*R1)/(R1+R2)))
        totalDV = DV1+DV2
        Time_Elapsed = math.pi * math.sqrt(((R1+R2)**3)/(8*MUU))
        # Returns (Delta_V (m/s), Time_Elapsed (seconds))
        return totalDV, Time_Elapsed
    
    ##########################################################################
    ##########################################################################
    
    # This section handles lifetime estimation. We initialise parameters.
    
    decay_alt = 0.0 # The decay amount in meters
    total_DV = 0.0 # Total Delta-V used (m/s)
    thrustcount = 0 # Counting the number of thrusts needed
    smaxis = orb_a * 1000 # Init the semi-major axis in meters
    meananom = np.deg2rad(orb_m) # Init the mean anomaly in radians
    
    # Now we need to setup the time sequence to run the decay simulation.
    # We will use a simple integrator without a need to propagate all 6 orbit
    # states, since solving for the radial distance is all all that is needed
    # to compute the atmospheric density.
    
    print("Running lifetime simulation now (this might take long). \n")
    
    # The life-time flag checks if satellite decays below ~86km Karman line,
    # or if the lifetime has exceeded 10 years...
    
    lifetime_flag = False
    total_duration = 0 # seconds
    
    while lifetime_flag == False and total_duration < 315360000:
        
        # Update the time step first. The lifetime computation will be sped
        # up to 10X the speed of the orbit maintenance simulation.
        tstart_dt = tstart_dt + (10 * tstep_dt)
        total_duration += 10 * tstep 
        
        # First, we need to find the mean motion of the orbit at this instant.
        kepler_period = OrbitPeriod(smaxis)
        meanmotion = (2*np.pi) / kepler_period
        
        # Next, we can compute the mean anomaly at the next time step. The 
        # complex-looking expression below is just meant to keep the mean
        # anomaly bounded within +/- 180 degrees rather than +/- 360.
        meananom = ( meananom + np.pi + ( meanmotion * 10 *tstep ) )
        meananom = meananom  % ( 2 * math.pi ) - math.pi
        
        # We can solve for the eccentric anomaly orb_E via Newton-Raphson
        # M = E - e*sin(E)
        eccnanom = kepler.solve_kepler(meananom, orb_e)
        
        # Now, we can substitute orb_E to find the satellite radial distance.
        rd = smaxis*(1-orb_e*np.cos(eccnanom))
        
        # Then, we can compute the decay rate (m/s), and apply it to the SMA.
        decay_rate = DecayRate(rd, smaxis,tstart_dt)
        decay_alt = decay_rate * (tstep)
        smaxis += decay_alt
        
        # Compute the current altitude
        current_altitude = rd/1000 - EARTHRADEQT
        
        # Check if the altitude hits below the Karman line.
        if current_altitude < 100.0:
            lifetime_flag = True
            lifetime_dt = tstart_dt
    
    # Check the status of the life time simulation.
    # lifetime_flag == True ==> orbit decay has occurred within 10 years.
    if lifetime_flag == True:
        
        # Construct the date-time string for the life time decay.
        lifetime_date = str(lifetime_dt).split()[0]
        lifetime_dd   = lifetime_date.split('-')[2]
        lifetime_mm   = lifetime_date.split('-')[1]
        lifetime_mm   = months_dict_inv[ lifetime_mm ] # Convert month string
        lifetime_yyyy = lifetime_date.split('-')[0]
        lifetime_time = str(lifetime_dt).split()[1] + '.000'
        
        # Compute the number of orbits before complete decay.
        lifetime_days = str(round( total_duration / 86400 ))
        lifetime_orbits = str(round(total_duration / OrbitPeriod(orb_a*1000)))
        
        # Construct the lifetime string to be printed out.
        lifetime_str = "Lifetime decay is estimated to be on "
        lifetime_str += lifetime_dd + ' '
        lifetime_str += lifetime_mm + ' '
        lifetime_str += lifetime_yyyy + ' '
        lifetime_str += lifetime_time + ' '
        lifetime_str += 'after ' + lifetime_orbits + ' orbits.'
        lifetime_str += 'The lifetime is ' + lifetime_days + ' days. \n'
        
        # Print the life time message.
        print(lifetime_str)
    
    # Check the status of the life time simulation.
    # lifetime_flag == False ==> orbit decay has NOT occurred within 10 years.
    if lifetime_flag == False:
        
        # Print the life time message.
        lifetime_str = 'Lifetime did not decay within the 10 year limit. \n'
        print(lifetime_str)
    
    # Reset the time.
    tstart_dt = timestr2datetime(tstart)
    total_duration = 0
    
    ##########################################################################
    ##########################################################################
    
    # Re-initialise all the parameters needed to do the run
    decay_alt = 0.0 # The decay amount in meters
    total_DV = 0.0 # Total Delta-V used (m/s)
    thrustcount = 0 # Counting the number of thrusts needed
    smaxis = orb_a * 1000 # Init the semi-major axis in meters
    meananom = np.deg2rad(orb_m) # Init the mean anomaly in radians
    
    # Initialise the lists for plotting
    sat_epochs, sat_altitude, sat_mean_sma = [], [], []
    
    # Now, we setup the same time sequence, all over again, but this time,
    # orbit maintenance manoeuvres will be included in the simulation.
    print("Running mission control sequence now (this might take long). \n")
    deltaV_file = open("output_manoeuvres.txt",'w')
    total_duration = 0 # seconds
    
    while tstart_dt <= tfinal_dt:
        
        # First, we need to find the mean motion of the orbit at this instant.
        kepler_period = OrbitPeriod(smaxis)
        meanmotion = (2*np.pi) / kepler_period
        
        # Next, we can compute the mean anomaly at the next time step.
        meananom = (meananom + (meanmotion * tstep)) % np.pi
        
        # We can solve for the eccentric anomaly orb_E via Newton-Raphson
        # M = E - e*sin(E)
        eccnanom = kepler.solve_kepler(meananom, orb_e)
        
        # Now, we can substitute orb_E to find the satellite radial distance.
        rd = smaxis*(1-orb_e*np.cos(eccnanom))
        
        # Then, we can compute the decay rate (m/s), and apply it to the SMA.
        decay_rate = DecayRate(rd, smaxis,tstart_dt)
        decay_alt = decay_rate * tstep
        smaxis += decay_alt
        
        # Append the data into the lists for plotting
        sat_epochs.append(total_duration)
        sat_altitude.append(rd/1000 - EARTHRADEQT)
        sat_mean_sma.append(smaxis/1000)
        
        # Update the time step first.
        tstart_dt = tstart_dt + tstep_dt
        total_duration += tstep
        
        # Now, we need to check if the orbit tolerance has been triggered.
        # If triggered, we will simply re-adjust the orbit as if a Hohmann
        # transfer had taken place.
        if smaxis < ((orb_a-maintenance_tolerance)*1000):
            DV, time_elapsed = HohmannTransferDV(rd,
                                                 maintenance_tolerance,
                                                 maintenance_fro)
            thrustcount += 1
            thrust_str = str(thrustcount) + " Maintain.Hohmann "
            thrust_str += str(tstart_dt) + " "
            thrust_str += str(DV) + " \n"
            deltaV_file.write(thrust_str)
            tstart_dt = tstart_dt+datetime.timedelta(seconds=int(time_elapsed))
            
            # Reset the semi-major axis
            if maintenance_fro == False:
                smaxis = orb_a*1000 
            else:
                smaxis = (orb_a + maintenance_tolerance) * 1000
            
            # Rack up the DV count.
            total_DV += DV
    
    # Close the Delta-V file
    deltaV_file.close()
    
    print("Mission successfully ran! Now extracting orbital data. \n")
    
    # Now, we need to start extracting relevant data.
    # Total impulse, inclusive of the margin of safety
    
    total_impulse = total_DV * sc_mass
    
    # Construct the summary information string objects
    impulse_str = "The total impulse needed: "
    impulse_str = impulse_str + str(total_impulse) + " \n"
    deltaV_str = "The total Delta-V (m/s) needed is "
    deltaV_str = deltaV_str + str(total_DV) + " \n"
    
    # Log the summary information
    summary_file = open("output_summary.txt",'w')
    summary_file.write(lifetime_str)
    summary_file.write(impulse_str)
    summary_file.write(deltaV_str)
    summary_file.close()
    
    return sat_epochs, sat_altitude, sat_mean_sma, total_DV, total_impulse
