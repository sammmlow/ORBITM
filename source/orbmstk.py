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
##    This function runs a full orbit maintenance scenario using STK         ##
##    Astrogator (either v10 or v11). Note that the user must have a valid   ##
##    license for STK Astrogator. The default atmospheric model used will    ##
##    be the Jacchia-Roberts model. A full orbit propagation would be        ##
##    performed, with STK GUI also being launched. Processing time will be   ##
##    a lot slower than in the native Sams mode.                             ##
##                                                                           ##
##    This code is quite difficult to understand even with my comments.      ##
##    If the reader wants to understand the code, it is best to first go     ##
##    through the STK Object Model API, documentation, and tutorial.         ##
##                                                                           ##
##    Important note, make sure you close any running instances of STK.      ##
##    (i.e. check your task manager)                                         ##
##                                                                           ##
##    Written by Samuel Y. W. Low.                                           ##
##    First created 12-Oct-2020 10:18 AM (+8 GMT)                            ##
##    Last modified 19-Sep-2021 22:27 PM (-7 GMT)                            ##
##                                                                           ##
###############################################################################
###############################################################################

# Import basic utilities
import os
import datetime
import comtypes
import numpy as np

# Import libraries needed to interact with COM
from comtypes.client import CreateObject
from comtypes.client import GetActiveObject

# To the User - for changes to program inputs, it is recommended to change the
# parameters in the ORBITM GUI, rather than through the `config.txt` file or
# even through the code below.

def orbmstk(orbm_mode, tstart, tfinal, sc_Cd, sc_area_d,
            orb_a, orb_e, orb_i, orb_R, orb_w, orb_m,
            maintenance_tolerance, maintenance_fro,
            sc_mass, isp_min, isp_max):
    '''This function runs a full orbit maintenance scenario using STK
    Astrogator (either v10 or v11). Note that the user must have a valid
    license for STK Astrogator. The default atmospheric model used will be
    the Jacchia-Roberts model. A full orbit propagation would be performed,
    with STK GUI also being launched. Processing time will be a lot slower
    than in the native Sams mode.
    
    Parameters
    ----------
    orbm_mode : int
        Orbit simulation program choice (2 for STK10, 3 for STK11)
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
    
    # The parameters below are not used, but will be set to defaults.
    thr_TankPressure    = 0.1 # tank pressure (Pa)
    thr_TankVolume      = 0.1 # tank volume (m^3)
    thr_FuelDensity     = 0.1 # fuel density (kg/m^3)
    thr_FuelMass        = 0.1 # fuel mass (kg)
    thr_MaximumFuelMass = 1.0 # max fuel mass (kg)
    
    # Check below if you want the STK GUI to open up too (default True)
    stk_gui = True
    
    # User, check if you are using STK 10 or 11. By default, the comtypes GUID
    # is using the STK10 GUID code that allows the Astrogator wrapper to load.
    # Without successful loading, the AgStkGatorLib file cannot be recognised
    # and created in the comtypes gen folder, and AstroGator cannot be run.
    
    # GUID for STK10: 90E096F9-9615-4BA8-BA23-680F8D236959
    GUID_STK10 = "{90E096F9-9615-4BA8-BA23-680F8D236959}"
    # GUID for STK11: 090D317C-31A7-4AF7-89CD-25FE18F4017C
    GUID_STK11 = "{090D317C-31A7-4AF7-89CD-25FE18F4017C}"
    
    # Replace below where necessary.
    if orbm_mode == 2:
        comtypes.client.GetModule((comtypes.GUID( GUID_STK10 ),1,0))
    elif orbm_mode == 3:
        comtypes.client.GetModule((comtypes.GUID( GUID_STK11 ),1,0))
    
    # As a rule of thumb, frozen repeat orbit maintenance generally takes
    # 02x as much Delta-V per thrust as regular altitude maintenance due to
    # the need for the thrusts to bring the SC above the reference to maintain
    # the eastward-to-westward ground track shift. However, the thrust 
    # frequency would only be half as frequent.
    
    ##########################################################################
    ##########################################################################
    ##                                                                      ##
    ##    TO THE USER: DO NOT CHANGE ANY OF THE CODE BELOW, AS THE CODE     ##
    ##    IS HIGHLY DEPENDENT ON INTERFACING WITH THE RIGHT POINTERS TO     ##
    ##        THE RIGHT CLASSES. EDIT THE CODE BELOW, AT YOUR RISK!         ##
    ##                                                                      ##
    ##########################################################################
    ##########################################################################
    
    # The program will now compute the total scenario time in seconds.
    months_dict = {'Jan':1, 'Feb':2, 'Mar':3, 'Apr':4,
                   'May':5, 'Jun':6, 'Jul':7, 'Aug':8,
                   'Sep':9, 'Oct':10,'Nov':11,'Dec':12}
    
    # Read the start epoch string as a datetime object
    tstart_dt = datetime.datetime(int(tstart[6:10]),
                                  int(months_dict[tstart[2:5]]),
                                  int(tstart[0]),
                                  int(tstart[11:13]),
                                  int(tstart[14:16]),
                                  int(tstart[17:19]))
    
    # Read the final epoch string as a datetime object
    tfinal_dt = datetime.datetime(int(tfinal[6:10]),
                                  int(months_dict[tfinal[2:5]]),
                                  int(tfinal[0]),
                                  int(tfinal[11:13]),
                                  int(tfinal[14:16]),
                                  int(tfinal[17:19]))
    
    # Read the duration between start and final as a datetime.timedelta object.
    tdelta_dt = tfinal_dt - tstart_dt
    
    # Compute the total scenario time in seconds
    tdelta = (tdelta_dt.days*86400) + tdelta_dt.seconds # int
    
    # The program will now compute what is the desired Delta-V per thrust
    # using a first order Taylor expansion of the Vis-Visa equation.
    
    # GM = 398.6004415e12 # gravity constant x Earth mass (m**3/s**2)
    velocity = ( ( 398.6004415e12 ) / ( orb_a * 1000 ) )**0.5
    delta_v = (0.25 * velocity * maintenance_tolerance) / (orb_a) # m/s
    
    # First, try to close any existing STK applications.
    print("Closing any pre-existing STK applications... \n")
    print("Check if you need to save existing scenarios? (Open the UI) \n")
    
    # Check if the user is running in STK 10
    if orbm_mode == 2:
        try:
            uiApp = GetActiveObject('STK10.Application')
            uiApp.Quit()
        except:
            print("Wrong STK version selected! Currently selected STK 10. \n")
            pass
    
    # Check if the user is running in STK 11
    elif orbm_mode == 3:
        try:
            uiApp = GetActiveObject('STK11.Application')
            uiApp.Quit()
        except:
            print("Wrong STK version selected! Currently selected STK 11. \n")
            pass
    
    ##########################################################################
    ##########################################################################
    ##                                                                      ##
    ##             THIS SECTION INITIALISES THE STK APPLICATION             ##
    ##                                                                      ##
    ##########################################################################
    ##########################################################################
    
    # Start STK10 Application
    print("Creating a new STK application. \n")
    if orbm_mode == 2:
        uiApp = CreateObject("STK10.Application")
    elif orbm_mode == 3:
        uiApp = CreateObject("STK11.Application")
    uiApp.Visible = stk_gui
    uiApp.UserControl = stk_gui
    stkRoot = uiApp.Personality2
    
    from comtypes.gen import STKObjects
    from comtypes.gen import STKUtil
    from comtypes.gen import AgStkGatorLib
    from comtypes.client import gen_dir
    
    print("Creating the STK scenario object. \n")
    stkRoot.NewScenario("Orbit_Maintenance")
    
    # Get a reference to the scenario object (null if no scenario loaded)
    scenario = stkRoot.CurrentScenario
    scenario2 = scenario.QueryInterface(STKObjects.IAgScenario)
    
    # Set the time period for the scenario.
    try:
        scenario2.StartTime = tstart
        scenario2.StopTime = tfinal
    except:
        scenario2.StopTime = tfinal
        scenario2.StartTime = tstart
    
    #Reset STK to the new start time
    stkRoot.Rewind()
    
    ##########################################################################
    ##########################################################################
    ##                                                                      ##
    ##    THIS SECTION CREATES A DUMMY SATELLITE FOR LIFE TIME ESTIMATION   ##
    ##                                                                      ##
    ##########################################################################
    ##########################################################################
    
    # This segment will create the life-time test satellite and propagate it.
    print("Creating the satellite (life-time) object. \n")
    sat = scenario.Children.New(STKObjects.eSatellite, 'Lifetime')
    sat2 = sat.QueryInterface(STKObjects.IAgSatellite)
    
    # You can gain access to the initial orbit state through the satellite's
    # propagator object. In the block below, get a pointer to the interface
    # IAgVePropagtorTwoBody. Then use that pointer to convert the orbit 
    # state into the classical representation, and obtain a pointer to the 
    # IAgOrbitStateClassical interface.
    
    sat2.SetPropagatorType(STKObjects.ePropagatorTwoBody)
    sat2prop  = sat2.Propagator
    sat2prop2 = sat2prop.QueryInterface(STKObjects.IAgVePropagatorTwoBody)
    sat2init = sat2prop2.InitialState.Representation
    sat2state = sat2init.ConvertTo(STKUtil.eOrbitStateClassical)
    sat2state2 = sat2state.QueryInterface(STKObjects.IAgOrbitStateClassical)
    
    # With the IAgOrbitStateClassical interface you will be able to set 
    # the values of the desired orbital elements.
    
    # The SizeShape property only provides a pointer to the interface 
    # IAgClassicalSizeShape. It does not immediately provide access to the 
    # semimajor axis or eccentricity values. To access those, you "cast" 
    # to the interface IAgClassicalSizeShapeSemimajorAxis provided by the 
    # object AgClassicalSizeShapeSemimajorAxis. 
    
    # Initialise all the pointers.
    pointer_sma = STKObjects.IAgClassicalSizeShapeSemimajorAxis # Same as ECC
    pointer_ecc = STKObjects.IAgClassicalSizeShapeSemimajorAxis # Same as SMA
    pointer_ran = STKObjects.IAgOrientationAscNodeRAAN
    pointer_muu = STKObjects.IAgClassicalLocationMeanAnomaly
    
    # Set the semi-major axis and eccentricity
    sat2state2.SizeShapeType = STKObjects.eSizeShapeSemimajorAxis
    sat2state2.SizeShape.QueryInterface( pointer_sma ).SemiMajorAxis = orb_a
    sat2state2.SizeShape.QueryInterface( pointer_ecc ).Eccentricity = orb_e
    
    # Set the inclination and argument of perigee
    sat2state2.Orientation.Inclination = orb_i
    sat2state2.Orientation.ArgOfPerigee = orb_w
    
    # For the RAAN, as in the case of the semi-major axis and eccentricity,
    # you must first specify the AscNodeType, then provide the value for 
    # the AscNode through the approriate interface.
    sat2state2.Orientation.AscNodeType = STKObjects.eAscNodeRAAN
    sat2state2.Orientation.AscNode.QueryInterface( pointer_ran ).Value = orb_R
    
    # Set the mean anomaly
    sat2state2.LocationType = STKObjects.eLocationMeanAnomaly
    sat2state2.Location.QueryInterface( pointer_muu ).Value = orb_m
    
    # Propagate the orbit
    sat2prop2.InitialState.Representation.Assign(sat2state2)
    sat2prop2.Propagate()
    
    # Prepare the STK Connect Command strings for the life-time computation.
    setLifeTime              = 'SetLifetime */Satellite/Lifetime '
    setLifeTimeDragCoeff     = setLifeTime + 'DragCoeff ' + str(sc_Cd)
    setLifeTimeReflectCoeff  = setLifeTime + 'ReflectCoeff ' + str(0.000001)
    setLifeTimeDragArea      = setLifeTime + 'DragArea ' + str(sc_area_d)
    setLifeTimeSunArea       = setLifeTime + 'SunArea ' + str(0.000001)
    setLifeTimeMass          = setLifeTime + 'Mass ' + str(sc_mass)
    setLifeTimeLimitType     = setLifeTime + 'LimitType Duration'
    setLifeTimeDurationLimit = setLifeTime + 'DurationLimit 3650'
    setLifeTimeDensityModel  = setLifeTime + 'DensityModel Jacchia70Lifetime'
    
    # Execute the STK Connect Command strings for life-time computation.
    stkRoot.ExecuteCommand( setLifeTimeDragCoeff     )
    stkRoot.ExecuteCommand( setLifeTimeReflectCoeff  )
    stkRoot.ExecuteCommand( setLifeTimeDragArea      )
    stkRoot.ExecuteCommand( setLifeTimeSunArea       )
    stkRoot.ExecuteCommand( setLifeTimeMass          )
    stkRoot.ExecuteCommand( setLifeTimeLimitType     )
    stkRoot.ExecuteCommand( setLifeTimeDurationLimit )
    stkRoot.ExecuteCommand( setLifeTimeDensityModel  )
    
    # Execute the STK Connect Command strings for life-time computation.
    resultsLifeTime = stkRoot.ExecuteCommand('Lifetime */Satellite/Lifetime')
    lifetime_str = resultsLifeTime.Item(0) + " \n"
    lifetime_str.replace('.The', '. The')
    print(lifetime_str)
    
    # Finally, remove the test satellite used to compute the life time.
    sat.Unload()
    
    ##########################################################################
    ##########################################################################
    ##                                                                      ##
    ##   THIS SECTION CREATES THE PRIMARY SATELLITE FOR ORBIT MAINTENANCE   ##
    ##   ASTROGATOR MISSION CONTROL SEQUENCE (MCS) WITH ALL PROPAGATE AND   ##
    ##    MANEUVER SEQUENCES WILL ALSO BE SET UP IN THE CODE BLOCK BELOW    ##
    ##                                                                      ##
    ##########################################################################
    ##########################################################################
    
    # This segment will create the test satellite and propagate it.
    print("Creating the satellite object with orbit maintenance. \n")
    satellite = scenario.Children.New(STKObjects.eSatellite, "Satellite")
    satellite2 = satellite.QueryInterface(STKObjects.IAgSatellite)
    satellite2.SetPropagatorType(STKObjects.ePropagatorAstrogator)
    
    # For AstroGator, we need to access a special class called IAgVADriverMCS.
    # Acquire an interface to the DriverMCS interface of Astrogator through 
    # the propagator object. This is the central interface from which to 
    # control the satellite via Astrogator.
    
    print("Creating the MCS interface object to Astrogator. \n")
    astg = satellite2.Propagator.QueryInterface(AgStkGatorLib.IAgVADriverMCS)
    mcs  = astg.MainSequence
    mcs.RemoveAll() # Clear all sequences
    
    # Next, we set the initial states of the satellite.
    # The respective arguments are the segment type, name of segment, and the
    # name of the segment where the segment of interest is inserted before.
    mcs.Insert(AgStkGatorLib.eVASegmentStateInitial,'Initial State','-')
    
    # Get the initial state and query its interface
    mcs_init  = mcs.Item('Initial State')
    mcs_init2 = mcs_init.QueryInterface(AgStkGatorLib.IAgVAMCSInitialState)
    
    # Set the orbit elements, get the elements and query its interface.
    mcs_init2.SetElementType(1) # Keplerian
    mcs_elem  = mcs_init2.Element
    mcs_init2.OrbitEpoch = tstart
    mcs_elem2 = mcs_elem.QueryInterface(AgStkGatorLib.IAgVAElementKeplerian)
    
    # Set the orbit elements
    print("Creating and setting the orbit elements. \n")
    mcs_elem2.ArgOfPeriapsis = orb_w
    mcs_elem2.Eccentricity   = orb_e
    mcs_elem2.Inclination    = orb_i
    mcs_elem2.RAAN           = orb_R
    mcs_elem2.SemiMajorAxis  = orb_a
    mcs_elem2.TrueAnomaly    = orb_m
    
    # Query the interface that allows setting of the spacecraft parameters.
    print("Creating and setting the spacecraft parameters. \n")
    pointer_mcs_scparams = AgStkGatorLib.IAgVASpacecraftParameters
    mcs_scparams = mcs_init2.SpacecraftParameters
    mcs_scparams2 = mcs_scparams.QueryInterface( pointer_mcs_scparams )
    
    # Set the spacecraft parameters
    mcs_scparams2.Cd                         = sc_Cd
    mcs_scparams2.Ck                         = 0.000001
    mcs_scparams2.Cr                         = 0.000001
    mcs_scparams2.DryMass                    = sc_mass
    mcs_scparams2.DragArea                   = sc_area_d
    mcs_scparams2.RadiationPressureArea      = 0.000001
    mcs_scparams2.SolarRadiationPressureArea = 0.000001
    
    print("Creating and setting spacecraft fuel tank parameters. \n")
    # Query the interface that allows setting of the fuel tank parameters.
    mcs_fueltank = mcs_init2.FuelTank
    mcs_fueltank2 = mcs_fueltank.QueryInterface(AgStkGatorLib.IAgVAFuelTank)
    
    # Set the fuel tank parameters
    mcs_fueltank2.TankPressure    = thr_TankPressure
    mcs_fueltank2.TankVolume      = thr_TankVolume
    mcs_fueltank2.FuelDensity     = thr_FuelDensity
    mcs_fueltank2.FuelMass        = thr_FuelMass
    mcs_fueltank2.MaximumFuelMass = thr_MaximumFuelMass
    
    # Now, we set the automatic sequence conditions for station-keeping.
    print("Creating the Automatic Sequence object. \n")
    acs = astg.AutoSequence
    acs.Add("Maintain")
    acs_main = acs.Item("Maintain")
    acs_main2 = acs_main.QueryInterface(AgStkGatorLib.IAgVAAutomaticSequence)
    acs_seq = acs_main2.Sequence
    
    # Let us initialise some pointers to the propagate and the maneuver
    # sequences used in STK AstroGator.
    pointer_seg_prop = AgStkGatorLib.eVASegmentTypePropagate
    pointer_seg_manv = AgStkGatorLib.eVASegmentTypeManeuver
    pointer_mcs_prop = AgStkGatorLib.IAgVAMCSPropagate
    pointer_mcs_manv = AgStkGatorLib.IAgVAMCSManeuver
    pointer_thr_vec  = AgStkGatorLib.eVAAttitudeControlThrustVector
    pointer_thr_imp  = AgStkGatorLib.IAgVAAttitudeControlImpulsiveThrustVector
    pointer_stopcond = AgStkGatorLib.IAgVAStoppingCondition
    
    # In the ACS, we add the propagate segment, change the propagate segment
    # to a UserSelect option, change the sequence to 'Maintain', and the
    # stopping condition of the UserSelect option should use a UserCalcObject
    # that is the Kozai-Iszak Mean SMA of the orbit. If the LEO crosses the
    # mean SMA threshold in the MCS, it will prompt the trigger of the ACS.
    
    # If the user is not maintaining an orbit for a frozen repeat, a single
    # thrust at the apogee should be good enough to raise its orbit.
    if maintenance_fro == False:
        
        # Begin inserting the propagation with an apogee thrust.
        Prop2Apo   = acs_seq.Insert(pointer_seg_prop,'Prop2Apo','-')
        ThrustApo  = acs_seq.Insert(pointer_seg_manv,'ThrustApo','-')
        
        # Now we query the interfaces for all of them.
        Prop2Apo2 = Prop2Apo.QueryInterface(pointer_mcs_prop)
        ThrustApo2 = ThrustApo.QueryInterface(pointer_mcs_manv)
        
        # We set the Prop2Apo segment to perform propagation to the apogee.
        Prop2Apo2_SC = Prop2Apo2.StoppingConditions
        Prop2Apo2_SC.Add('Apoapsis')
        Prop2Apo2_SC.Cut('Duration') # Not sure why remove doesn't work
        
        # Then, we set the thrust at the apogee using a fixed thrust.
        ThrustApo2Man = ThrustApo2.Maneuver
        ThrustApo2Man.SetAttitudeControlType(pointer_thr_vec)
        ThrustApo2ManAtt = ThrustApo2Man.AttitudeControl
        ThrustApo2_Vector = ThrustApo2ManAtt.QueryInterface(pointer_thr_imp)
        ThrustApo2_Vector.DeltaVVector.AssignCartesian(delta_v/1000, 0.0, 0.0)
    
    # If the user is maintaining an orbit for a frozen repeat, a second
    # thrust at perigee is needed to minimise changes to the eccentricity.
    if maintenance_fro == True:
        
        # Insert the propagation to perigee with the perigee thrust.
        Prop2Peri  = acs_seq.Insert(pointer_seg_prop,'Prop2Peri','-')
        ThrustPeri = acs_seq.Insert(pointer_seg_manv,'ThrustPeri','-')
        
        # Now we query the interfaces for all of them.
        Prop2Peri2 = Prop2Peri.QueryInterface(pointer_mcs_prop)
        ThrustPeri2 = ThrustPeri.QueryInterface(pointer_mcs_manv)
        
        # We set the Prop2Peri segment to perform propagation to the perigee.
        Prop2Peri2_SC = Prop2Peri2.StoppingConditions
        Prop2Peri2_SC.Add('Periapsis')
        Prop2Peri2_SC.Cut('Duration') # Not sure why remove doesn't work
        
        # Then, we set the thrust at the perigee using a fixed thrust.
        ThrustPeri2Man = ThrustPeri2.Maneuver
        ThrustPeri2Man.SetAttitudeControlType(pointer_thr_vec)
        ThrustPeri2ManAtt = ThrustPeri2Man.AttitudeControl
        ThrustPeri2_Vector = ThrustPeri2ManAtt.QueryInterface(pointer_thr_imp)
        ThrustPeri2_Vector.DeltaVVector.AssignCartesian(delta_v/1000, 0.0, 0.0)
        
        # Begin inserting the propagation with an apogee thrust.
        Prop2Apo   = acs_seq.Insert(pointer_seg_prop,'Prop2Apo','-')
        ThrustApo  = acs_seq.Insert(pointer_seg_manv,'ThrustApo','-')
        
        # Now we query the interfaces for all of them.
        Prop2Apo2 = Prop2Apo.QueryInterface(pointer_mcs_prop)
        ThrustApo2 = ThrustApo.QueryInterface(pointer_mcs_manv)
        
        # We set the Prop2Apo segment to perform propagation to the apogee.
        Prop2Apo2_SC = Prop2Apo2.StoppingConditions
        Prop2Apo2_SC.Add('Apoapsis')
        Prop2Apo2_SC.Cut('Duration') # Not sure why remove doesn't work
        
        # Then, we set the thrust at the apogee using a fixed thrust.
        ThrustApo2Man = ThrustApo2.Maneuver
        ThrustApo2Man.SetAttitudeControlType(pointer_thr_vec)
        ThrustApo2ManAtt = ThrustApo2Man.AttitudeControl
        ThrustApo2_Vector = ThrustApo2ManAtt.QueryInterface(pointer_thr_imp)
        ThrustApo2_Vector.DeltaVVector.AssignCartesian(delta_v/1000, 0.0, 0.0)
    
    # At this stage, the automatic sequence oject has been successfully built.
    # We just need to know how to call the automatic sequence whenever stop
    # conditions are met (i.e. when the satellite crosses the threshold)
    
    print("Setting the MCS segments and piecing everything together... \n")
    
    # For the MCS, it needs only a propagate segment with a duration = tdelta.
    PropMCS  = mcs.Insert(pointer_seg_prop,'PropMCS','-')
    PropMCS2 = PropMCS.QueryInterface(pointer_mcs_prop)
    PropMCS2_Stop_Properties = PropMCS2.StoppingConditions.Item(0).Properties
    PropMCS2_Stop_Properties.QueryInterface(pointer_stopcond).Trip = tdelta
    
    # We will add the mean semi-major axis as a stopping condition.
    PropMCS2_Stop_SMA = PropMCS2.StoppingConditions.Add('UserSelect')
    PropMCS2_Stop_SMA_P  = PropMCS2_Stop_SMA.Properties
    PropMCS2_Stop_SMA_P2 = PropMCS2_Stop_SMA_P.QueryInterface(pointer_stopcond)
    PropMCS2_Stop_SMA_P2.UserCalcObjectName = 'Mean Semimajor Axis'
    PropMCS2_Stop_SMA_P2.Trip = orb_a - maintenance_tolerance
    PropMCS2_Stop_SMA_P2.RepeatCount = 1
    PropMCS2_Stop_SMA_P2.Sequence = 'Maintain'
    
    # AGIers say: You can set multiple stopping conditions for a propagate
    # segment. Astrogator stops propagating the satellite when it meets one 
    # of the stopping conditions.
    
    # Run the MCS
    print("Running mission control sequence (this might take long)... \n")
    astg.RunMCS()
    
    print("Mission successfully ran! Now extracting orbital data. \n")
    
    ##########################################################################
    ##########################################################################
    ##                                                                      ##
    ##   THIS SECTION EXTRACTS THE DATA PROVIDERS FROM THE PROPAGATED SAT   ##
    ##    (NOTE - DATA PROVIDERS ARE SET TO EPOCH-SECONDS FOR STABILITY)    ##
    ##                                                                      ##
    ##########################################################################
    ##########################################################################
    
    # Now, we need to extract relevant data. The data we will need are:
    sat_epochs = []   # Data Provider 0: Epoch seconds of STK scenario
    sat_mean_sma = [] # Data Provider 1: Mean semimajor axes values
    sat_altitude = [] # Data Provider 2: Altitude values
    
    # Change DateFormat dimension to epoch seconds for ease of data handling
    stkRoot.UnitPreferences.Item('DateFormat').SetCurrentUnit('EpSec')
    
    # Get the pointer to Data Provider Groups
    pointer_data = STKObjects.IAgDataProviderGroup
    pointer_tvar = STKObjects.IAgDataPrvTimeVar
    
    # Get the Kozai-Izsak Mean data providers pointer and interface to it.
    SatMean  = satellite.DataProviders.Item("Kozai-Izsak Mean")
    SatMean2 = SatMean.QueryInterface( pointer_data )
    SatMean2_ICRF  = SatMean2.Group.Item("ICRF")
    SatMean2_ICRF2 = SatMean2_ICRF.QueryInterface( pointer_tvar )
    SatMean2_ICRF2_Get = SatMean2_ICRF2.Exec(scenario2.StartTime,
                                              scenario2.StopTime,
                                              3600)
    SatMean2_ICRF2_SMA = np.array(SatMean2_ICRF2_Get.DataSets.ToArray())
    
    # Get altitude from the LLA State-Fixed pointer and interface to it.
    SatAlt  = satellite.DataProviders.Item("LLA State")
    SatAlt2 = SatAlt.QueryInterface( pointer_data )
    SatAlt2_Fixed  = SatAlt2.Group.Item("Fixed")
    SatAlt2_Fixed2 = SatAlt2_Fixed.QueryInterface( pointer_tvar )
    SatAlt2_Fixed2_Get = SatAlt2_Fixed2.Exec(scenario2.StartTime,
                                              scenario2.StopTime,
                                              3600)
    SatAlt2_Fixed2_Alt = np.array(SatAlt2_Fixed2_Get.DataSets.ToArray())
    
    # Check if the epochs of two data providers are not the same (debug)
    if len(SatAlt2_Fixed2_Alt) != len(SatMean2_ICRF2_SMA):
        print("Warning! Something went wrong with the data provider parsing.")
        print("Length of ALT and SMA arrays do not match! Code broken? \n")
    
    # Extract the mean semi-major axes and altitude values
    for epoch in range(0,len(SatMean2_ICRF2_SMA)):
        sat_epochs.append(float(SatMean2_ICRF2_SMA[epoch][0]))
        sat_mean_sma.append(float(SatMean2_ICRF2_SMA[epoch][1]))
        sat_altitude.append(float(SatAlt2_Fixed2_Alt[epoch][3]))
    
    # Get the maneuver summary data providers pointer and interface to it.
    sat_deltaV  = satellite.DataProviders.Item("Maneuver Summary")
    sat_deltaV2 = sat_deltaV.QueryInterface(STKObjects.IAgDataPrvInterval)
    sat_deltaV2_Data = sat_deltaV2.Exec(scenario2.StartTime,
                                        scenario2.StopTime)
    sat_deltaV2_Array = np.array(sat_deltaV2_Data.DataSets.ToArray())
    
    # Extract the Delta-V values into a text file.
    deltaV_file = open("output_manoeuvres.txt",'w')
    for thrust in sat_deltaV2_Array:
        thrust_str  = thrust[0] + ' '
        thrust_str += thrust[1] + ' '
        thrust_str += thrust[2] + ' '
        thrust_str += thrust[5] + ' \n'
        deltaV_file.write(thrust_str)
    deltaV_file.close()
    
    # Compute the total Delta-V (m/s) and the total impulse needed.
    total_DV = len(sat_deltaV2_Array) * delta_v
    
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
