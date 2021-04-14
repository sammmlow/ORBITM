# -*- coding: utf-8 -*-

###############################################################################
###############################################################################
##                                                                           ##
##     _____ ___  ____  ___  _____       ______                              ##
##    |  _  | _ \|  _ \|_ _||_   _|     |      |                             ##
##    | |_| |   <|  _ < | |   | |       | \  / |  _                          ##
##    |_____|_|\_|____/|___|  |_|       |_|\/|_| |_|                         ##
##                                                     v 1.0                 ##
##                                                                           ##
##    FILE DESCRIPTION:                                                      ##
##                                                                           ##
##    This code is quite difficult to understand even with my comments.      ##
##    If the reader wants to understand the code, it is best to first go     ##
##    through the STK Object Model API, documentation, and tutorial.         ##
##                                                                           ##
##    Important note, make sure you close any running instances of STK.      ##
##    (i.e. check your task manager)                                         ##
##                                                                           ##
##    Written by Samuel Y. W. Low.                                           ##
##    First created 12-10-2020 10:18 AM (+8 GMT)                             ##
##    Last modified 30-03-2021 08:33 PM (+8 GMT)                             ##
##                                                                           ##
###############################################################################
###############################################################################

# Import basic utilities
import os
import datetime
import comtypes
import numpy as np
import matplotlib.pyplot as plt

# Needed to interact with COM
from comtypes.client import CreateObject
from comtypes.client import GetActiveObject

def orbm_run_stk(orbm_mode, tstart, tfinal,
                 sc_Cd, sc_area_d, sc_Ck, sc_area_a, sc_Cr, sc_area_r,
                 orb_a, orb_e, orb_i, orb_R, orb_w, orb_m,
                 maintenance_tolerance,
                 maintenance_margin,
                 maintenance_fro,
                 sc_mass, isp_min, isp_max):

    # The parameters below are not used, but will be set to defaults.
    thr_TankPressure    = 0.1 # tank pressure (Pa)
    thr_TankVolume      = 0.1 # tank volume (m^3)
    thr_FuelDensity     = 0.1 # fuel density (kg/m^3)
    thr_FuelMass        = 0.1 # fuel mass (kg)
    thr_MaximumFuelMass = 1.0 # max fuel mass (kg)
    
    # For thruster sizing and plotting, what range of Isp is needed?
    plot_Isp_Min = isp_min # s
    plot_Isp_Max = isp_max # s
    
    # Check below if you want the STK GUI to open up too (default True)
    stk_gui = True
    
    # User, check if you are using STK10 or STK11. By default, the comtypes GUID
    # is using the STK10 GUID code that allows the Astrogator wrapper to load.
    # Without successful loading, the AgStkGatorLib file cannot be recognised
    # and created in the comtypes gen folder, and AstroGator cannot be run.
    
    # GUID for STK10: 90E096F9-9615-4BA8-BA23-680F8D236959
    # GUID for STK11: 090D317C-31A7-4AF7-89CD-25FE18F4017C
    
    # Replace below where necessary.
    if orbm_mode == 2:
        comtypes.client.GetModule((comtypes.GUID("{90E096F9-9615-4BA8-BA23-680F8D236959}"),1,0))
    elif orbm_mode == 3:
        comtypes.client.GetModule((comtypes.GUID("{090D317C-31A7-4AF7-89CD-25FE18F4017C}"),1,0))
    
    # As a rule of thumb, frozen repeat orbit maintenance generally takes about
    # 02x as much Delta-V per thrust as regular altitude maintenance due to the 
    # need for the thrusts to bring the SC above the reference to maintain the 
    # eastward-to-westward ground track shift.
    
    """ #######################################################################
    
    TO THE USER: DO NOT CHANGE ANY OF THE CODE BELOW, AS THE CODE IS HIGHLY
    DEPENDENT ON INTERFACING WITH THE RIGHT POINTERS TO THE RIGHT CLASSES.
    EDIT THE CODE BELOW, AT YOUR RISK, AND ONLY IF YOU KNOW WHAT YOU ARE DOING!
    
    ####################################################################### """
    
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
    
    # Read the time delta between start and final as a datetime-timedelta object
    tdelta_dt = tfinal_dt - tstart_dt
    
    # Compute the total scenario time in seconds
    tdelta = (tdelta_dt.days*86400) + tdelta_dt.seconds # int
    
    ############################################################################
    ############################################################################
    
    # The program will now compute what is the desired Delta-V per thrust
    # using a first order Taylor expansion of the Vis-Visa equation.
    
    GM = 398.6004415e12 # gravity constant x Earth mass (m**3/s**2)
    velocity = ((398.6004415e12)/(orb_a*1000))**0.5
    delta_v = (0.25*velocity*maintenance_tolerance)/(orb_a * 1000) # km/s
    
    ############################################################################
    ############################################################################
    
    # First, try to close any existing STK applications.
    print("Closing any pre-existing STK applications... \n")
    print("Check if you need to save your existing scenarios? (Open the UI) \n")
    # Check if the user is running in STK 10
    if orbm_mode == 2:
        try:
            uiApp = GetActiveObject('STK10.Application')
            uiApp.Quit()
        except:
            pass
    
    # Check if the user is running in STK 11
    elif orbm_mode == 3:
        try:
            uiApp = GetActiveObject('STK11.Application')
            uiApp.Quit()
        except:
            pass
    
    ############################################################################
    ############################################################################
    
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
    scenario2.SetTimePeriod( tstart, tfinal )
    
    #Reset STK to the new start time
    stkRoot.Rewind()
    
    ############################################################################
    ############################################################################
    
    # This segment will create the life-time test satellite and propagate it.
    print("Creating the satellite (life-time) object. \n")
    sat = scenario.Children.New(STKObjects.eSatellite, 'Lifetime')
    sat2 = sat.QueryInterface(STKObjects.IAgSatellite)
    
    # You can gain access to the initial orbit state through the satellite's
    # propagator object. In the block below, get a pointer to the interface
    # IAgVePropagtorTwoBody. Then use that pointer to convert the orbit state
    # into the classical representation, and obtain a pointer to the interface
    # IAgOrbitStateClassical.
    sat2.SetPropagatorType(STKObjects.ePropagatorTwoBody)
    sat2prop = sat2.Propagator.QueryInterface(STKObjects.IAgVePropagatorTwoBody)
    sat2init = sat2prop.InitialState.Representation
    sat2state = sat2init.ConvertTo(STKUtil.eOrbitStateClassical)
    sat2state2 = sat2state.QueryInterface(STKObjects.IAgOrbitStateClassical)
    
    # With the IAgOrbitStateClassical interface you will be able to set the values
    # of the desired orbital elements.
    
    # The SizeShape property only provides a pointer to the IAgClassicalSizeShape
    # interface, which does not immediately provide access to the semimajor axis
    # or eccentricity values. To access those, you "cast" to the interface
    # IAgClassicalSizeShapeSemimajorAxis provided by the object
    # AgClassicalSizeShapeSemimajorAxis. 
    sat2state2.SizeShapeType = STKObjects.eSizeShapeSemimajorAxis
    sat2state2.SizeShape.QueryInterface(STKObjects.IAgClassicalSizeShapeSemimajorAxis).SemiMajorAxis = orb_a
    sat2state2.SizeShape.QueryInterface(STKObjects.IAgClassicalSizeShapeSemimajorAxis).Eccentricity = orb_e
    
    # Set the inclination and argument of perigee
    sat2state2.Orientation.Inclination = orb_i
    sat2state2.Orientation.ArgOfPerigee = orb_w
    
    # For the RAAN, much as in the case of the semi-major axis and eccentricity,
    # you must first specify the AscNodeType, then provide the value for the
    # AscNode through the approriate interface.
    sat2state2.Orientation.AscNodeType = STKObjects.eAscNodeRAAN
    sat2state2.Orientation.AscNode.QueryInterface(STKObjects.IAgOrientationAscNodeRAAN).Value = orb_R
    
    # Set the mean anomaly
    sat2state2.LocationType = STKObjects.eLocationMeanAnomaly
    sat2state2.Location.QueryInterface(STKObjects.IAgClassicalLocationMeanAnomaly).Value = orb_m
    
    # Propagate the orbit
    sat2prop.InitialState.Representation.Assign(sat2state2)
    sat2prop.Propagate()
    
    # Prepare the STK Connect Command strings for the life-time computation.
    setLifeTime              = 'SetLifetime */Satellite/Lifetime '
    setLifeTimeDragCoeff     = setLifeTime + 'DragCoeff ' + str(sc_Cd)
    setLifeTimeReflectCoeff  = setLifeTime + 'ReflectCoeff ' + str(sc_Cr)
    setLifeTimeDragArea      = setLifeTime + 'DragArea ' + str(sc_area_d)
    setLifeTimeSunArea       = setLifeTime + 'SunArea ' + str(sc_area_r)
    setLifeTimeMass          = setLifeTime + 'Mass ' + str(sc_mass)
    setLifeTimeLimitType     = setLifeTime + 'LimitType Duration'
    setLifeTimeDurationLimit = setLifeTime + 'DurationLimit 3650'
    setLifeTimeDensityModel  = setLifeTime + 'DensityModel Jacchia70Lifetime'
    
    # Execute the STK Connect Command strings for life-time computation settings.
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
    print(lifetime_str)
    
    # Finally, remove the test satellite used to compute the life time.
    sat.Unload()
    
    ############################################################################
    ############################################################################
    
    # This segment will create the test satellite and propagate it.
    print("Creating the satellite object with orbit maintenance. \n")
    satellite = scenario.Children.New(STKObjects.eSatellite, "Satellite")
    #print("Querying the IAgStkObject interface of the satellite. \n")
    satellite2 = satellite.QueryInterface(STKObjects.IAgSatellite)
    satellite2.SetPropagatorType(STKObjects.ePropagatorAstrogator) # Astrogator
    
    # For AstroGator, we need to access a special class called the IAgVADriverMCS.
    # Acquire an interface to the DriverMCS interface of Astrogator through the
    # propagator object. This is the central interface from which to control the
    # satellite via Astrogator.
    
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
    
    # Set the orbit elements, get the elements and query its interface
    mcs_init2.SetElementType(1) # Keplerian
    mcs_elem  = mcs_init2.Element
    mcs_init2.OrbitEpoch = tstart
    mcs_elem2 = mcs_elem.QueryInterface(AgStkGatorLib.IAgVAElementKeplerian)
    
    print("Creating and setting the orbit elements. \n")
    # Set the orbit elements
    mcs_elem2.ArgOfPeriapsis = orb_w
    mcs_elem2.Eccentricity   = orb_e
    mcs_elem2.Inclination    = orb_i
    mcs_elem2.RAAN           = orb_R
    mcs_elem2.SemiMajorAxis  = orb_a
    mcs_elem2.TrueAnomaly    = orb_m
    
    print("Creating and setting the spacecraft parameters. \n")
    # Query the interface that allows setting of the spacecraft parameters
    mcs_scparams = mcs_init2.SpacecraftParameters
    mcs_scparams2 = mcs_scparams.QueryInterface(AgStkGatorLib.IAgVASpacecraftParameters)
    
    # Set the spacecraft parameters
    mcs_scparams2.Cd                         = sc_Cd
    mcs_scparams2.Ck                         = sc_Ck
    mcs_scparams2.Cr                         = sc_Cr
    mcs_scparams2.DryMass                    = sc_mass
    mcs_scparams2.DragArea                   = sc_area_d
    mcs_scparams2.RadiationPressureArea      = sc_area_r
    mcs_scparams2.SolarRadiationPressureArea = sc_area_a
    
    print("Creating and setting spacecraft fuel tank parameters. \n")
    # Query the interface that allows setting of the fuel tank parameters
    mcs_fueltank = mcs_init2.FuelTank
    mcs_fueltank2 = mcs_fueltank.QueryInterface(AgStkGatorLib.IAgVAFuelTank)
    
    # Set the fuel tank parameters
    mcs_fueltank2.TankPressure    = thr_TankPressure
    mcs_fueltank2.TankVolume      = thr_TankVolume
    mcs_fueltank2.FuelDensity     = thr_FuelDensity
    mcs_fueltank2.FuelMass        = thr_FuelMass
    mcs_fueltank2.MaximumFuelMass = thr_MaximumFuelMass
    
    # Now we are going to set the automatic sequence conditions for station-keeping
    print("Creating the Automatic Sequence object. \n")
    acs = astg.AutoSequence
    acs.Add("Maintain")
    acs_main = acs.Item("Maintain")
    acs_main2 = acs_main.QueryInterface(AgStkGatorLib.IAgVAAutomaticSequence)
    acs_seq = acs_main2.Sequence
    
    # In the ACS, we add the propagate segment, change the propagate segment to 
    # a UserSelect option, change the sequence to 'Maintain', and the stopping
    # condition of the UserSelect option should use a UserCalcObject that is the
    # Kozai-Iszak Mean SMA of the orbit. If the LEO crosses the mean SMA threshold
    # in the MCS, it will prompt the trigger of the ACS.
    
    # If the user is not maintaining an orbit for a frozen repeat, a single
    # thrust at the apogee should be good enough to raise its orbit.
    if maintenance_fro == False:
        
        # Begin inserting the propagation with an apogee thrust.
        Prop2Apo   = acs_seq.Insert(AgStkGatorLib.eVASegmentTypePropagate,'Prop2Apo','-')
        ThrustApo  = acs_seq.Insert(AgStkGatorLib.eVASegmentTypeManeuver,'ThrustApo','-')
        
        # Now we query the interfaces for all of them.
        Prop2Apo2 = Prop2Apo.QueryInterface(AgStkGatorLib.IAgVAMCSPropagate)
        ThrustApo2 = ThrustApo.QueryInterface(AgStkGatorLib.IAgVAMCSManeuver)
        
        # We set the Prop2Apo segment to perform an orbit propagation to the apogee.
        Prop2Apo2_SC = Prop2Apo2.StoppingConditions
        Prop2Apo2_SC.Add('Apoapsis')
        Prop2Apo2_SC.Cut('Duration') # Not sure why remove doesn't work
        
        # Then, we set the thruster fire levels at the apogee using a fixed thrust.
        ThrustApo2.Maneuver.SetAttitudeControlType(AgStkGatorLib.eVAAttitudeControlThrustVector)
        ThrustApo2_Vector = ThrustApo2.Maneuver.AttitudeControl.QueryInterface(AgStkGatorLib.IAgVAAttitudeControlImpulsiveThrustVector)
        ThrustApo2_Vector.DeltaVVector.AssignCartesian(delta_v, 0.0, 0.0);
    
    # If the user is maintaining an orbit for a frozen repeat, a second perigee
    # thrust is needed in order to minimise disrupting the osculating eccentricity.
    if maintenance_fro == True:
        
        # Insert the propagation to perigee with the perigee thrust.
        Prop2Peri  = acs_seq.Insert(AgStkGatorLib.eVASegmentTypePropagate,'Prop2Peri','-')
        ThrustPeri = acs_seq.Insert(AgStkGatorLib.eVASegmentTypeManeuver,'ThrustPeri','-')
        
        # Now we query the interfaces for all of them.
        Prop2Peri2 = Prop2Peri.QueryInterface(AgStkGatorLib.IAgVAMCSPropagate)
        ThrustPeri2 = ThrustPeri.QueryInterface(AgStkGatorLib.IAgVAMCSManeuver)
        
        # We set the Prop2Peri segment to perform an orbit propagation to the perigee.
        Prop2Peri2_SC = Prop2Peri2.StoppingConditions
        Prop2Peri2_SC.Add('Periapsis')
        Prop2Peri2_SC.Cut('Duration') # Not sure why remove doesn't work
        
        # Then, we set the thruster fire levels at the perigee using a fixed thrust.
        ThrustPeri2.Maneuver.SetAttitudeControlType(AgStkGatorLib.eVAAttitudeControlThrustVector)
        ThrustPeri2_Vector = ThrustPeri2.Maneuver.AttitudeControl.QueryInterface(AgStkGatorLib.IAgVAAttitudeControlImpulsiveThrustVector)
        ThrustPeri2_Vector.DeltaVVector.AssignCartesian(delta_v, 0.0, 0.0);
        
        # Begin inserting the propagation with an apogee thrust.
        Prop2Apo   = acs_seq.Insert(AgStkGatorLib.eVASegmentTypePropagate,'Prop2Apo','-')
        ThrustApo  = acs_seq.Insert(AgStkGatorLib.eVASegmentTypeManeuver,'ThrustApo','-')
        
        # Now we query the interfaces for all of them.
        Prop2Apo2 = Prop2Apo.QueryInterface(AgStkGatorLib.IAgVAMCSPropagate)
        ThrustApo2 = ThrustApo.QueryInterface(AgStkGatorLib.IAgVAMCSManeuver)
        
        # We set the Prop2Apo segment to perform an orbit propagation to the apogee.
        Prop2Apo2_SC = Prop2Apo2.StoppingConditions
        Prop2Apo2_SC.Add('Apoapsis')
        Prop2Apo2_SC.Cut('Duration') # Not sure why remove doesn't work
        
        # Then, we set the thruster fire levels at the apogee using a fixed thrust.
        ThrustApo2.Maneuver.SetAttitudeControlType(AgStkGatorLib.eVAAttitudeControlThrustVector)
        ThrustApo2_Vector = ThrustApo2.Maneuver.AttitudeControl.QueryInterface(AgStkGatorLib.IAgVAAttitudeControlImpulsiveThrustVector)
        ThrustApo2_Vector.DeltaVVector.AssignCartesian(delta_v, 0.0, 0.0);
    
    # At this stage, the automatic sequence oject has been successfully built.
    # We just need to know how to call the automatic sequence whenever the stop
    # conditions have been met (i.e. when the satellite crosses the threshold)
    
    print("Setting the MCS segments and piecing everything together... \n")
    
    # For the MCS, it needs only a propagate segment, with a duration = tdelta.
    PropMain  = mcs.Insert(AgStkGatorLib.eVASegmentTypePropagate,'PropMain','-')
    PropMain2 = PropMain.QueryInterface(AgStkGatorLib.IAgVAMCSPropagate)
    PropMain2_StopCon_Dura = PropMain2.StoppingConditions.Item(0)
    PropMain2_StopCon_Dura_Prop2 = PropMain2_StopCon_Dura.Properties.QueryInterface(AgStkGatorLib.IAgVAStoppingCondition)
    PropMain2_StopCon_Dura_Prop2.Trip = tdelta
    
    # We will add the mean semi-major axis as a stopping condition.
    PropMain2_StopCon_SMA = PropMain2.StoppingConditions.Add('UserSelect')
    PropMain2_StopCon_SMA_Properties   = PropMain2_StopCon_SMA.Properties
    PropMain2_StopCon_SMA_Properties2  = PropMain2_StopCon_SMA_Properties.QueryInterface(AgStkGatorLib.IAgVAStoppingCondition)
    PropMain2_StopCon_SMA_Properties2.UserCalcObjectName = 'Mean Semimajor Axis'
    PropMain2_StopCon_SMA_Properties2.Trip = orb_a - maintenance_tolerance
    PropMain2_StopCon_SMA_Properties2.RepeatCount = 1
    PropMain2_StopCon_SMA_Properties2.Sequence = 'Maintain'
    
    # AGIers say: You can set multiple stopping conditions for a propagate
    # segment. Astrogator stops propagating the satellite when it meets one of
    # the stopping conditions.
    
    # Run the MCS
    print("Running the mission control sequence now (this might take long)... \n")
    astg.RunMCS()
    
    print("Mission successfully ran! Now extracting orbital data. \n")
    # Now, we need to start extracting relevant data. The data we will need are:
    
    # 1 - Kozai-Iszak mean semimajor axes values.
    sat_epochs = []
    sat_mean_sma = []
    
    # 2 - Altitude values.
    sat_altitude = []
    
    # Get the Kozai-Izsak Mean data providers pointer and interface to it.
    sat_ki_mean  = satellite.DataProviders.Item("Kozai-Izsak Mean")
    sat_ki_mean2 = sat_ki_mean.QueryInterface(STKObjects.IAgDataProviderGroup)
    sat_ki_mean2_ICRF  = sat_ki_mean2.Group.Item("ICRF")
    sat_ki_mean2_ICRF2 = sat_ki_mean2_ICRF.QueryInterface(STKObjects.IAgDataPrvTimeVar)
    sat_ki_mean2_ICRF2_Data = sat_ki_mean2_ICRF2.Exec(scenario2.StartTime, scenario2.StopTime, 3600)
    sat_ki_mean_sma = np.array(sat_ki_mean2_ICRF2_Data.DataSets.ToArray())
    
    time_error_flag = False
    
    # Extract the Kozai-Izsak Mean semi major axis
    for epoch in range(0,len(sat_ki_mean_sma)):
        
        epochstr = sat_ki_mean_sma[epoch][0] # Read raw epoch string
        epochlis = epochstr.split()
        mean_sma = sat_ki_mean_sma[epoch][1] # Read raw mean SMA string
        
        yy = int(epochlis[2])
        mm = int(months_dict[epochlis[1]])
        dd = int(epochlis[0])
        hh = int(epochlis[3][0:2])
        mn = int(epochlis[3][3:5])
        ss = int(epochlis[3][6:8])
        
        # If STK outputs end-of-denominator numbers for some weird reason.
        if ss == 60:
            ss = 0
            time_error_flag = True
        
        epoch_dt = datetime.datetime(yy,mm,dd,hh,mn,ss)
        if time_error_flag == True:
            epoch_dt = epoch_dt + datetime.timedelta(seconds=60)
            time_error_flag = False
        
        sat_epochs.append(epoch_dt)
        sat_mean_sma.append(float(mean_sma))
    
    # Get the altitude values from the LLA State-Fixed pointer and interface to it.
    sat_alt  = satellite.DataProviders.Item("LLA State")
    sat_alt2 = sat_alt.QueryInterface(STKObjects.IAgDataProviderGroup)
    sat_alt2_Fixed  = sat_alt2.Group.Item("Fixed")
    sat_alt2_Fixed2 = sat_alt2_Fixed.QueryInterface(STKObjects.IAgDataPrvTimeVar)
    sat_alt2_Fixed2_Data = sat_alt2_Fixed2.Exec(scenario2.StartTime, scenario2.StopTime, 3600)
    sat_alt_data_final = np.array(sat_alt2_Fixed2_Data.DataSets.ToArray())
    
    if len(sat_alt_data_final) != len(sat_ki_mean_sma):
        print("Warning! Something went wrong with the data provider parsing.")
        print("Length of altitude and SMA arrays do not match! Code broken? \n")
    
    # Extract the altitude values
    for epoch in range(0,len(sat_alt_data_final)):
        
        altitude = sat_alt_data_final[epoch][3] # Read raw mean SMA string
        sat_altitude.append(float(altitude))
    
    # Get the maneuver summary data providers pointer and interface to it.
    sat_deltaV  = satellite.DataProviders.Item("Maneuver Summary")
    sat_deltaV2 = sat_deltaV.QueryInterface(STKObjects.IAgDataPrvInterval)
    sat_deltaV2_Data = sat_deltaV2.Exec(scenario2.StartTime, scenario2.StopTime)
    sat_deltaV2_Array = np.array(sat_deltaV2_Data.DataSets.ToArray())
    
    # Extract the Delta-V values
    deltaV_file = open("output_manoeuvre_stk.txt",'w')
    for thrust in sat_deltaV2_Array:
        thrust_str  = thrust[0] + ' '
        thrust_str += thrust[1] + ' '
        thrust_str += thrust[2] + ' '
        thrust_str += thrust[5] + ' \n'
        deltaV_file.write(thrust_str)
    deltaV_file.close()
    
    # Compute the total Delta-V and the total impulse needed.
    total_deltaV = len(sat_deltaV2_Array) * delta_v
    
    # Total impulse, inclusive of the margin of safety
    total_impulse = total_deltaV * 1000 * sc_mass * maintenance_margin
    Isp = np.linspace(plot_Isp_Min, plot_Isp_Max, 500) # Isp axis, in (s)
    Mf = total_impulse / (Isp*9.81)
    
    # Construct the summary information string objects
    impulse_str = "The total impulse needed: "
    impulse_str = impulse_str + str(total_impulse) + " \n"
    deltaV_str = "The total Delta-V (m/s) needed is "
    deltaV_str = deltaV_str + str(total_deltaV * 1000) + " \n"
    
    # Log the summary information
    summary_file = open("output_summary_stk.txt",'w')
    summary_file.write(lifetime_str)
    summary_file.write(impulse_str)
    summary_file.write(deltaV_str)
    summary_file.close()
    
    # Print the impulse and Delta-V requirements statement.
    print(impulse_str)
    print(deltaV_str)
    
    # Plotting of altitudes
    plt.figure(1)
    plt.title("Plot of Satellite Altitude (km) Against Date-Time")
    plt.ylabel('Altitude (km)')
    plt.xlabel('Date-Time')
    plt.scatter(sat_epochs,sat_altitude,s=4,alpha=0.3)
    plt.grid()
    
    # Plotting of Kozai-Izsak mean semi-major axes
    plt.figure(2)
    plt.title("Plot of Kozai-Izsak Mean Semimajor Axis (km) Against Date-Time")
    plt.ylabel('Kozai-Izsak Mean Semimajor Axis (km)')
    plt.xlabel('Date-Time')
    plt.plot(sat_epochs,sat_mean_sma)
    plt.grid()
    
    # Thruster sizing profile of Isp Against Mass
    plt.figure(3)
    plt.title("Thruster Profiling for Feasible ISPs (s) Against Fuel Mass (kg)")
    plt.ylabel('Mass of Fuel Required (kg)')
    plt.xlabel('Specific Impulse (s)')
    plt.plot(Isp, Mf)
    plt.grid()
    
    """ #########################################################################
    
    Notes: This part of the code reads the text file of all the shortlisted,
    thrusters and then plots them along the Isp-to-fuel-mass sizing chart.
    
    ######################################################################### """
    
    # Now we compare the mission propulsion requirements against thrusters.
    try:
        thr_file = open("thruster_shortlist.txt","r")
    except:
        
        # Otherwise, generate the file
        thr_file = open("thruster_shortlist.txt","w")
        thr_file.write("COMPANY         ")
        thr_file.write("MODEL           ")
        thr_file.write("ISP_S           ")
        thr_file.write("FUEL_MASS_KG    ")
        thr_file.write("THRUST_N        ")
        thr_file.write("END \n")
        thr_file.write("ALIENA          ")
        thr_file.write("MUSIC           ")
        thr_file.write("1000            ")
        thr_file.write("3.000           ")
        thr_file.write("0.004           ")
        thr_file.write("END \n")
        thr_file.close()
        
        # Now, try to open the file
        thr_file = open("thruster_shortlist.txt","r")
        
    
    thr_compn = []
    thr_model = []
    thr_isp_s = []
    thr_fuelm = []
    thr_force = []
    
    for line in thr_file:
        line_split = line.split()
        if line_split[0] != "COMPANY":
            thr_compn.append(str(line_split[0]))
            thr_model.append(str(line_split[1]))
            thr_isp_s.append(float(line_split[2]))
            thr_fuelm.append(float(line_split[3]))
            thr_force.append(str(line_split[4]))
    thr_file.close()
    
    # plot_Isp_Min = 200.0 # N s
    # plot_Isp_Max = 1250.0 # N s
    bwidth = (plot_Isp_Max - plot_Isp_Min)/50
    
    barchart = plt.bar(thr_isp_s, thr_fuelm, width = bwidth, color='green')
    
    # Then, we label each thruster accordingly.
    barcount = 0
    for rect in barchart:
        height = rect.get_height()
        bartext = thr_compn[barcount] + '\n'
        bartext = bartext + thr_model[barcount] + '\n'
        bartext = bartext + thr_force[barcount] + 'N'
        plt.text(rect.get_x() + rect.get_width()/2.0,
                 rect.get_height(),
                 bartext,
                 ha='center', va='bottom')
        barcount += 1
    
    return None
    # END OF SCRIPT
