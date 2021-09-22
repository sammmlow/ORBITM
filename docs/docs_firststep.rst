..
   ###########################################################################
   ###########################################################################
   ##                                                                       ##
   ##     _____ ___  ___  ___  _____      __  __                            ##
   ##    |  _  | _ \| _ \|_ _||_   _|    |  \/  |                           ##
   ##    | |_| |   <| _ < | |   | |   _  | \  / |                           ##
   ##    |_____|_|\_|___/|___|  |_|  |_| |_|\/|_|                           ##
   ##                                                v 1.1                  ##
   ##                                                                       ##
   ###########################################################################
   ###########################################################################

.. image:: /_static/orbitm_logo.png

|

First Steps
===========

ORBITM can be launched simply by running the Python file **orbitm.py** in the main directory (equivalent to the directory on the master branch on ORBITM's GitHub). You should see a GUI that looks like the one below, pop up:

.. figure:: /_images/orbmgui_blank.jpg
    :align: center
    
    Figure 2.1 - ORBITM Graphical User Interface

As a first step, you can simply load the default input parameters by clicking "Load Config". The parameters are actually saved in the **config.txt** file. Configuration parameters should be loaded through the GUI, and not manually updated through the **config.txt** file since then the software can't check for string formatting errors (e.g. a negative drag coefficient that was typed into the config file manually and by accident would crash the program).

.. warning:: An orbit eccentricity of exactly zero is not allowed, to prevent potential divide-by-zero issues during calculations. For near-perfectly circular orbits, set an arbitrarily small eccentricity number (e.g., 0.00001).

Second, you can choose which lifetime and maintenance simulator you would like to use. By default, the native **"Sam's Simulator"** will be selected. If STK is selected, remember that **a valid STK license for STK Integration and STK Astrogator is required.**

Third, proceed to set the remainder of your orbit and spacecraft parameters. These parameters are self-explanatory. Some parameters are specific to ORBITM and are defined below:

:Maintenance Tolerance Band: Threshold altitude that triggers an orbit raise (km).

:Frozen Repeat Ground Track: Option to consider frozen repeat orbit maintenance.

.. note:: The altitude value used as the reference with the **Maintenance Tolerance Band** is defined as the geodetic altitude. In **"Sam's Simulator"** mode, this altitude is the Keplerian semi-major axis minus Earth's equatorial radius. However, STK-mode accounts for the oblateness of the Earth, and the presence of other perturbing forces. Thus, in STK mode, the Earth-centered altitude often oscillates with the latitude and with time, even for a "circular orbit".
    
	In order to prevent triggering a manoeuvre on the wrong rising or falling edge of the altitude profile, the STK mode triggers manoeuvres based on the **Kozai-Izsak Mean Semi-Major Axis** rather than the nominal altitude. It is acknowledged that this changes the definition for the "threshold altitude", especially for highly elliptical orbits running in STK.

If the option for a **Frozen Repeat Ground Track** was selected, then the orbit raise will happen from the negative threshold below its initial altitude to a positive threshold above the initial altitude, and the thrusts will be done in perigee-apogee pairs. This minimises secular perturbations to the orbit mean eccentricity and the orbit nodal period. For example, if the **Frozen Repeat Ground Track** option was selected, with an initial altitude of 500km, and a threshold of 5km, the orbit raise will be triggered at 495km, and raised to 505km in two thrust pairs. Else, in the non-frozen repeat orbit maintenance case, the orbit raise will be triggered at 495km, and raised back to the original 500km.

Finally, before running the simulation (**Run ORBITM**), check the **"thrusters.txt"** file on the main ORBITM directory, and add thruster specifications that you wish to size your missions against. Some thrusters are written inside this shortlist as an example. The purpose of this shortlist is to provide the program with the thruster specifications, so that the mission planner can graphically compare its Isp and fuel capacity to the suitability of the mission (as an output of the program).

Now, you can run ORBITM. If you are a licensed STK user and you would like to use ORBITM to automate an orbit maintenance simulation in STK for you, make sure there are **no** running instances of STK before you run ORBITM in the STK 10 or 11 mode.

.. note:: ORBITM considers only atmospheric drag and not albedo and radiation effects.
