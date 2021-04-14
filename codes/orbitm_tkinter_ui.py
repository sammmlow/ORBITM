#!/usr/bin/env python3

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
##    This file contains the GUI class which based on Python tkinter.        ##
##    The class will be called in the main OrbitM python file.               ##
##    (No inputs and outputs, this file only holds the GUI class object.     ##
##                                                                           ##
##    Written by Samuel Y. W. Low.                                           ##
##    Last modified 25-11-2020 08:48 AM (+8 GMT)                             ##
##                                                                           ##
###############################################################################
###############################################################################


# IMPORT PUBLIC LIBRARIES
import tkinter as tk
from PIL import Image, ImageTk
from os.path import dirname, abspath

# IMPORT THE MAIN ROUTINES
from codes.orbitm_run_offline import orbm_run_offline
from codes.orbitm_run_stk import orbm_run_stk

class run_gui:

    def __init__(self, master):
        
        # Create the main frame and window.
        master.title('Orbit Maintenance and Propulsion Sizing')
        master.geometry('960x960')
        
        # Initialise the basic text labels (found in the configuration file):
        self.txt0  = 'Choose your orbit simulation program:'
        self.txt1  = 'Start Epoch (e.g. 1-Jan-2012-12:00:00.000):'
        self.txt2  = 'Final Epoch (e.g. 1-Jan-2015-12:00:00.000):'
        self.txt3  = 'Spacecraft Atmospheric Drag Coefficient (Cd):'
        self.txt4  = 'Spacecraft Atmospheric Drag Surface Area (m^2):'
        self.txt5  = 'Spacecraft Albedo Pressure Coefficient (Ck):'
        self.txt6  = 'Spacecraft Albedo Pressure Surface Area (m^2):'
        self.txt7  = 'Spacecraft Radiation Pressure Coefficient (Cr):'
        self.txt8  = 'Spacecraft Radiation Pressure Surface Area (m^2):'
        self.txt9  = 'Orbit Semi-Major Axis (km):'
        self.txt10 = 'Orbit Eccentricity (no units):'
        self.txt11 = 'Orbit Inclination (degrees):'
        self.txt12 = 'Orbit Right Asc. Node (degrees):'
        self.txt13 = 'Orbit Arg. Perigee (degrees):'
        self.txt14 = 'Orbit Mean Anomaly (degrees):'
        self.txt15 = 'Maintenance Tolerance Band (km):'
        self.txt16 = 'Maintenance Mission Margin (1.0=100%):'
        self.txt17 = 'Maintenance for Frozen Repeat Ground Track?'
        self.txt18 = 'Wet Mass of the Spacecraft (kg)'
        self.txt19 = 'For plotting, input x-axis minimum for Isp (s)'
        self.txt20 = 'For plotting, input x-axis maximum for Isp (s)'
        
        # Initialise tkinter variables for the entries corresponding to above.
        self.var00 = tk.IntVar()    # Orbit simulation program choice.
        self.var01 = tk.StringVar() # Start Epoch String
        self.var02 = tk.StringVar() # Final Epoch String
        self.var03 = tk.DoubleVar() # Atmospheric Drag Coefficient (Cd)
        self.var04 = tk.DoubleVar() # Atmospheric Drag Surface Area (m^2)
        self.var05 = tk.DoubleVar() # Albedo Pressure Coefficient (Ck)
        self.var06 = tk.DoubleVar() # Albedo Pressure Surface Area (m^2)
        self.var07 = tk.DoubleVar() # Radiation Pressure Coefficient (Cr)
        self.var08 = tk.DoubleVar() # Radiation Pressure Surface Area (m^2)
        self.var09 = tk.DoubleVar() # Orbit Semi-Major Axis (km)
        self.var10 = tk.DoubleVar() # Orbit Eccentricity (no units)
        self.var11 = tk.DoubleVar() # Orbit Inclination (degrees)
        self.var12 = tk.DoubleVar() # Orbit Right Asc. Node (degrees)
        self.var13 = tk.DoubleVar() # Orbit Arg. Perigee (degrees)
        self.var14 = tk.DoubleVar() # Orbit Mean Anomaly (degrees)
        self.var15 = tk.DoubleVar() # Maintenance Tolerance Band (km)
        self.var16 = tk.DoubleVar() # Maintenance Mission Margin (1.0 = 100%)
        self.var17 = tk.IntVar() # Toggle for normal or FRO Computation
        self.var18 = tk.DoubleVar() # Spacecraft wet mass (kg)
        self.var19 = tk.DoubleVar() # Isp minimum (s) for x-axis
        self.var20 = tk.DoubleVar() # Isp maximum (s) for x-axis
        
        # Define the path to the LEOGPS logo file.
        orbitm_logo = dirname(dirname(abspath(__file__)))
        orbitm_logo = orbitm_logo + '\gui\orbm_logo.png'
        
        # Configure the background image and load the logo.
        image = Image.open( orbitm_logo )
        photo = ImageTk.PhotoImage(image)
        self.logo = tk.Label(image=photo)
        self.logo.image = photo
        self.logo.grid(row=0, column=0, padx=20, pady=20)
        
        # Add a button to read default entries from 'config.txt'.
        self.cfgR = tk.Button(master, text='Load Config', command=self.cfg_R)
        self.cfgR.grid(row=0, column=1, padx=20, pady=2)
        self.cfgR.configure(bg="light blue")
        
        # Add a button to save entries into 'config.txt'.
        self.cfgW = tk.Button(master, text='Save Config', command=self.cfg_W)
        self.cfgW.grid(row=0, column=2, padx=20, pady=2)
        self.cfgW.configure(bg="light blue")
        
        # Add a button to run LEOGPS.
        self.cfgW = tk.Button(master, text='Run ORBITM', command=self.run)
        self.cfgW.grid(row=0, column=3, padx=20, pady=2)
        self.cfgW.configure(bg="light blue")
        
        # Choose the program to run the orbit simulation
        self.label00 = tk.Label(master, text=self.txt0 )
        self.label00.grid(row=1, column=0, padx=20, pady=2, sticky='w')
        self.entry00a = tk.Radiobutton(master, text="Sam's",
                                       variable=self.var00, value=1)
        self.entry00a.grid(row=1, column=1, padx=5, pady=2, sticky='w')
        self.entry00b = tk.Radiobutton(master, text='STK10',
                                       variable=self.var00, value=2)
        self.entry00b.grid(row=1, column=2, padx=5, pady=2, sticky='w')
        self.entry00c = tk.Radiobutton(master, text='STK11',
                                       variable=self.var00, value=3)
        self.entry00c.grid(row=1, column=3, padx=5, pady=2, sticky='w')
        self.errtx00 = tk.Label(master, text='', fg='red')
        self.errtx00.grid(row=1, column=4, padx=5, pady=2, sticky='w')
        
        # Start Epoch String
        self.label01 = tk.Label(master, text=self.txt1 )
        self.label01.grid(row=2, column=0, padx=20, pady=2, sticky='w')
        self.entry01 = tk.Entry(master, width=25, textvariable=self.var01)
        self.entry01.grid(row=2, column=1, padx=5, pady=2, sticky='w',
                          columnspan=2)
        self.errtx01 = tk.Label(master, text='', fg='red' )
        self.errtx01.grid(row=2, column=3, padx=5, pady=2, sticky='w')
        
        # Final Epoch String
        self.label02 = tk.Label(master, text=self.txt2 )
        self.label02.grid(row=3, column=0, padx=20, pady=2, sticky='w')
        self.entry02 = tk.Entry(master, width=25, textvariable=self.var02)
        self.entry02.grid(row=3, column=1, padx=5, pady=2, sticky='w',
                          columnspan=2)
        self.errtx02 = tk.Label(master, text='', fg='red' )
        self.errtx02.grid(row=3, column=3, padx=5, pady=2, sticky='w')
        
        # Atmospheric Drag Coefficient (Cd)
        self.label03 = tk.Label(master, text=self.txt3 )
        self.label03.grid(row=4, column=0, padx=20, pady=2, sticky='w')
        self.entry03 = tk.Entry(master, width=10, textvariable=self.var03)
        self.entry03.grid(row=4, column=1, padx=5, pady=2, sticky='w')
        self.errtx03 = tk.Label(master, text='', fg='red' )
        self.errtx03.grid(row=4, column=3, padx=5, pady=2, sticky='w')
        
        # Atmospheric Drag Surface Area (m^2)
        self.label04 = tk.Label(master, text=self.txt4 )
        self.label04.grid(row=5, column=0, padx=20, pady=2, sticky='w')
        self.entry04 = tk.Entry(master, width=10, textvariable=self.var04)
        self.entry04.grid(row=5, column=1, padx=5, pady=2, sticky='w')
        self.errtx04 = tk.Label(master, text='', fg='red' )
        self.errtx04.grid(row=5, column=3, padx=5, pady=2, sticky='w')
        
        # Albedo Pressure Coefficient (Ck)
        self.label05 = tk.Label(master, text=self.txt5 )
        self.label05.grid(row=6, column=0, padx=20, pady=2, sticky='w')
        self.entry05 = tk.Entry(master, width=10, textvariable=self.var05)
        self.entry05.grid(row=6, column=1, padx=5, pady=2, sticky='w')
        self.errtx05 = tk.Label(master, text='', fg='red' )
        self.errtx05.grid(row=6, column=3, padx=5, pady=2, sticky='w')
        
        # Albedo Pressure Surface Area (m^2)
        self.label06 = tk.Label(master, text=self.txt6 )
        self.label06.grid(row=7, column=0, padx=20, pady=2, sticky='w')
        self.entry06 = tk.Entry(master, width=10, textvariable=self.var06)
        self.entry06.grid(row=7, column=1, padx=5, pady=2, sticky='w')
        self.errtx06 = tk.Label(master, text='', fg='red' )
        self.errtx06.grid(row=7, column=3, padx=5, pady=2, sticky='w')
        
        # Radiation Pressure Coefficient (Cr)
        self.label07 = tk.Label(master, text=self.txt7 )
        self.label07.grid(row=8, column=0, padx=20, pady=2, sticky='w')
        self.entry07 = tk.Entry(master, width=10, textvariable=self.var07)
        self.entry07.grid(row=8, column=1, padx=5, pady=2, sticky='w')
        self.errtx07 = tk.Label(master, text='', fg='red' )
        self.errtx07.grid(row=8, column=3, padx=5, pady=2, sticky='w')
        
        # Radiation Pressure Surface Area (m^2)
        self.label08 = tk.Label(master, text=self.txt8 )
        self.label08.grid(row=9, column=0, padx=20, pady=2, sticky='w')
        self.entry08 = tk.Entry(master, width=10, textvariable=self.var08)
        self.entry08.grid(row=9, column=1, padx=5, pady=2, sticky='w')
        self.errtx08 = tk.Label(master, text='', fg='red' )
        self.errtx08.grid(row=9, column=3, padx=5, pady=2, sticky='w')
        
        # Orbit Semi-Major Axis (km)
        self.label09 = tk.Label(master, text=self.txt9 )
        self.label09.grid(row=10, column=0, padx=20, pady=2, sticky='w')
        self.entry09 = tk.Entry(master, width=10, textvariable=self.var09)
        self.entry09.grid(row=10, column=1, padx=5, pady=2, sticky='w')
        self.errtx09 = tk.Label(master, text='', fg='red' )
        self.errtx09.grid(row=10, column=3, padx=5, pady=2, sticky='w')
        
        # Orbit Eccentricity (no units)
        self.label10 = tk.Label(master, text=self.txt10 )
        self.label10.grid(row=11, column=0, padx=20, pady=2, sticky='w')
        self.entry10 = tk.Entry(master, width=10, textvariable=self.var10)
        self.entry10.grid(row=11, column=1, padx=5, pady=2, sticky='w')
        self.errtx10 = tk.Label(master, text='', fg='red' )
        self.errtx10.grid(row=11, column=3, padx=5, pady=2, sticky='w')
        
        # Orbit Inclination (degrees)
        self.label11 = tk.Label(master, text=self.txt11 )
        self.label11.grid(row=12, column=0, padx=20, pady=2, sticky='w')
        self.entry11 = tk.Entry(master, width=10, textvariable=self.var11)
        self.entry11.grid(row=12, column=1, padx=5, pady=2, sticky='w')
        self.errtx11 = tk.Label(master, text='', fg='red' )
        self.errtx11.grid(row=12, column=3, padx=5, pady=2, sticky='w')
        
        # Orbit Right Asc. Node (degrees)
        self.label12 = tk.Label(master, text=self.txt12 )
        self.label12.grid(row=13, column=0, padx=20, pady=2, sticky='w')
        self.entry12 = tk.Entry(master, width=10, textvariable=self.var12)
        self.entry12.grid(row=13, column=1, padx=5, pady=2, sticky='w')
        self.errtx12 = tk.Label(master, text='', fg='red' )
        self.errtx12.grid(row=13, column=3, padx=5, pady=2, sticky='w')
        
        # Orbit Arg. Perigee (degrees)
        self.label13 = tk.Label(master, text=self.txt13 )
        self.label13.grid(row=14, column=0, padx=20, pady=2, sticky='w')
        self.entry13 = tk.Entry(master, width=10, textvariable=self.var13)
        self.entry13.grid(row=14, column=1, padx=5, pady=2, sticky='w')
        self.errtx13 = tk.Label(master, text='', fg='red' )
        self.errtx13.grid(row=14, column=3, padx=5, pady=2, sticky='w')
        
        # Orbit Mean Anomaly (degrees)
        self.label14 = tk.Label(master, text=self.txt14 )
        self.label14.grid(row=15, column=0, padx=20, pady=2, sticky='w')
        self.entry14 = tk.Entry(master, width=10, textvariable=self.var14)
        self.entry14.grid(row=15, column=1, padx=5, pady=2, sticky='w')
        self.errtx14 = tk.Label(master, text='', fg='red' )
        self.errtx14.grid(row=15, column=3, padx=5, pady=2, sticky='w')
        
        # Maintenance Tolerance Band (km)
        self.label15 = tk.Label(master, text=self.txt15 )
        self.label15.grid(row=16, column=0, padx=20, pady=2, sticky='w')
        self.entry15 = tk.Entry(master, width=10, textvariable=self.var15)
        self.entry15.grid(row=16, column=1, padx=5, pady=2, sticky='w')
        self.errtx15 = tk.Label(master, text='', fg='red' )
        self.errtx15.grid(row=16, column=3, padx=5, pady=2, sticky='w')
        
        # Maintenance Mission Margin (1.0 = 100%)
        self.label16 = tk.Label(master, text=self.txt16 )
        self.label16.grid(row=17, column=0, padx=20, pady=2, sticky='w')
        self.entry16 = tk.Entry(master, width=10, textvariable=self.var16)
        self.entry16.grid(row=17, column=1, padx=5, pady=2, sticky='w')
        self.errtx16 = tk.Label(master, text='', fg='red' )
        self.errtx16.grid(row=17, column=3, padx=5, pady=2, sticky='w')
        
        # Toggle for normal or FRO Computation
        self.label17 = tk.Label(master, text=self.txt17 )
        self.label17.grid(row=18, column=0, padx=20, pady=2, sticky='w')
        self.entry17 = tk.Checkbutton(master, text='True/False',
                                      variable=self.var17)
        self.entry17.grid(row=18, column=1, padx=5, pady=2, sticky='w')
        self.errtx17 = tk.Label(master, text='', fg='red' )
        self.errtx17.grid(row=18, column=3, padx=5, pady=2, sticky='w')
        
        # Spacecraft wet mass (kg)
        self.label18 = tk.Label(master, text=self.txt18 )
        self.label18.grid(row=19, column=0, padx=20, pady=2, sticky='w')
        self.entry18 = tk.Entry(master, width=10, textvariable=self.var18)
        self.entry18.grid(row=19, column=1, padx=5, pady=2, sticky='w')
        self.errtx18 = tk.Label(master, text='', fg='red' )
        self.errtx18.grid(row=19, column=3, padx=5, pady=2, sticky='w')
        
        # Specific impulse Isp minimum for x-axis (s)
        self.label19 = tk.Label(master, text=self.txt19 )
        self.label19.grid(row=20, column=0, padx=20, pady=2, sticky='w')
        self.entry19 = tk.Entry(master, width=15, textvariable=self.var19)
        self.entry19.grid(row=20, column=1, padx=5, pady=2, sticky='w')
        self.errtx19 = tk.Label(master, text='', fg='red' )
        self.errtx19.grid(row=20, column=3, padx=5, pady=2, sticky='w')
        
        # Specific impulse Isp maximum for x-axis (s)
        self.label20 = tk.Label(master, text=self.txt20 )
        self.label20.grid(row=21, column=0, padx=20, pady=2, sticky='w')
        self.entry20 = tk.Entry(master, width=15, textvariable=self.var20)
        self.entry20.grid(row=21, column=1, padx=5, pady=2, sticky='w')
        self.errtx20 = tk.Label(master, text='', fg='red' )
        self.errtx20.grid(row=21, column=3, padx=5, pady=2, sticky='w')
    
    # Method to load default values from the configuration file.
    def cfg_R(self):
        
        cwd = dirname(dirname(abspath(__file__))) # Current working directory
        iwd = cwd + '\config\config.txt' # Inputs files
        inputfile = open(iwd,'r') # Open the config.txt file
        inps = {} # Create a dictionary to store all the input 
        integers = [ 'orbsim' ]
        floats = ['sc_Cd','sc_Ad','sc_Ck','sc_Ak','sc_Cr','sc_Ar',
                  'orb_a','orb_e','orb_i','orb_R','orb_w','orb_m',
                  'orbm_tolr','orbm_marg','sc_mass','isp_min','isp_max']
        
        # Now we parse through the config.txt file.
        for line in inputfile:
            
            # Check for input entry with an 'I', then split and format.
            if line[0] == 'I':
                line_inp = line[3:].split()
                
                # Now, let's try to parse parameters meant to be integers.
                if line_inp[0] in integers:
                    
                    try:
                        inps[ line_inp[0] ] = int(line_inp[1])
                    except ValueError:
                        print('Warning, error when reading '+ line_inp[0] +'!')
                        inps[ line_inp[0] ] = line_inp[1]
                
                # then we parse parameters meant to be floats.
                elif line_inp[0] in floats:
                    
                    try:
                        inps[ line_inp[0] ] = float(line_inp[1])
                    except ValueError:
                        print('Warning, error when reading '+ line_inp[0] +'!')
                        inps[ line_inp[0] ] = line_inp[1]
                        
                # For all other parameters, just log them down as strings.
                else:
                    inps[ line_inp[0] ] = line_inp[1]
        
        # Close the file when done
        inputfile.close()
        
        # Now, we will check for errors in the inputs.
        
        # First, check for the orbit simulation program.
        if inps['orbsim'] == 1:
            self.entry00a.select()
            self.entry00b.deselect()
            self.entry00c.deselect()
            self.errtx00.configure(text='')
        elif inps['orbsim'] == 2:
            self.entry00a.deselect()
            self.entry00b.select()
            self.entry00c.deselect()
            self.errtx00.configure(text='')
        elif inps['orbsim'] == 3:
            self.entry00a.deselect()
            self.entry00b.deselect()
            self.entry00c.select()
            self.errtx00.configure(text='')
        else:
            self.errtx00.configure(text='No program?')
        
        # Check for the first epoch string.            
        if len(inps['tstart']) > 24 or inps['tstart'].count('-') != 3:
            self.errtx01.configure(text='Invalid! Check the epoch string!')
        else:
            self.var01.set(inps['tstart'])
            self.errtx01.configure(text='')
        
        # Check for the final epoch string
        if len(inps['tfinal']) > 24 or inps['tfinal'].count('-') != 3:
            self.errtx02.configure(text='Invalid! Check the epoch string!')
        else:
            self.var02.set(inps['tfinal'])
            self.errtx02.configure(text='')
        
        # Check for Spacecraft Atmospheric Drag Coefficient (Cd)
        if inps['sc_Cd'] > 0.0:
            self.var03.set(inps['sc_Cd'])
            self.errtx03.configure(text='')
        else:
            self.errtx03.configure(text='Error! Must be a positive float!')
        
        # Check for Spacecraft Atmospheric Drag Surface Area (m^2)
        if inps['sc_Ad'] > 0.0:
            self.var04.set(inps['sc_Ad'])
            self.errtx04.configure(text='')
        else:
            self.errtx04.configure(text='Error! Must be a positive float!')
        
        # Input the Spacecraft Albedo Pressure Coefficient (Ck)
        if inps['sc_Ck'] > 0.0:
            self.var05.set(inps['sc_Ck'])
            self.errtx05.configure(text='')
        else:
            self.errtx05.configure(text='Error! Must be a positive float!') 
        
        # Check for Spacecraft Albedo Pressure Surface Area (m^2)
        if inps['sc_Ak'] > 0.0:
            self.var06.set(inps['sc_Ak'])
            self.errtx06.configure(text='')
        else:
            self.errtx06.configure(text='Error! Must be a positive float!') 
        
        # Check for Spacecraft Radiation Pressure Surface Area (m^2)
        if inps['sc_Cr'] > 0.0:
            self.var07.set(inps['sc_Cr'])
            self.errtx07.configure(text='')
        else:
            self.errtx07.configure(text='Error! Must be a positive float!') 
        
        # Check for Spacecraft Radiation Pressure Surface Area (m^2)
        if inps['sc_Ar'] > 0.0:
            self.var08.set(inps['sc_Ar'])
            self.errtx08.configure(text='')
        else:
            self.errtx08.configure(text='Error! Must be a positive float!') 
        
        # Check for Orbit Semi-Major Axis (km)
        if inps['orb_a'] > 6378.140:
            self.var09.set(inps['orb_a'])
            self.errtx09.configure(text='')
        else:
            self.errtx09.configure(text='Satellite is below Earth surface!') 
        
        # Check for Orbit Eccentricity (no units)
        if inps['orb_e'] > 0.0 and inps['orb_e'] < 1.0:
            self.var10.set(inps['orb_e'])
            self.errtx10.configure(text='')
        else:
            self.errtx10.configure(text='Eccentricity must be 0 < e < 1!') 
        
        # Check for Orbit Inclination (degrees)
        if inps['orb_i'] > -180.0 and inps['orb_i'] < 180.0:
            self.var11.set(inps['orb_i'])
            self.errtx11.configure(text='')
        else:
            self.errtx11.configure(text='Inclination must be within +/- 180!')
        
        # Check for Orbit Right Asc. Node (degrees)
        if inps['orb_R'] > 0.0 and inps['orb_R'] <= 360.0:
            self.var12.set(inps['orb_R'])
            self.errtx12.configure(text='')
        else:
            self.errtx12.configure(text='RAAN must be 0 < R <= 360!') 
        
        # Check for Orbit Arg. Perigee (degrees)
        if inps['orb_w'] > 0.0 and inps['orb_w'] <= 360.0:
            self.var13.set(inps['orb_w'])
            self.errtx13.configure(text='')
        else:
            self.errtx13.configure(text='ARGP must be 0 < W <= 360!') 
        
        # Check for Orbit Mean Anomaly (degrees)
        if inps['orb_m'] > 0.0 and inps['orb_m'] <= 360.0:
            self.var14.set(inps['orb_m'])
            self.errtx14.configure(text='')
        else:
            self.errtx14.configure(text='Mean Anomaly must be 0 < W <= 360!') 
        
        # Check for Maintenance Tolerance Band (km)
        if inps['orbm_tolr'] > 0.0:
            self.var15.set(inps['orbm_tolr'])
            self.errtx15.configure(text='')
        else:
            self.errtx15.configure(text='Error! Must be a positive float!') 
        
        # Check for Maintenance Mission Margin (1.0 = 100%)
        if inps['orbm_marg'] >= 1.0:
            self.var16.set(inps['orbm_marg'])
            self.errtx16.configure(text='')
        else:
            self.errtx16.configure(text='Warning! Below 100% mission margin!') 
        
        # Check for Maintenance Mission Margin (1.0 = 100%)
        if inps['orbm_fro_flag'] == 'True':
            self.entry17.select()
            self.errtx17.configure(text='')
        elif inps['orbm_fro_flag'] == 'False':
            self.entry17.deselect()
            self.errtx17.configure(text='')
        else:
            self.errtx17.configure(text='Error! Input must be True/False!')
        
        # Check for the spacecraft wet mass (kg)
        if inps['sc_mass'] > 0.0:
            self.var18.set(inps['sc_mass'])
            self.errtx18.configure(text='')
        else:
            self.errtx18.configure(text='Error! Must be a positive float!')
            
        # Check for the user input Isp x-axis minimum (s)
        if inps['isp_min'] > 0.0:
            self.var19.set(inps['isp_min'])
            self.errtx19.configure(text='')
        else:
            self.errtx19.configure(text='Error! Must be a positive float!')
        
        # Check for the user input Isp x-axis maximum (s)
        if inps['isp_max'] > 0.0:
            self.var20.set(inps['isp_max'])
            self.errtx20.configure(text='')
        else:
            self.errtx20.configure(text='Error! Must be a positive float!')
        
        return None
    
    # Method for writing the entries into the config.txt file.
    def cfg_W(self):
        
        cwd = dirname(dirname(abspath(__file__))) # Current working directory
        iwd = cwd + '\config\config.txt' # Inputs files
        input_r = open(iwd,'r') # Open the config.txt file
        record = [] # Array to record the strings
        error_flag = True # Check if user-defined variables is correct.
        
        # Variables to be recorded based on tkinter entries.
        var_arr = [self.var00, self.var01, self.var02, self.var03, self.var04,
                   self.var05, self.var06, self.var07, self.var08, self.var09,
                   self.var10, self.var11, self.var12, self.var13, self.var14,
                   self.var15, self.var16, self.var17, self.var18, self.var19,
                   self.var20]
        
        # Key values (to be referred).
        key_arr = ['orbsim', 'tstart', 'tfinal',
                   'sc_Cd', 'sc_Ad', 'sc_Ck', 'sc_Ak', 'sc_Cr', 'sc_Ar',
                   'orb_a', 'orb_e', 'orb_i', 'orb_R', 'orb_w', 'orb_m',
                   'orbm_tolr', 'orbm_marg', 'orbm_fro_flag','sc_mass',
                   'isp_min', 'isp_max']
        
        t_f_arr = ['orbm_fro_flag']
        
        # First, check for the orbit simulation program.
        orbprog = self.var00.get()
        if orbprog == 1 or orbprog == 2 or orbprog == 3:
            self.errtx00.configure(text='')
        else:
            self.entry00a.deselect()
            self.entry00b.deselect()
            self.entry00c.deselect()
            self.errtx00.configure(text='No program?')
            error_flag = False
        
        # Check for the first epoch string.
        if len(self.var01.get()) > 24 or self.var01.get().count('-') != 3:
            self.errtx01.configure(text='Invalid! Check the epoch string!')
            error_flag = False
        else:
            self.errtx01.configure(text='')
        
        # Check for the final epoch string
        if len(self.var02.get()) > 24 or self.var02.get().count('-') != 3:
            self.errtx02.configure(text='Invalid! Check the epoch string.')
            error_flag = False
        else:
            self.errtx02.configure(text='')
        
        # Check for Spacecraft Atmospheric Drag Coefficient (Cd)
        if type(self.var03.get()) == float and self.var03.get() > 0.0:
            self.errtx03.configure(text='')
        else:
            self.errtx03.configure(text='Error! Must be a positive float!')
            error_flag = False
        
        # Check for Spacecraft Atmospheric Drag Surface Area (m^2)
        if self.var04.get() > 0.0:
            self.errtx04.configure(text='')
        else:
            self.errtx04.configure(text='Error! Must be a positive float!')
            error_flag = False
        
        # Input the Spacecraft Albedo Pressure Coefficient (Ck)
        if self.var05.get() > 0.0:
            self.errtx05.configure(text='')
        else:
            self.errtx05.configure(text='Error! Must be a positive float!')
            error_flag = False
        
        # Check for Spacecraft Albedo Pressure Surface Area (m^2)
        if self.var06.get() > 0.0:
            self.errtx06.configure(text='')
        else:
            self.errtx06.configure(text='Error! Must be a positive float!')
            error_flag = False
        
        # Check for Spacecraft Radiation Pressure Surface Area (m^2)
        if self.var07.get() > 0.0:
            self.errtx07.configure(text='')
        else:
            self.errtx07.configure(text='Error! Must be a positive float!')
            error_flag = False
        
        # Check for Spacecraft Radiation Pressure Surface Area (m^2)
        if self.var08.get() > 0.0:
            self.errtx08.configure(text='')
        else:
            self.errtx08.configure(text='Error! Must be a positive float!')
            error_flag = False
        
        # Check for Orbit Semi-Major Axis (km)
        if self.var09.get() > 6378.140:
            self.errtx09.configure(text='')
        else:
            self.errtx09.configure(text='Satellite is below Earth surface!')
            error_flag = False
        
        # Check for Orbit Eccentricity (no units)
        if self.var10.get() > 0.0 and self.var10.get() < 1.0:
            self.errtx10.configure(text='')
        else:
            self.errtx10.configure(text='Eccentricity must be 0 < e < 1!')
            error_flag = False
        
        # Check for Orbit Inclination (degrees)
        if self.var11.get() > -180.0 and self.var11.get() < 180.0:
            self.errtx11.configure(text='')
        else:
            self.errtx11.configure(text='Inclination must be within +/- 180!')
            error_flag = False
        
        # Check for Orbit Right Asc. Node (degrees)
        if self.var12.get() > 0.0 and self.var12.get() <= 360.0:
            self.errtx12.configure(text='')
        else:
            self.errtx12.configure(text='RAAN must be 0 < R <= 360!')
            error_flag = False
        
        # Check for Orbit Arg. Perigee (degrees)
        if self.var13.get() > 0.0 and self.var13.get() <= 360.0:
            self.errtx13.configure(text='')
        else:
            self.errtx13.configure(text='ARGP must be 0 < W <= 360!')
            error_flag = False
        
        # Check for Orbit Mean Anomaly (degrees)
        if self.var14.get() > 0.0 and self.var14.get() <= 360.0:
            self.errtx14.configure(text='')
        else:
            self.errtx14.configure(text='Mean Anomaly must be 0 < W <= 360!')
            error_flag = False
        
        # Check for Maintenance Tolerance Band (km)
        if self.var15.get() > 0.0:
            self.errtx15.configure(text='')
        else:
            self.errtx15.configure(text='Error! Must be a positive float!')
            error_flag = False
        
        # Check for Maintenance Mission Margin (1.0 = 100%)
        if self.var16.get() >= 1.0:
            self.errtx16.configure(text='')
        else:
            self.errtx16.configure(text='Warning! Below 100% mission margin!')
        
        # Check for Maintenance Mission Margin (1.0 = 100%)
        if self.var17.get() == 1 or self.var17.get() == 0:
            self.errtx17.configure(text='')
        else:
            self.errtx17.configure(text='Error! Should be 1/0 (T/F)!')
            error_flag = False
        
        # Check for Spacecraft Wet Mass (kg)
        if self.var18.get() > 0.0:
            self.errtx18.configure(text='')
        else:
            self.errtx18.configure(text='Error! Mass must be positive!')
            error_flag = False
            
        # Input Isp minimum for x-axis plot (s)
        if self.var19.get() > 0.0:
            self.errtx19.configure(text='')
        else:
            self.errtx19.configure(text='Error! Isp must be positive!')
            error_flag = False
            
        # Input Isp maximum for x-axis plot (s)
        if self.var20.get() > 0.0:
            self.errtx20.configure(text='')
        else:
            self.errtx20.configure(text='Error! Isp must be positive!')
            error_flag = False
            
        # If error_flag == 1, then no errors found, proceed with overwriting.
        if error_flag == True:
            
            # Now we parse through the config.txt file.
            for line in input_r:
                
                if line[0] == 'I':
                    words = line.split() # Split string into list of words.
                    key   = words[1] # Get the key from config.txt
                    value = words[2] # Get the value from config.txt
                    
                    # Get the updated value based on the index from key list.
                    value_new = str(var_arr[ key_arr.index(key) ].get())
                    
                    if key in t_f_arr:
                        if value_new == '1':
                            value_new = 'True'
                        if value_new == '0':
                            value_new = 'False'
                    
                    line_new  = line.replace(value, value_new)
                
                else:
                    line_new = line
                
                # Now, record the entries.
                record.append(line_new)
            
            # Close the file when done
            input_r.close()
            
            # Now, we open and overwrite the config.txt file.
            input_w = open(iwd,'w') # Open the config.txt file
            
            for text in record:
                input_w.write(text)
            
            input_w.close()
        
        return error_flag
    
    def run(self):
                
        err_flag = self.cfg_W() # Write and save to config file before running
        
        # Initialize all the input parameters with sensible variable names
        orbm_mode = self.var00.get()
        tstart_t = self.var01.get().replace('-',' ')
        tfinal_t = self.var02.get().replace('-',' ')
        sc_Cd, sc_area_d = self.var03.get(), self.var04.get()
        sc_Ck, sc_area_a = self.var05.get(), self.var06.get()
        sc_Cr, sc_area_r = self.var07.get(), self.var08.get()
        orb_a, orb_e = self.var09.get(), self.var10.get()
        orb_i, orb_R = self.var11.get(), self.var12.get()
        orb_w, orb_m = self.var13.get(), self.var14.get()
        maintenance_tolerance = self.var15.get()
        maintenance_margin = self.var16.get()
        maintenance_fro = self.var17.get()
        sc_mass = self.var18.get()
        isp_min = self.var19.get()
        isp_max = self.var20.get()
        
        
        # If using Sam's orbit maintenance offline program
        if orbm_mode == 1 and err_flag == True:
            orbm_run_offline(tstart_t, tfinal_t,
                             sc_Cd, sc_area_d,
                             sc_Ck, sc_area_a,
                             sc_Cr, sc_area_r,
                             orb_a, orb_e, orb_i, orb_R, orb_w, orb_m,
                             maintenance_tolerance,
                             maintenance_margin,
                             maintenance_fro,
                             sc_mass, isp_min, isp_max)
        
        # If using STK 10 or STK 11:
        elif (orbm_mode == 2 or orbm_mode == 3) and err_flag == True:
            orbm_run_stk(orbm_mode, tstart_t, tfinal_t,
                         sc_Cd, sc_area_d,
                         sc_Ck, sc_area_a,
                         sc_Cr, sc_area_r,
                         orb_a, orb_e, orb_i, orb_R, orb_w, orb_m,
                         maintenance_tolerance,
                         maintenance_margin,
                         maintenance_fro,
                         sc_mass, isp_min, isp_max)
        
        return None
