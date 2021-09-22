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
##    This file contains the GUI class which based on Python tkinter.        ##
##    The class will be called in the main OrbitM python file.               ##
##    (No inputs and outputs, this file only holds the GUI class object.     ##
##                                                                           ##
##    Written by Samuel Y. W. Low.                                           ##
##    First created 27-Nov-2020 12:00 PM (+8 GMT)                            ##
##    Last modified 19-Sep-2021 22:27 PM (-7 GMT)                            ##
##                                                                           ##
###############################################################################
###############################################################################

# Import global libraries
import datetime
import numpy as np
import tkinter as tk
from PIL import Image, ImageTk
from os.path import dirname, abspath
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.backends.backend_tkagg import NavigationToolbar2Tk

# Import local libraries
from source.orbmrun import orbmrun
from source.orbmstk import orbmstk

class RunGUI:
    
    '''This class represents the entire ORBITM GUI, as a TKinter object.
    The constructor takes in a single tkinter.Tk() object as the root GUI.
    All buttons in the GUI are linked to the methods described below.
    
    Methods
    -------
    cfg_R( self )
        This method does two things. First, this method checks that all inputs
        in config.txt are correct. Second, it copies the input parameters into
        the GUI (as TKinter variables).
    cfg_W( self )
        This method does two things. First, this method checks that all inputs
        in the GUI are correct. Second, it copies the GUI parameters into the
        config.txt file, overwriting it.
    clr( self )
        Clears all existing plots in the ORBITM GUI.
    run( self )
        Run the ORBITM program using the orbmrun.py script (or the orbmstk.py
        script if the user chooses and has a valid STK Astrogator license.
    '''

    def __init__(self, master):
        
        '''
        Loads the GUI class object, taking a tkinter.Tk() object as input.
        
        Example initialisation:
        >> root = tkinter.Tk()
        >> root_gui = run_gui( root )
        >> root.mainloop()
        '''
        
        # Create the main frame and window.
        master.title('ORBIT.M v1.1 - Orbit Maintenance and Propulsion Sizing')
        screen_height = master.winfo_screenheight()
        screen_width = master.winfo_screenwidth()
        master_geometry  = str(int(screen_width*0.85)) + 'x'
        master_geometry += str(int(screen_height*0.85))
        master.geometry(master_geometry)
        
        #####################################################################
        #####################################################################
        ###                                                               ###
        ###         Initialisation of text labels and variables           ###
        ###                                                               ###
        #####################################################################
        #####################################################################
        
        # Initialise the basic text labels (found in the configuration file):
        self.txt00 = 'Choose the orbit simulation program'
        self.txt01 = 'Epoch start (e.g. 1-Jan-2012-12:00:00)'
        self.txt02 = 'Epoch final (e.g. 1-Jan-2015-12:00:00)'
        self.txt03 = 'Spacecraft Drag Coefficient (Cd)'
        self.txt04 = 'Spacecraft Drag Surface Area (m2)'
        self.txt05 = 'Orbit Semi-Major Axis (km)'
        self.txt06 = 'Orbit Eccentricity (no units)'
        self.txt07 = 'Orbit Inclination (degrees)'
        self.txt08 = 'Orbit Right Asc. Node (degrees)'
        self.txt09 = 'Orbit Arg. Perigee (degrees)'
        self.txt10 = 'Orbit Mean Anomaly (degrees)'
        self.txt11 = 'Maintenance Tolerance Band (km)'
        self.txt12 = 'Maintain for Repeat Ground Track?'
        self.txt13 = 'Wet Mass of the Spacecraft (kg)'
        self.txt14 = 'Minimum X-Axis Scale for Isp (s) Plot'
        self.txt15 = 'Maximum X-Axis Scale for Isp (s) Plot'
        
        # Initialise tkinter variables for the entries corresponding to above.
        self.var00 = tk.IntVar()    # Orbit simulation program choice.
        self.var01 = tk.StringVar() # Start Epoch String
        self.var02 = tk.StringVar() # Final Epoch String
        self.var03 = tk.DoubleVar() # Atmospheric Drag Coefficient (Cd)
        self.var04 = tk.DoubleVar() # Atmospheric Drag Surface Area (m^2)
        self.var05 = tk.DoubleVar() # Orbit Semi-Major Axis (km)
        self.var06 = tk.DoubleVar() # Orbit Eccentricity (no units)
        self.var07 = tk.DoubleVar() # Orbit Inclination (degrees)
        self.var08 = tk.DoubleVar() # Orbit Right Asc. Node (degrees)
        self.var09 = tk.DoubleVar() # Orbit Arg. Perigee (degrees)
        self.var10 = tk.DoubleVar() # Orbit Mean Anomaly (degrees)
        self.var11 = tk.DoubleVar() # Maintenance Tolerance Band (km)
        self.var12 = tk.IntVar()    # Toggle for normal or FRO Computation
        self.var13 = tk.DoubleVar() # Spacecraft wet mass (kg)
        self.var14 = tk.DoubleVar() # Isp minimum (s) for x-axis
        self.var15 = tk.DoubleVar() # Isp maximum (s) for x-axis
        
        #####################################################################
        #####################################################################
        ###                                                               ###
        ###             Configure the software logo display               ###
        ###                                                               ###
        #####################################################################
        #####################################################################
        
        # Define the path to the ORBITM logo file.
        orbitm_logo = dirname(dirname(abspath(__file__)))
        orbitm_logo = orbitm_logo + '\docs\_static\orbitm_logo.png'
        
        # Get the users current screen height.
        screen_height = master.winfo_screenheight()
        
        # Configure the background image and load the logo.
        image = Image.open( orbitm_logo )
        photo = ImageTk.PhotoImage(image)
        image_h = photo.height()
        image_w = photo.width()
        logo_scale = image_w / image_h
        logo_height = int(screen_height/6)
        logo_width = int(logo_height*logo_scale)
        image_resize = image.resize(( logo_width, logo_height ))
        image_logo = ImageTk.PhotoImage(image_resize)
        self.logo = tk.Label(image=image_logo)
        self.logo.image = image_logo
        self.logo.grid(row=0, column=0, padx=20, pady=20, sticky='w', 
                       rowspan=2, columnspan=4)
        
        #####################################################################
        #####################################################################
        ###                                                               ###
        ###         Add basic buttons for LOAD, SAVE, CLEAR, RUN          ###
        ###                                                               ###
        #####################################################################
        #####################################################################
        
        # Add a button to read default entries from 'config.txt'.
        self.cfgR = tk.Button(master, text='Load Config', command=self.cfg_R)
        self.cfgR.grid(row=1, column=5, padx=20, pady=2)
        self.cfgR.configure(bg="light blue")
        
        # Add a button to save entries into 'config.txt'.
        self.cfgW = tk.Button(master, text='Save Config', command=self.cfg_W)
        self.cfgW.grid(row=1, column=6, padx=20, pady=2)
        self.cfgW.configure(bg="light blue")
        
        # Add a button to clear the relative orbit plots.
        self.clrBtn = tk.Button(master, text='Clear Plots', command=self.clr)
        self.clrBtn.grid(row=1, column=7, padx=20, pady=5)
        self.clrBtn.configure(bg="light blue")
        
        # Add a button to run ORBITM.
        self.cfgW = tk.Button(master, text='Run ORBITM', command=self.run)
        self.cfgW.grid(row=1, column=8, padx=20, pady=2)
        self.cfgW.configure(bg="light blue")
        
        #####################################################################
        #####################################################################
        ###                                                               ###
        ###               Add labels and main data entries                ###
        ###                                                               ###
        #####################################################################
        #####################################################################
        
        # Choose the program to run the orbit simulation
        self.label00 = tk.Label(master, text=self.txt00 )
        self.label00.grid(row=2, column=0, padx=40, pady=2, sticky='w')
        self.entry00a = tk.Radiobutton(master, text="Sam's",
                                       variable=self.var00, value=1)
        self.entry00a.grid(row=2, column=1, padx=3, pady=2, sticky='w')
        self.entry00b = tk.Radiobutton(master, text='STK10',
                                       variable=self.var00, value=2)
        self.entry00b.grid(row=2, column=2, padx=3, pady=2, sticky='w')
        self.entry00c = tk.Radiobutton(master, text='STK11',
                                       variable=self.var00, value=3)
        self.entry00c.grid(row=2, column=3, padx=3, pady=2, sticky='w')
        self.errtx00 = tk.Label(master, text='', fg='red' )
        self.errtx00.grid(row=2, column=4, padx=5, pady=2, sticky='w')
        
        # Start Epoch String
        self.label01 = tk.Label(master, text=self.txt01 )
        self.label01.grid(row=3, column=0, padx=40, pady=2, sticky='w')
        self.entry01 = tk.Entry(master, width=30, textvariable=self.var01)
        self.entry01.grid(row=3, column=1, padx=5, pady=2, sticky='w',
                          columnspan=3)
        self.errtx01 = tk.Label(master, text='', fg='red' )
        self.errtx01.grid(row=3, column=4, padx=5, pady=2, sticky='w')
        
        # Final Epoch String
        self.label02 = tk.Label(master, text=self.txt02 )
        self.label02.grid(row=4, column=0, padx=40, pady=2, sticky='w')
        self.entry02 = tk.Entry(master, width=30, textvariable=self.var02)
        self.entry02.grid(row=4, column=1, padx=5, pady=2, sticky='w',
                          columnspan=3)
        self.errtx02 = tk.Label(master, text='', fg='red' )
        self.errtx02.grid(row=4, column=4, padx=5, pady=2, sticky='w')
        
        # Atmospheric Drag Coefficient (Cd)
        self.label03 = tk.Label(master, text=self.txt03 )
        self.label03.grid(row=5, column=0, padx=40, pady=2, sticky='w')
        self.entry03 = tk.Entry(master, width=8, textvariable=self.var03)
        self.entry03.grid(row=5, column=1, padx=5, pady=2, sticky='w')
        self.errtx03 = tk.Label(master, text='', fg='red' )
        self.errtx03.grid(row=5, column=4, padx=5, pady=2, sticky='w')
        
        # Atmospheric Drag Surface Area (m^2)
        self.label04 = tk.Label(master, text=self.txt04 )
        self.label04.grid(row=6, column=0, padx=40, pady=2, sticky='w')
        self.entry04 = tk.Entry(master, width=8, textvariable=self.var04)
        self.entry04.grid(row=6, column=1, padx=5, pady=2, sticky='w')
        self.errtx04 = tk.Label(master, text='', fg='red' )
        self.errtx04.grid(row=6, column=4, padx=5, pady=2, sticky='w')
        
        # Orbit Semi-Major Axis (km)
        self.label05 = tk.Label(master, text=self.txt05 )
        self.label05.grid(row=7, column=0, padx=40, pady=2, sticky='w')
        self.entry05 = tk.Entry(master, width=13, textvariable=self.var05)
        self.entry05.grid(row=7, column=1, padx=5, pady=2, sticky='w',
                          columnspan = 2)
        self.errtx05 = tk.Label(master, text='', fg='red' )
        self.errtx05.grid(row=7, column=4, padx=5, pady=2, sticky='w')
        
        # Orbit Eccentricity (no units)
        self.label06 = tk.Label(master, text=self.txt06 )
        self.label06.grid(row=8, column=0, padx=40, pady=2, sticky='w')
        self.entry06 = tk.Entry(master, width=13, textvariable=self.var06)
        self.entry06.grid(row=8, column=1, padx=5, pady=2, sticky='w',
                          columnspan = 2)
        self.errtx06 = tk.Label(master, text='', fg='red' )
        self.errtx06.grid(row=8, column=4, padx=5, pady=2, sticky='w')
        
        # Orbit Inclination (degrees)
        self.label07 = tk.Label(master, text=self.txt07 )
        self.label07.grid(row=9, column=0, padx=40, pady=2, sticky='w')
        self.entry07 = tk.Entry(master, width=13, textvariable=self.var07)
        self.entry07.grid(row=9, column=1, padx=5, pady=2, sticky='w',
                          columnspan = 2)
        self.errtx07 = tk.Label(master, text='', fg='red' )
        self.errtx07.grid(row=9, column=4, padx=5, pady=2, sticky='w')
        
        # Orbit Right Asc. Node (degrees)
        self.label08 = tk.Label(master, text=self.txt08 )
        self.label08.grid(row=10, column=0, padx=40, pady=2, sticky='w')
        self.entry08 = tk.Entry(master, width=13, textvariable=self.var08)
        self.entry08.grid(row=10, column=1, padx=5, pady=2, sticky='w',
                          columnspan = 2)
        self.errtx08 = tk.Label(master, text='', fg='red' )
        self.errtx08.grid(row=10, column=4, padx=5, pady=2, sticky='w')
        
        # Orbit Arg. Perigee (degrees)
        self.label09 = tk.Label(master, text=self.txt09 )
        self.label09.grid(row=11, column=0, padx=40, pady=2, sticky='w')
        self.entry09 = tk.Entry(master, width=13, textvariable=self.var09)
        self.entry09.grid(row=11, column=1, padx=5, pady=2, sticky='w',
                          columnspan = 2)
        self.errtx09 = tk.Label(master, text='', fg='red' )
        self.errtx09.grid(row=11, column=4, padx=5, pady=2, sticky='w')
        
        # Orbit Mean Anomaly (degrees)
        self.label10 = tk.Label(master, text=self.txt10 )
        self.label10.grid(row=12, column=0, padx=40, pady=2, sticky='w')
        self.entry10 = tk.Entry(master, width=13, textvariable=self.var10)
        self.entry10.grid(row=12, column=1, padx=5, pady=2, sticky='w',
                          columnspan = 2)
        self.errtx10 = tk.Label(master, text='', fg='red' )
        self.errtx10.grid(row=12, column=4, padx=5, pady=2, sticky='w')
        
        # Maintenance Tolerance Band (km)
        self.label11 = tk.Label(master, text=self.txt11 )
        self.label11.grid(row=13, column=0, padx=40, pady=2, sticky='w')
        self.entry11 = tk.Entry(master, width=8, textvariable=self.var11)
        self.entry11.grid(row=13, column=1, padx=5, pady=2, sticky='w')
        self.errtx11 = tk.Label(master, text='', fg='red' )
        self.errtx11.grid(row=13, column=4, padx=5, pady=2, sticky='w')
        
        # Toggle for normal or FRO Computation
        self.label12 = tk.Label(master, text=self.txt12 )
        self.label12.grid(row=14, column=0, padx=40, pady=2, sticky='w')
        self.entry12 = tk.Checkbutton(master, text='True/False',
                                      variable=self.var12)
        self.entry12.grid(row=14, column=1, padx=5, pady=2, sticky='w')
        self.errtx12 = tk.Label(master, text='', fg='red' )
        self.errtx12.grid(row=14, column=4, padx=5, pady=2, sticky='w')
        
        # Spacecraft wet mass (kg)
        self.label13 = tk.Label(master, text=self.txt13 )
        self.label13.grid(row=15, column=0, padx=40, pady=2, sticky='w')
        self.entry13 = tk.Entry(master, width=8, textvariable=self.var13)
        self.entry13.grid(row=15, column=1, padx=5, pady=2, sticky='w')
        self.errtx13 = tk.Label(master, text='', fg='red' )
        self.errtx13.grid(row=15, column=4, padx=5, pady=2, sticky='w')
        
        # Specific impulse Isp minimum for x-axis (s)
        self.label14 = tk.Label(master, text=self.txt14 )
        self.label14.grid(row=16, column=0, padx=40, pady=2, sticky='w')
        self.entry14 = tk.Entry(master, width=8, textvariable=self.var14)
        self.entry14.grid(row=16, column=1, padx=5, pady=2, sticky='w')
        self.errtx14 = tk.Label(master, text='', fg='red' )
        self.errtx14.grid(row=16, column=4, padx=5, pady=2, sticky='w')
        
        # Specific impulse Isp maximum for x-axis (s)
        self.label15 = tk.Label(master, text=self.txt15 )
        self.label15.grid(row=17, column=0, padx=40, pady=2, sticky='w')
        self.entry15 = tk.Entry(master, width=8, textvariable=self.var15)
        self.entry15.grid(row=17, column=1, padx=5, pady=2, sticky='w')
        self.errtx15 = tk.Label(master, text='', fg='red' )
        self.errtx15.grid(row=17, column=4, padx=5, pady=2, sticky='w')
        
        #####################################################################
        #####################################################################
        ###                                                               ###
        ###              Configure the plot area in the GUI               ###
        ###                                                               ###
        #####################################################################
        #####################################################################
        
        # Now, we add a sub-frame in the tkinter GUI so that we can embed the
        # the output plots for orbit maintenance and decay
        self.toolbarFrame = tk.Frame(master)
        self.toolbarFrame.grid(row=2, column=5, padx=20, pady=10,
                               columnspan=4, rowspan=18, sticky='s')
        
        # Create the 2D matplotlib figure, with three subplots.
        self.Fig = Figure(figsize=(7,5), dpi = master.winfo_fpixels('2.0c'),
                          linewidth=8, edgecolor="#DDDDDD")
        self.Axis211 = self.Fig.add_subplot(211) # Plot propulsion sizing
        self.Axis223 = self.Fig.add_subplot(223) # Plot altitude profile
        self.Axis224 = self.Fig.add_subplot(224) # Plot SMA profile
        self.Fig.set_tight_layout(True)
        self.FigPlot = FigureCanvasTkAgg(self.Fig, self.toolbarFrame)
        self.FigPlot.get_tk_widget().pack(expand=True)
        
        # Plotting of altitudes (titles and axis only)
        self.Axis223.set_ylabel('Altitude (km)')
        self.Axis223.set_xlabel('Date-Time')
        
        # Plotting of Kozai-Izsak mean semi-major axes (titles and axis only)
        self.Axis224.set_ylabel('Mean Semimajor Axis (km)')
        self.Axis224.set_xlabel('Date-Time')
        
        # Thruster sizing profile of Isp Against Mass (titles and axis only)
        self.Axis211.set_ylabel('Mass of Fuel Required (kg)')
        self.Axis211.set_xlabel('Specific Impulse (s)')
        
        # At this point, you can insert plots if you want. For example,
        # self.orbAxis.scatter([1,2,3],[1,2,3])
        self.FigPlot.draw()
        
        # Add the matplotlib navigation toolbar.
        self.toolbar = NavigationToolbar2Tk(self.FigPlot, self.toolbarFrame)
        self.toolbar.update()
        
        #####################################################################
        #####################################################################
        ###                                                               ###
        ###    Finally, define string containers for error and warning    ###
        ###    messages to inform the user if input conditions violate    ###
        ###    formatting or physical principles. If the length of this   ###
        ###    string variable > 0, then it triggers an error message.    ###
        ###                                                               ###
        #####################################################################
        #####################################################################
        
        self.error_msgprint = '' # Error message to print.
    
    #########################################################################
    #########################################################################
    ###                                                                   ###
    ###          Method to load default values from config.txt            ###
    ###                                                                   ###
    #########################################################################
    #########################################################################
    
    def cfg_R(self):
        
        # First, ask the user if he/she wishes to proceed.
        cfg_R_msg = 'Load parameters from the "config.txt" file? \n'
        cfg_R_msg += '(This will overwrite existing inputs in the GUI!)'
        cfg_R_ask = tk.messagebox.askyesno('Load Config', cfg_R_msg)
        if cfg_R_ask == False:
            return None
        
        # Else, continue with loading the configuration file.
        cwd = dirname(dirname(abspath(__file__))) # Current working directory
        iwd = cwd + '\config\config.txt' # Inputs files
        inputfile = open(iwd,'r') # Open the config.txt file
        inps = {} # Create a dictionary to store all the input 
        integers = [ 'orbsim' ]
        floats = ['sc_Cd','sc_Ad','orb_a','orb_e','orb_i','orb_R','orb_w',
                  'orb_m','orbm_tolr','sc_mass','isp_min','isp_max']
        
        #####################################################################
        #####################################################################
        ###                                                               ###
        ###   Parsing through the config.txt file to extract parameters   ###
        ###                                                               ###
        #####################################################################
        #####################################################################
        
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
                        errmsg = 'Error, expected an integer when reading '
                        errmsg = errmsg + line_inp[0] + ' in config.txt! \n'
                        print(errmsg)
                        self.error_msgprint += errmsg
                        inps[ line_inp[0] ] = 'Invalid'
                
                # then we parse parameters meant to be floats.
                elif line_inp[0] in floats:
                    
                    try:
                        inps[ line_inp[0] ] = float(line_inp[1])
                    except ValueError:
                        errmsg = 'Error, expected a float when reading '
                        errmsg = errmsg + line_inp[0] + ' in config.txt! \n'
                        print(errmsg)
                        self.error_msgprint += errmsg
                        inps[ line_inp[0] ] = 'Invalid'
                        
                # For all other parameters, just log them down as strings.
                else:
                    inps[ line_inp[0] ] = line_inp[1]
        
        # Close the file when done
        inputfile.close()
        
        # Prepare a dictionary to convert month strings into integers.
        months_dict = {'Jan':1, 'Feb':2, 'Mar':3, 'Apr':4,
                       'May':5, 'Jun':6, 'Jul':7, 'Aug':8,
                       'Sep':9, 'Oct':10,'Nov':11,'Dec':12}
        
        #####################################################################
        #####################################################################
        ###                                                               ###
        ###  Parsing through the inputs dictionary to verify parameters   ###
        ###                                                               ###
        #####################################################################
        #####################################################################
        
        # 0. First, check for the orbit simulation program (no `errtx00`).
        
        if inps['orbsim'] == 1:
            errmsg = ''
            self.entry00a.select()
            self.entry00b.deselect()
            self.entry00c.deselect()
            self.errtx00.configure(text='')
        elif inps['orbsim'] == 2:
            errmsg = ''
            self.entry00a.deselect()
            self.entry00b.select()
            self.entry00c.deselect()
            self.errtx00.configure(text='')
        elif inps['orbsim'] == 3:
            errmsg = ''
            self.entry00a.deselect()
            self.entry00b.deselect()
            self.entry00c.select()
            self.errtx00.configure(text='')
        else:
            errmsg = 'Invalid simulation option! Check config.txt! \n'
            self.entry00a.deselect()
            self.entry00b.deselect()
            self.entry00c.deselect()
            self.errtx00.configure(text='!')
        self.error_msgprint += errmsg
        
        #####################################################################
        #####################################################################
        
        # 1. Check for the first epoch string.
        
        self.var01.set(inps['tstart'])
        if inps['tstart'].count('-') == 3:
            
            inp_v01s = inps['tstart'].split('-')
            
            # Check if it can be converted into a datetime object.
            try:
                time_str1 = inp_v01s[3].split(':')
                inp_v01d = datetime.datetime(int(inp_v01s[2]),
                                             int(months_dict[inp_v01s[1]]),
                                             int(inp_v01s[0]),
                                             int(time_str1[0]),
                                             int(time_str1[1]),
                                             int(float(time_str1[2])))
                errmsg = ''
                self.errtx01.configure(text='')
            
            # If not, throw an exception and add it to the error log.
            except:
                errmsg = 'Error! Invalid date and time parameters! \n'
                self.errtx01.configure(text='!')
        
        # Else, throw a formatting error.
        else:
            errmsg = 'Error! Invalid date time format in config.txt! \n'
            self.errtx01.configure(text='!')
        self.error_msgprint += errmsg
        
        #####################################################################
        #####################################################################
        
        # 2. Check for the final epoch string
        
        self.var02.set(inps['tfinal'])
        if inps['tfinal'].count('-') == 3:
            
            inp_v02s = inps['tfinal'].split('-')
            
            # Check if it can be converted into a datetime object.
            try:
                time_str2 = inp_v02s[3].split(':')
                inp_v02d = datetime.datetime(int(inp_v02s[2]),
                                             int(months_dict[inp_v02s[1]]),
                                             int(inp_v02s[0]),
                                             int(time_str2[0]),
                                             int(time_str2[1]),
                                             int(float(time_str2[2])))
                
                # Check if the final epoch is after the initial epoch.
                if inp_v02d <= inp_v01d:
                    errmsg = 'Error! The epoch final is before start! \n'
                    self.errtx02.configure(text='!')
                else:
                    errmsg = ''
                    self.errtx02.configure(text='')
            
            # If not, throw an exception and add it to the error log.
            except:
                errmsg = 'Error! Invalid date and time parameters! \n'
                self.errtx02.configure(text='!')
        
        # Else, throw a formatting error.
        else:
            errmsg = 'Error! Invalid date time format in config.txt! \n'
            self.errtx02.configure(text='!')
        self.error_msgprint += errmsg
        
        #####################################################################
        #####################################################################
        
        # 3. Check for Spacecraft Atmospheric Drag Coefficient (Cd)
        
        self.var03.set(inps['sc_Cd'])
        if type(inps['sc_Cd']) == float and inps['sc_Cd'] > 0.0:
            errmsg = ''
            self.errtx03.configure(text='')
        else:
            errmsg = 'Error! Drag coefficient must be a positive float! \n'
            self.errtx03.configure(text='!')
        self.error_msgprint += errmsg
                
        #####################################################################
        #####################################################################
        
        # 4. Check for Spacecraft Atmospheric Drag Surface Area (m^2)
        
        self.var04.set(inps['sc_Ad'])
        if type(inps['sc_Ad']) == float and inps['sc_Ad'] > 0.0:
            errmsg = ''
            self.errtx04.configure(text='')
        else:
            errmsg = 'Error! Drag surface area must be a positive float! \n'
            self.errtx04.configure(text='!')
        self.error_msgprint += errmsg
        
        #####################################################################
        #####################################################################
        
        # 5. Check for Orbit Semi-Major Axis (km)
        
        self.var05.set(inps['orb_a'])
        if type(inps['orb_a']) != float:
            errmsg = 'Error! Semi-major axis must be a float! \n'
            self.errtx05.configure(text='!')
        elif inps['orb_a'] < 6378.14:
            errmsg = 'Error! Semi-major axis below Earth surface! \n'
            self.errtx05.configure(text='!')
        elif inps['orb_a'] > 385000.0:
            errmsg = 'Error! Semi-major axis beyond Earth orbit! \n'
            self.errtx05.configure(text='!')
        else:
            errmsg = ''
            self.errtx05.configure(text='')
        self.error_msgprint += errmsg
        
        #####################################################################
        #####################################################################
        
        # 6. Check for Orbit Eccentricity (no units)
        
        self.var06.set(inps['orb_e'])
        if type(inps['orb_e']) != float:
            errmsg = 'Error! Eccentricity must be a float! \n'
            self.errtx06.configure(text='!')
        elif inps['orb_e'] < 0:
            errmsg = 'Error! Eccentricity cannot be < 0! \n'
            self.errtx06.configure(text='!')
        elif inps['orb_e'] >= 1.0:
            errmsg = 'Error! Eccentricity cannot be >= 1! \n'
            self.errtx06.configure(text='!')
        elif type(inps['orb_a']) == float:
            if ( ( 1 - inps['orb_e'] ) * inps['orb_a'] ) < 6378.14:
                errmsg = 'Error! Perigee altitude below Earth surface! \n'
                self.errtx06.configure(text='!')
            else:
                errmsg = ''
                self.errtx06.configure(text='')
        else:
            errmsg = ''
            self.errtx06.configure(text='')
        self.error_msgprint += errmsg
        
        #####################################################################
        #####################################################################
        
        # 7. Check for Orbit Inclination (degrees)
        
        self.var07.set(inps['orb_i'])
        if type(inps['orb_i']) != float or abs(inps['orb_i']) > 180.0:
            errmsg = 'Error! Inclination must be a float between +/- 180! \n'
            self.errtx07.configure(text='!')
        else:
            errmsg = ''
            self.errtx07.configure(text='')
        self.error_msgprint += errmsg
        
        #####################################################################
        #####################################################################
        
        # 8. Check for Orbit Right Asc. Node (degrees)
        
        self.var08.set(inps['orb_R'])
        if type(inps['orb_R']) != float or abs(inps['orb_R']) > 180.0:
            errmsg = 'Error! RAAN must be a float between +/- 180! \n'
            self.errtx08.configure(text='!')
        else:
            errmsg = ''
            self.errtx08.configure(text='')
        self.error_msgprint += errmsg
        
        #####################################################################
        #####################################################################
        
        # 9. Check for Orbit Arg. Perigee (degrees)
        
        self.var09.set(inps['orb_w'])
        if type(inps['orb_w']) != float or abs(inps['orb_w']) > 180.0:
            errmsg = 'Error! Perigee Arg must be a float between +/- 180! \n'
            self.errtx09.configure(text='!')
        else:
            errmsg = ''
            self.errtx09.configure(text='')
        self.error_msgprint += errmsg
        
        
        #####################################################################
        #####################################################################
        
        # 10. Check for Orbit Mean Anomaly (degrees)
        
        self.var10.set(inps['orb_m'])
        if type(inps['orb_m']) != float or abs(inps['orb_m']) > 180.0:
            errmsg = 'Error! Mean anomaly must be a float between +/- 180! \n'
            self.errtx10.configure(text='!')
        else:
            errmsg = ''
            self.errtx10.configure(text='')
        self.error_msgprint += errmsg
        
        #####################################################################
        #####################################################################
        
        # 11. Check for Maintenance Tolerance Band (km)
        
        self.var11.set(inps['orbm_tolr'])
        if type(inps['orbm_tolr']) == float and inps['orbm_tolr'] > 0.0:
            errmsg = ''
            self.errtx11.configure(text='')
        else:
            errmsg = 'Error! Tolerance must be a positive float (km)! \n'
            self.errtx11.configure(text='!')
        self.error_msgprint += errmsg
        
        #####################################################################
        #####################################################################
        
        # 12. Check if the orbit maintenance should be performed for FRO.
        
        if inps['orbm_fro_flag'] == 'True':
            errmsg = ''
            self.entry12.select()
            self.errtx12.configure(text='')
        elif inps['orbm_fro_flag'] == 'False':
            errmsg = ''
            self.entry12.deselect()
            self.errtx12.configure(text='')
        else:
            errmsg = 'Error! Repeat orbit option must be 1 or 0 (T/F)! \n'
            self.entry12.deselect()
            self.errtx12.configure(text='!')
        self.error_msgprint += errmsg
        
        #####################################################################
        #####################################################################
        
        # 13. Check for the spacecraft wet mass (kg)
        
        self.var13.set(inps['sc_mass'])
        if type(inps['sc_mass']) == float and inps['sc_mass'] > 0.0:
            errmsg = ''
            self.errtx13.configure(text='')
        else:
            errmsg = 'Error! Spacecraft mass must be a positive float! \n'
            self.errtx13.configure(text='!')
        self.error_msgprint += errmsg
        
        #####################################################################
        #####################################################################
        
        # 14. Check for the user input Isp x-axis minimum (s)
        
        self.var14.set(inps['isp_min'])
        if type(inps['isp_min']) == float and inps['isp_min'] > 0.0:
            errmsg = ''
            self.errtx14.configure(text='')
        else:
            errmsg = 'Error! Axis scale (min) must be a positive float! \n'
            self.errtx14.configure(text='!')
        self.error_msgprint += errmsg
        
        #####################################################################
        #####################################################################
        
        # 15. Check for the user input Isp x-axis maximum (s)
        
        self.var15.set(inps['isp_max'])
        if type(inps['isp_max']) == float and inps['isp_max'] > 0.0:
            if inps['isp_max'] > inps['isp_min']:
                errmsg = ''
                self.errtx15.configure(text='')
            else:
                errmsg = 'Error! ISP axis max must be greater than min! \n'
                self.errtx15.configure(text='!')
        else:
            errmsg = 'Error! Axis scale (max) must be a positive float! \n'
            self.errtx15.configure(text='!')
        self.error_msgprint += errmsg
        
        #####################################################################
        #####################################################################
        
        # Finally, display an error textbox if there are any error messages.
        
        if len(self.error_msgprint) > 0:
            tk.messagebox.showerror("Error with Configuration File!",
                                    self.error_msgprint)
            self.error_msgprint = '' # Reset error message
        
        return None
    
    #########################################################################
    #########################################################################
    ###                                                                   ###
    ###     Method for writing the entries into the config.txt file.      ###
    ###                                                                   ###
    #########################################################################
    #########################################################################
    
    def cfg_W(self):
        
        # First, ask the user if he/she wishes to proceed.
        cfg_W_msg = 'Save GUI parameters into the "config.txt" file? \n'
        cfg_W_msg += '(This will overwrite parameters in "config.txt"!)'
        cfg_W_ask = tk.messagebox.askyesno('Save Config', cfg_W_msg)
        if cfg_W_ask == False:
            return None
        
        # Else, continue with saving the configurations.
        cwd = dirname(dirname(abspath(__file__))) # Current working directory
        iwd = cwd + '\config\config.txt' # Inputs files
        input_r = open(iwd,'r') # Open the config.txt file
        record = [] # Array to record the strings
        
        # Variables to be recorded based on tkinter entries.
        var_arr = [self.var00, self.var01, self.var02, self.var03,
                   self.var04, self.var05, self.var06, self.var07,
                   self.var08, self.var09, self.var10, self.var11,
                   self.var12, self.var13, self.var14, self.var15]
        
        # Key values (to be referred).
        key_arr = ['orbsim', 'tstart', 'tfinal', 'sc_Cd', 'sc_Ad',
                   'orb_a', 'orb_e', 'orb_i', 'orb_R', 'orb_w', 'orb_m',
                   'orbm_tolr', 'orbm_fro_flag','sc_mass', 
                   'isp_min', 'isp_max']
        
        t_f_arr = ['orbm_fro_flag']
        
        # Prepare a dictionary to convert month strings into integers.
        months_dict = {'Jan':1, 'Feb':2, 'Mar':3, 'Apr':4,
                       'May':5, 'Jun':6, 'Jul':7, 'Aug':8,
                       'Sep':9, 'Oct':10,'Nov':11,'Dec':12}
        
        #####################################################################
        #####################################################################
        
        # 0. First, check for the orbit simulation program.
        
        try:
            _v00 = self.var00.get() # Exception raised if entry is erroneous
            if _v00 == 1 or _v00 == 2 or _v00 == 3:
                errmsg = ''
                self.errtx00.configure(text='')
            else:
                errmsg = 'Error! Expected a valid orbit mode (integer)! \n'
                self.errtx00.configure(text='!')
        except:
            _v00 = 1
            errmsg = 'Error! Unreadable orbit mode (integer)! \n'
            self.errtx01.configure(text='!')
            self.entry00a.select()
            self.entry00b.deselect()
            self.entry00c.deselect()
        finally:
            self.error_msgprint += errmsg
        
        #####################################################################
        #####################################################################
        
        # 1. Check for the first epoch string.
        
        try:
            _v01 = self.var01.get() # Exception raised if entry is erroneous
            
            # Check if the date and time are separated by hyphens
            if _v01.count('-') == 3:
                inp_v01s = _v01.split('-')
                
                # Check if it can be converted into a datetime object.
                try:
                    time_str1 = inp_v01s[3].split(':')
                    inp_v01d = datetime.datetime(int(inp_v01s[2]),
                                                 int(months_dict[inp_v01s[1]]),
                                                 int(inp_v01s[0]),
                                                 int(time_str1[0]),
                                                 int(time_str1[1]),
                                                 int(float(time_str1[2])))
                    errmsg = ''
                    self.errtx01.configure(text='')
                
                # If not, throw an exception and add it to the error log.
                except:
                    errmsg = 'Error! Invalid initial epoch parameters! \n'
                    self.errtx01.configure(text='!')
                
            # Else, raise an error message.
            else:
                errmsg = 'Error! Invalid initial epoch parameters! \n'
                self.errtx01.configure(text='!')
        except:
            errmsg = 'Error! Invalid initial epoch parameters! \n'
            self.errtx01.configure(text='!')
        finally:
            self.error_msgprint += errmsg
        
        #####################################################################
        #####################################################################
        
        # 2. Check for the final epoch string
        
        try:
            _v02 = self.var02.get() # Exception raised if entry is erroneous
            
            # Check if the date and time are separated by hyphens
            if _v02.count('-') == 3:
                inp_v02s = _v02.split('-')
                
                # Check if it can be converted into a datetime object.
                try:
                    time_str2 = inp_v02s[3].split(':')
                    inp_v02d = datetime.datetime(int(inp_v02s[2]),
                                                 int(months_dict[inp_v02s[1]]),
                                                 int(inp_v02s[0]),
                                                 int(time_str2[0]),
                                                 int(time_str2[1]),
                                                 int(float(time_str2[2])))
                    
                    # Check if the final epoch is after the initial epoch.
                    if inp_v02d <= inp_v01d:
                        errmsg = 'Error! The epoch final is before start! \n'
                        self.errtx02.configure(text='!')
                    else:
                        errmsg = ''
                        self.errtx02.configure(text='')
                
                # If not, throw an exception and add it to the error log.
                except:
                    errmsg = 'Error! Invalid final epoch parameters! \n'
                    self.errtx02.configure(text='!')
                
            # Else, raise an error message.
            else:
                errmsg = 'Error! Invalid final epoch parameters! \n'
                self.errtx02.configure(text='!')
            
        except:
            errmsg = 'Error! Invalid final epoch parameters! \n'
            self.errtx02.configure(text='!')
        finally:
            self.error_msgprint += errmsg
        
        #####################################################################
        #####################################################################
        
        # 3. Check for Spacecraft Atmospheric Drag Coefficient (Cd)
        
        try:
            _v03 = self.var03.get() # Exception raised if entry is erroneous
            if _v03 > 0.0:
                errmsg = ''
                self.errtx03.configure(text='')
            else:
                errmsg = 'Error! Drag coefficient must be a positive float! \n'
                self.errtx03.configure(text='!')
        except:
            errmsg = 'Error! Drag coefficient must be a positive float! \n'
            self.errtx03.configure(text='!')
        finally:
            self.error_msgprint += errmsg
        
        #####################################################################
        #####################################################################
        
        # 4. Check for Spacecraft Atmospheric Drag Surface Area (m^2)
        
        try:
            _v04 = self.var04.get() # Exception raised if entry is erroneous
            if _v04 > 0.0:
                errmsg = ''
                self.errtx04.configure(text='')
            else:
                errmsg = 'Error! Drag area must be a positive float! \n'
                self.errtx04.configure(text='!')
        except:
            errmsg = 'Error! Drag area must be a positive float! \n'
            self.errtx04.configure(text='!')
        finally:
            self.error_msgprint += errmsg
        
        #####################################################################
        #####################################################################
        
        # 5. Check for Orbit Semi-Major Axis (km)
        
        try:
            _v05 = self.var05.get() # Exception raised if entry is erroneous
            if type(_v05) != float:
                errmsg = 'Error! Semi-major axis must be a float! \n'
                self.errtx05.configure(text='!')
            elif _v05 < 6378.14:
                errmsg = 'Error! Semi-major axis below Earth surface! \n'
                self.errtx05.configure(text='!')
            elif _v05 > 385000.0:
                errmsg = 'Error! Semi-major axis beyond Earth orbit! \n'
                self.errtx05.configure(text='!')
            else:
                errmsg = ''
                self.errtx05.configure(text='')
        except:
            errmsg = 'Error! Semi-major axis must be a float! \n'
            self.errtx05.configure(text='!')
        finally:
            self.error_msgprint += errmsg
        
        #####################################################################
        #####################################################################
        
        # 6. Check for Orbit Eccentricity (no units)
        
        try:
            _v06 = self.var06.get() # Exception raised if entry is erroneous
            if type(_v06) != float:
                errmsg = 'Error! Eccentricity must be a float! \n'
                self.errtx06.configure(text='!')
            elif _v06 < 0:
                errmsg = 'Error! Eccentricity cannot be < 0! \n'
                self.errtx06.configure(text='!')
            elif _v06 >= 1.0:
                errmsg = 'Error! Eccentricity cannot be >= 1! \n'
                self.errtx06.configure(text='!')
            elif type(_v06) == float:
                if ( ( 1 - _v06 ) * _v05 ) < 6378.14:
                    errmsg = 'Error! Perigee altitude below Earth surface! \n'
                    self.errtx06.configure(text='!')
                else:
                    errmsg = ''
                    self.errtx06.configure(text='')
            else:
                errmsg = ''
                self.errtx06.configure(text='')
        except:
            errmsg = 'Error! Eccentricity must be a float! \n'
            self.errtx06.configure(text='!')
        finally:
            self.error_msgprint += errmsg
        
        #####################################################################
        #####################################################################
        
        # 7. Check for Orbit Inclination (degrees)
        
        try:
            _v07 = self.var07.get() # Exception raised if entry is erroneous
            if type(_v07) != float or abs(_v07) > 180.0:
                errmsg = 'Error! Inclination must be float between +/- 180! \n'
                self.errtx07.configure(text='!')
            else:
                errmsg = ''
                self.errtx07.configure(text='')
        except:
            errmsg = 'Error! Inclination must be float between +/- 180! \n'
            self.errtx07.configure(text='!')
        finally:
            self.error_msgprint += errmsg
        
        #####################################################################
        #####################################################################
        
        # 8. Check for Orbit Right Asc. Node (degrees)
        
        try:
            _v08 = self.var08.get() # Exception raised if entry is erroneous
            if type(_v08) != float or abs(_v08) > 180.0:
                errmsg = 'Error! RAAN must be a float between +/- 180! \n'
                self.errtx08.configure(text='!')
            else:
                errmsg = ''
                self.errtx08.configure(text='')
        except:
            errmsg = 'Error! RAAN must be a float between +/- 180! \n'
            self.errtx08.configure(text='!')
        finally:
            self.error_msgprint += errmsg
        
        #####################################################################
        #####################################################################
        
        # 9. Check for Orbit Arg. Perigee (degrees)
        
        try:
            _v09 = self.var09.get() # Exception raised if entry is erroneous
            if type(_v09) != float or abs(_v09) > 180.0:
                errmsg = 'Error! Perigee Arg must be float between +/- 180! \n'
                self.errtx09.configure(text='!')
            else:
                errmsg = ''
                self.errtx09.configure(text='')
        except:
            errmsg = 'Error! Perigee Arg must be float between +/- 180! \n'
            self.errtx09.configure(text='!')
        finally:
            self.error_msgprint += errmsg
        
        #####################################################################
        #####################################################################
        
        # 10. Check for Orbit Mean Anomaly (degrees)
        
        try:
            _v10 = self.var10.get() # Exception raised if entry is erroneous
            if type(_v10) != float or abs(_v10) > 180.0:
                errmsg = 'Error! Mean Anomaly must be between +/- 180.0! \n'
                self.errtx10.configure(text='!')
            else:
                errmsg = ''
                self.errtx10.configure(text='')
        except:
            errmsg = 'Error! Mean Anomaly must be between +/- 180.0! \n'
            self.errtx10.configure(text='!')
        finally:
            self.error_msgprint += errmsg
        
        #####################################################################
        #####################################################################
        
        # 11. Check for Maintenance Tolerance Band (km)
        
        try:
            _v11 = self.var11.get() # Exception raised if entry is erroneous
            if type(_v11) == float and _v11 > 0.0:
                errmsg = ''
                self.errtx11.configure(text='')
            else:
                errmsg = 'Error! Tolerance must be a positive float (km)! \n'
                self.errtx11.configure(text='!')
        except:
            errmsg = 'Error! Tolerance must be a positive float (km)! \n'
            self.errtx11.configure(text='!')
        finally:
            self.error_msgprint += errmsg
        
        #####################################################################
        #####################################################################
        
        # 12. Check if the orbit maintenance should be performed for FRO.
        
        try:
            _v12 = self.var12.get() # Exception raised if entry is erroneous
            if _v12 == 1 or _v12 == 0:
                errmsg = ''
                self.errtx12.configure(text='')
            else:
                _v12 = 0 # Default
                errmsg = 'Error! Invalid option for repeat orbit option! \n'
                self.errtx12.configure(text='!')
                self.var12.set(_v12)
                self.entry12.deselect()
        except:
            _v12 = 0 # Default
            errmsg = 'Error! Invalid option for repeat orbit option! \n'
            self.errtx12.configure(text='!')
            self.var12.set(_v12)
            self.entry12.deselect()
        finally:
            self.error_msgprint += errmsg
        
        #####################################################################
        #####################################################################
        
        # 13. Check for the spacecraft wet mass (kg)
        
        try:
            _v13 = self.var13.get() # Exception raised if entry is erroneous
            if _v13 > 0.0:
                errmsg = ''
                self.errtx13.configure(text='')
            else:
                errmsg = 'Error! Spacecraft mass must be a positive float! \n'
                self.errtx13.configure(text='!')
        except:
            errmsg = 'Error! Spacecraft mass must be a positive float! \n'
            self.errtx13.configure(text='!')
        finally:
            self.error_msgprint += errmsg
        
        #####################################################################
        #####################################################################
        
        # 14. Check for the user input Isp x-axis minimum (s)
        
        try:
            _v14 = self.var14.get() # Exception raised if entry is erroneous
            if _v14 > 0.0:
                errmsg = ''
                self.errtx14.configure(text='')
            else:
                errmsg = 'Error! Axis scale (min) must be a positive float! \n'
                self.errtx14.configure(text='!')
        except:
            errmsg = 'Error! Axis scale (min) must be a positive float! \n'
            self.errtx14.configure(text='!')
        finally:
            self.error_msgprint += errmsg
        
        #####################################################################
        #####################################################################
        
        # 15. Check for the user input Isp x-axis maximum (s)
        
        try:
            _v15 = self.var15.get() # Exception raised if entry is erroneous
            if _v15 > 0.0:
                if _v15 > _v14:
                    errmsg = ''
                    self.errtx15.configure(text='')
                else:
                    errmsg = 'Error! ISP axis max must be greater than min! \n'
                    self.errtx15.configure(text='!')
            else:
                errmsg = 'Error! Axis scale (max) must be a positive float! \n'
                self.errtx15.configure(text='!')
        except:
            errmsg = 'Error! Axis scale (max) must be a positive float! \n'
            self.errtx15.configure(text='!')
        finally:
            self.error_msgprint += errmsg
        
        #####################################################################
        #####################################################################
                
        # Finally, display an error textbox if there are any error messages.
        if len(self.error_msgprint) > 0:
            tk.messagebox.showerror("Error with Inputs!",
                                    self.error_msgprint)
            self.error_msgprint = '' # Reset error message
            return False
        
        else:
            
            # Now we parse through the config.txt file.
            for line in input_r:
                
                if line[0] == 'I':
                    words = line.split() # Split string into list of words.
                    key   = words[1] # Get the key from config.txt
                    value = words[2] # Get the value from config.txt
                    
                    # Get the updated value based on the index from key list.
                    value_new = str(var_arr[ key_arr.index(key) ].get())
                    
                    if key in t_f_arr:
                        if value_new == '1': # This is from Tkinter
                            value_new = 'True'
                        if value_new == '0': # This is from Tkinter
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
        
        return True
    
    #########################################################################
    #########################################################################
    ###                                                                   ###
    ###    Clears all existing relative orbit plots in the ORBITM GUI.    ###
    ###                                                                   ###
    #########################################################################
    #########################################################################
    
    def clr(self):
        
        # First, ask the user if he/she wishes to proceed.
        clr_msg = 'Clear the plot in the GUI?'
        clr_ask = tk.messagebox.askyesno('Clear Plots', clr_msg)
        if clr_ask == False:
            return None
        
        # Else, continue with clearing the plots.
        self.Axis211.clear() # Clear propulsion sizing
        self.Axis223.clear() # Clear altitude profile
        self.Axis224.clear() # Clear SMA profile
        self.FigPlot.draw()  # Update
        
        return None
    
    #########################################################################
    #########################################################################
    ###                                                                   ###
    ###        Run ORBITM by clicking on the RUN button in the GUI.       ###
    ###                                                                   ###
    #########################################################################
    #########################################################################
    
    def run(self):
        
        # First, ask the user if he/she wishes to proceed.
        run_msg = 'Save Parameters and Run ORBITM?'
        run_ask = tk.messagebox.askyesno('Run ORBITM?', run_msg)
        if run_ask == False:
            return None
        
        ######################################################################
        ######################################################################
        
        # Try to save the inputs, get the inputs, and run ORBITM...
        
        try:
            
            # Write the GUI config and run.
            error_flag = self.cfg_W()
            
            # Initialize all the input parameters
            mode = self.var00.get() # Orbit simulation program choice.
            ts   = self.var01.get().replace('-',' ') # Start Epoch String
            tf   = self.var02.get().replace('-',' ') # Final Epoch String
            Cd   = self.var03.get() # Atmospheric Drag Coefficient (Cd)
            Ad   = self.var04.get() # Atmospheric Drag Surface Area (m^2)
            oa   = self.var05.get() # Orbit Semi-Major Axis (km)
            oe   = self.var06.get() # Orbit Eccentricity (no units)
            oi   = self.var07.get() # Orbit Inclination (degrees)
            oR   = self.var08.get() # Orbit Right Asc. Node (degrees)
            ow   = self.var09.get() # Orbit Arg. Perigee (degrees)
            oM   = self.var10.get() # Orbit Mean Anomaly (degrees)
            mt   = self.var11.get() # Maintenance Tolerance Band (km)
            mf   = self.var12.get() # Toggle for normal or FRO Computation
            kg   = self.var13.get() # Spacecraft wet mass (kg)
            imin = self.var14.get() # Isp minimum (s) for x-axis
            imax = self.var15.get() # Isp maximum (s) for x-axis
            
            ##################################################################
            ##################################################################
            
            # If using Sam's orbit maintenance offline program
            if mode == 1 and error_flag == True:
                
                # Outputs five objects, with N = total number of samples.
                # epoch -> 1xN list comprising time of simulation in seconds
                # alt   -> 1xN list comprising geodetic altitude of satellite
                # sma   -> 1xN list comprising mean semi-major axis values
                # dv    -> Float value for total Delta-V required
                # imp   -> Float value for total impulse required
                
                epoch, alt, sma, dv, imp = orbmrun( ts, tf, Cd, Ad,
                                                    oa, oe, oi, oR, ow, oM,
                                                    mt, mf, kg, imin, imax )
            
            ##################################################################
            ##################################################################
            
            # If using STK 10 or STK 11
            elif (mode == 2 or mode == 3) and error_flag == True:
                
                # Outputs five objects, with N = total number of samples.
                # epoch -> 1xN list comprising time of simulation in seconds
                # alt   -> 1xN list comprising geodetic altitude of satellite
                # sma   -> 1xN list comprising mean semi-major axis values
                # dv    -> Float value for total Delta-V required
                # imp   -> Float value for total impulse required
                
                epoch, alt, sma, dv, imp = orbmstk( mode, ts, tf, Cd, Ad,
                                                    oa, oe, oi, oR, ow, oM,
                                                    mt, mf, kg, imin, imax )
            
            ##################################################################
            ##################################################################
            
            # Else... something went wrong?
            else:
                
                print('Invalid orbit mode and/or error flag encountered! \n')
                return None
            
            ##################################################################
            ##################################################################
            
            # Plotting of altitudes
            self.Axis223.set_ylabel('Geodetic Altitude (km)')
            self.Axis223.set_xlabel('Time of Simulation (s)')
            self.Axis223.scatter( epoch, alt, s=4, alpha=0.3 )
            self.Axis223.grid(True)
            
            # Plotting of Kozai-Izsak mean semi-major axes
            self.Axis224.set_ylabel('Mean Semimajor Axis (km)')
            self.Axis224.set_xlabel('Time of Simulation (s)')
            self.Axis224.plot( epoch, sma )
            self.Axis224.grid(True)
            
            # Thruster sizing profile of Isp Against Mass
            Isp = np.linspace( imin, imax, 500 ) # Isp axis, in (s)
            Mf  = imp / ( Isp * 9.81 )
            self.Axis211.set_ylabel('Mass of Fuel Required (kg)')
            self.Axis211.set_xlabel('Specific Impulse (s)')
            self.Axis211.plot(Isp, Mf)
            self.Axis211.grid(True)
            
            # Read the "thrusters.txt" text file, and compare specs against
            # the mission propulsion requirements for orbit maintenance.
            
            try:
                thr_file = open("thrusters.txt","r")
            except:
                
                # Otherwise, generate the file
                thr_file = open("thrusters.txt","w")
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
                thr_file = open("thrusters.txt","r")
            
            # Prepare arrays for bar chart plots
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
            bwidth = ( imin - imax ) / 50
            
            # Plot the thruster performances.
            barchart = self.Axis211.bar(thr_isp_s,
                                        thr_fuelm,
                                        width = bwidth,
                                        color='green')
            
            # Then, we label each thruster accordingly.
            barcount = 0
            for rect in barchart:
                bartext = thr_compn[barcount] + '\n'
                bartext = bartext + thr_model[barcount] + '\n'
                bartext = bartext + thr_force[barcount] + 'N'
                self.Axis211.text(rect.get_x() + rect.get_width()/2.0,
                                  rect.get_height(),
                                  bartext,
                                  ha='center', va='bottom')
                barcount += 1
            
            # Finally, update the plot.
            self.FigPlot.draw()
            
            # Print the output impulse and DV requirement message of ORBITM.
            imp_str = "The total impulse needed is "
            imp_str = imp_str + str(imp) + "kg m/s \n"
            dV_str  = "The total Delta-V (m/s) needed is "
            dV_str  = dV_str + str(dv) + " \n"
            tk.messagebox.showinfo(title="Complete!", message=imp_str+dV_str)
            
            ##################################################################
            ##################################################################
            
        except Exception as excpt:
            print('Error in running!')
            print(excpt)
            pass
        
        return None
    
