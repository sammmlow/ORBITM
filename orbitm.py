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
##    Runs the GUI that interfaces with the ORBITM source code stack.        ##
##                                                                           ##
##    Written by Samuel Y. W. Low.                                           ##
##    First created 27-Nov-2020 12:00 PM (+8 GMT)                            ##
##    Last modified 19-Sep-2021 22:27 PM (-7 GMT).                           ##
##                                                                           ##
###############################################################################
###############################################################################

# Import our GUI libraries.
import tkinter

# Import local libraries.
from source import orbmgui

# Initialise the GUI.
root = tkinter.Tk()
root_gui = orbmgui.RunGUI( root )
root.mainloop()