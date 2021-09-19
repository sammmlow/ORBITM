###############################################################################
###############################################################################
##                                                                           ##
##     _____ ___  ___  ___  _____      __  __                                ##
##    |  _  | _ \| _ \|_ _||_   _|    |  \/  |                               ##
##    | |_| |   <| _ < | |   | |   _  | \  / |                               ##
##    |_____|_|\_|___/|___|  |_|  |_| |_|\/|_|                               ##
##                                                     v 1.0                 ##
##                                                                           ##
##    FILE DESCRIPTION:                                                      ##
##                                                                           ##
##    Runs the GUI that interfaces with the ORBITM source code stack.        ##
##                                                                           ##
##    Written by Samuel Y. W. Low.                                           ##
##    Last modified 26-10-2020.                                              ##
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