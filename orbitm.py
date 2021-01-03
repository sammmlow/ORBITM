###############################################################################
###############################################################################
##                                                                           ##
##     _____ ___  ____  ___  _____       ______                              ##
##    |  _  | _ \|  _ \|_ _||_   _|     |      |                             ##
##    | |_| |   <|  _ < | |   | |       | \  / |  _                          ##
##    |_____|_|\_|____/|___|  |_|       |_|\/|_| |_|                         ##
##                                                     v 0.1                 ##
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
from codes import orbitm_tkinter_ui

# Initialise the GUI.
root = tkinter.Tk()
root_gui = orbitm_tkinter_ui.run_gui( root )
root.mainloop()