# -*- coding: utf-8 -*-
"""
Created on Tue Dec 13 17:52:21 2022

@author: yagol
"""

import subprocess
import os

# =============================================================================
# Parameters of the Scripts
# =============================================================================

# Here we choose the plots we want to do and also if we want to send an email or
# not. The variable Scripts should contain Limit_Flux, FOV_GW, CL_Coverage or
#  email, and its just one str separated by blanks.

Scripts = 'FOV_GW CL_Coverage Limit_Flux email'

# GCN_mode should be chosen no for running over /Data_FITS

GCN_mode = 'no'

# if we want to run over a specific file we can choose GCN_mode = "yes", then
# we should include the GCN_ID, which can simply chose *** if this information
# is not relevant. With this method the three plots are produced and the email
# is sent. If we don´t want all plots or the circular email, we can run the 
# GCN_mode = 'no'  mode with just the file we want in the folder /Data_FITS

GCN_ID   = "***"
FILENAME = "GW150914_skymap.fits.gz"

'''
The Separate_Lobes function in Contour.py compares the average distance between 
points in a lobeand the distance of a new point with respect to the last one, 
when following the ordered points array. So parameter_lobes is just the number 
of times greater than the distance to a new point should be with respect to the
average in order to split and create a new lobe.
 
In GCN_mode it is automatically set to 8, which on average works fine
'''

parameter_lobes = 8

'''
This is the resolution from the CL_Coverage plot, reducing this parameter 
will lead to a better resolution on the plot at the cost of a longer run time.
This parameter is set to 1. automatically in GCN_mode = "yes"
'''

#SCAN_RESOLUTION = 1./2.  # deg     # Scan every 2 minutes
#SCAN_RESOLUTION  = 1.    # deg     # Scan every 4 minutes 
SCAN_RESOLUTION = 15.    # deg     # Scan every hour
#SCAN_RESOLUTION = 360.   # deg     # 1 point only at the time of the event

# =============================================================================
# Run the Scripts
# =============================================================================

currentdir = os.getcwd()

FILENAME = currentdir + "/Data_FITS/" + FILENAME


if GCN_mode == "no":

    process = subprocess.run(["python3", "Contour.py",GCN_mode, Scripts, str(parameter_lobes), str(SCAN_RESOLUTION)]) 

else:
    
    process = subprocess.run(["python3", "Contour.py",GCN_mode, FILENAME, GCN_ID])