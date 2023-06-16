# -*- coding: utf-8 -*-
"""
Created on Tue Dec 13 17:52:21 2022

Alternative script to obtain particular plots from the automation of GCN circulars
from a particular FITS file. This Script is prepared to run over all FITS included
in the Data_FITS folder, or alternatively, it can be launched over one particular
FITS file if we give the name.

"""

import subprocess
import os

# =============================================================================
# Parameters of the Scripts
# =============================================================================

'''
Scripts is a variable that includes what kind of plots we want and also if we
want to receive an email or not. It is just a str, that should include what we 
want from the following 4 things

 Limit_Flux  FOV_GW  CL_Coverage  email
 
for instance, if we want just the FOV_GW plot, and we want to receive an email
we use Scripts = 'FOV_GW email' as you can see, everything between the comas and
separated by blanks. if we want all, then Scripts = 'Limit_Flux FOV_GW CL_Coverage email'
'''


Scripts = 'FOV_GW'


'''
GCN_mode should be chosen 'no' for running over /Data_FITS
if we want to run over a specific file we can choose GCN_mode = "yes", then
we should include the GCN_ID, which can simply be chosen *** if this information
is not relevant. With this method the three plots are produced and the email
is sent. If we don´t want all plots or the circular email, we can run the 
GCN_mode = 'no'  mode with just the file we want in the folder /Data_FITS
'''

GCN_mode = 'no'
GCN_ID   = "GW190517"
FILENAME = "bayestar.fits.gz,0"

'''
The Separate_Lobes function in Contour.py compares the average distance between 
points in a lobe and the distance of a new point with respect to the last one, 
when following the ordered points array. So, parameter_lobes is just the number 
of times greater than the distance to a new point should be with respect to the
average in order to split and create a new lobe.
 
In GCN_mode it is automatically set to 8, which on average works fine, but for
optimal performance can be tuned here
'''

parameter_lobes = 8

'''
This is the resolution from the CL_Coverage plot, reducing this parameter 
will lead to a better resolution on the plot at the cost of a longer run time.
This parameter is set to 1. automatically in GCN_mode = "yes"
'''

#SCAN_RESOLUTION = 1./2.  # deg     # Scan every 2 minutes
SCAN_RESOLUTION  = 30.    # deg     # Scan every 4 minutes 
#SCAN_RESOLUTION = 15.    # deg     # Scan every hour 15
#SCAN_RESOLUTION = 360.   # deg     # 1 point only at the time of the event


'''
T parameter is a 'yes' or 'no' parameter which asks if we want, in the FOV plot,
to have AR - t_source (AR minus the sidereal time of the source, for this T = 'yes') 
or just AR (for this T = 'no'). In GCN_mode it is set to 'no'

Mollweide works similarly, if we want Mollweide projection 'yes', otherwise, 'no'
In GCN_mode it is set to 'yes'.

'''

T = 'no'

Mollweide = 'yes'

# =============================================================================
# Run the Scripts
# =============================================================================

parameters = str(parameter_lobes) + ' ' + str(SCAN_RESOLUTION) + ' ' + T + ' ' + Mollweide

currentdir = os.getcwd()

FILENAME = currentdir + "/Data_FITS/" + FILENAME

'''
In the case of GCN_mode = 'no' we basically pass the subprocess a str called parameters 
which contains all our above preferences like the Scripts variable. Of course
we pass the Scripts variable too.

In the case of GCN_mode = 'yes' we pase the Filename and the GCN_ID
'''


if GCN_mode == "no":

    process = subprocess.run(["python3", "Contour.py",GCN_mode, Scripts, parameters]) 

else:
    
    process = subprocess.run(["python3", "Contour.py",GCN_mode, FILENAME, GCN_ID])