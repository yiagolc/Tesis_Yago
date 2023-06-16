# -*- coding: utf-8 -*-
"""
This Script has the aim of plotting the percentage of the 90% CL region of a GW
in coincidence with the FOV of the Pierre Auger observatory (DGH, DGL, ES and total)
during an entire day after the detection time, and also the probability of the
position in the sky of the GW contained in this FOV as a function of the time.
The Script just needs the FITS_file for working.

"""

# =============================================================================
# Imports
# =============================================================================

import sys
import healpy as hp
import numpy as np
import astropy.coordinates
import astropy.time
import astropy.units as u
import os

import matplotlib.pyplot as plt
import matplotlib as mpl
import _pickle as pcl 

import SendTheMail as SDM
from string import Template

import subprocess

import time

# =============================================================================
# Manual parameters to the program in the case of using it independently
# =============================================================================

#print(os.getpid())

#scriptStartTime = time.time()
'''

SUFFIX = '.fits.gz'
PREFIX = os.getcwd() + "/Data_FITS/"

DECREASE_RESOLUTION = False
DECREASE_NSIDE      = 512 # corresponds to 0.0131139632064 deg^2 per pixel ~ (0.115 deg)^2

#SCAN_RESOLUTION = 1./2.  # deg     # Scan every 2 minutes
SCAN_RESOLUTION  = 1.     # deg     # Scan every 4 minutes (1)
#SCAN_RESOLUTION = 15.    # deg     # Scan every hour
#SCAN_RESOLUTION = 360.  # deg     # 1 point only at the time of the event


#PREFIX = 'Auger/alerts/'
SUFFIX = '.fits.gz'
PREFIX = "Data_FITS/"
PICKLE_PREFIX = "Pickle/"


PICKLE_FILENAME = PICKLE_PREFIX + FILENAME


if PREFIX not in FILENAME:
  FILENAME = PREFIX + FILENAME
  
if SUFFIX not in FILENAME:
  SUFFIX='.fits' # assuming you used the unpacked fits file
  
'''
# =============================================================================
# Functions
# =============================================================================


def get_ConfReg_minimum_prob(prob,CL):

    '''
    prob : 1D array containing probability in each pixel
    CL : confidence level
    
    This function takes prob and sums all the elements from the largest to the 
    smallest until the sum is greater than CL, then returns the prob element 
    that made the sum greater than CL. Since it is expected that sum(prob) = 1 
    and CL < 1, then an error warning is given in case the probabilities sum
    do not exceed CL
    
    '''    
    
    tempsum = 0
    
    # np.sort(prob)[::-1] ordena las probabilidades de mayor a menor
    
    for i in np.sort(prob)[::-1]:
        tempsum+=i # np.cumsum()
        if tempsum > CL:
            return i
    print('THIS SHOULD NOT HAPPEN! STILL, THE PROBABILITIES DID NOT SUM UP TO',CL)
    sys.exit()
    
def print_and_write(outputText,outfile):

    outfile.write(outputText+'\n')
    
# =============================================================================
# Constant definitions
# =============================================================================

SKY = 4 * 180**2 / np.pi
AUGERLAT, AUGERLONG, AUGERHEIGHT = -35.20666735*u.deg, -69.315833*u.deg, 1400*u.m
CL  = 0.9
regions = [ 'DGL', 'DGH', 'ES' ]
colors  =  [ 'g',   'b',   'r' ]

observatory = astropy.coordinates.EarthLocation( lat=AUGERLAT, lon=AUGERLONG, height=AUGERHEIGHT )

# =============================================================================
#%% Directory (only useful if we want to run over /FITS_files)
# =============================================================================


'''
script_dir = os.getcwd()

path = os.path.join(script_dir, "Data_FITS")

#text_in_file = "bayestar.fits.gz,0"
#text_in_file = "GW170814_skymap.fits.gz"
text_in_file = ".gz"

list_path_file = []   
list_file      = [] 

for dirpath, dirnames, filenames in os.walk(path):
    for filename in [f for f in filenames if text_in_file in f]:
        
        path_with_filename_e=os.path.join(dirpath, filename)
        list_path_file.append(path_with_filename_e)
        list_file.append(filename)
    list_path_file.sort() 
    list_file.sort()  

'''



# =============================================================================
#%% Start by reading input data
# =============================================================================

if(len(sys.argv) < 2) :
    print('You have not include a fits file in the input')
    sys.exit()

else:
    
    GCN_mode = sys.argv[1]
    FILENAME = sys.argv[2]
    SCAN_RESOLUTION = float(sys.argv[4])
    
    print(sys.argv)
    
    if GCN_mode == "yes":

        GCN_ID  = sys.argv[3]
        Scripts = "FOV_GW CL_Coverage Limit_Flux email"
        
    else:
        Scripts = sys.argv[3]
        
    
if "CL_Coverage" not in Scripts.split():
    
         # If CL_Coverage is not included between the plots that we want to do
         # then we don't want to run this Script, then we continue
        
         
         # We need the GWname
         
         prob, header = hp.read_map( FILENAME, h=True)
         
         time_index = -1
         for i in range(len(header)):
          try:
            header[i].index('MJD-OBS') 
            time_index = i
            
          except:
              pass
            
        # Este bucle simplemente encuentra la posicion en la cabecera en la que está 
        # MJD-OBS y OBJECT y la asigna a time_index y object_index, respectivamente
        
         if time_index > -1 : 
             GWtime = astropy.time.Time( header[time_index][1], format='mjd' )
            
             GWname = "GW" + str(GWtime.datetime)[2:4] + str(GWtime.datetime)[5:7] + str(GWtime.datetime)[8:10] + "_" + str(GWtime.datetime)[11:13]+str(GWtime.datetime)[14:16]
            
             if GCN_mode == 'yes':
            
                GWname = GCN_ID
         else:               
             print('time_index invalid:',time_index,'could not find the MJD from the header file')
             sys.exit()
             
         # If we want to execute Limit_Flux, then we do it here 
         
         if GCN_mode == 'no':
             GCN_ID = GWname
            
         if "Limit_Flux" in Scripts:
             print('Executing Limit_Flux.py for', GWname)
             process = subprocess.run([ 'python3', 'Limit_Flux.py','no', FILENAME, Scripts, GCN_ID])
         
         # If we want the email, it is sent here, but in this case without the
         # Circular template, because we are not in the GCN_mode            
            
         print(Scripts.split())
            
         if "email" in Scripts.split():
             
                    # The main text
                    
                    MT = 'Resultados'
                    
                    print("email") # To notify that the email is going to be send
                    
                    attachments = []
        
                    # We send just the wanted plots
                    
                    if "FOV_GW" in Scripts.split():
                        attachments.append("Plot_FOV_Contours/fov_Auger_{0}_mollweide.png".format(GWname))
                    if "CL_Coverage" in Scripts.split():
                        attachments.append("Plot_CL_Coverage/GW_confidence_region_coverage" + GWname + ".png")
                    if "Limit_Flux" in Scripts.split():
                        attachments.append("Plots_Flux_limit/" + GWname + ".png")
        
                        
                    SDM.SendTheMail(SDM.sender, SDM.sendermail, SDM.receivers, SDM.receivermails, MT, attachments)
                    
        
         sys.exit()# if CL_Coverage is not in the Scripts str, we don´t want to run this Script    
        
                    
            
         
SUFFIX = '.fits.gz' # to take into account files .fits and .fits.gz
if SUFFIX not in FILENAME:
    SUFFIX = '.fits'
PREFIX = os.getcwd() + "/Data_FITS/"
    
#print(FILENAME[len(PREFIX):-len(SUFFIX)])
#print("Data_CL_Coverage/"+FILENAME[len(PREFIX):-len(SUFFIX)]+'_out.txt')

# This is a file where we are going to store data from this script

outfile = open("Data_CL_Coverage/"+FILENAME[len(PREFIX):-len(SUFFIX)]+'_out.txt','w')


# =============================================================================
# Read Sky map
# =============================================================================

# Download sky map. Not implemented now
#~ subprocess.check_call([    'curl', '-O', '--netrc',    'https://gracedb.ligo.org/apiweb/events/G298389/files/LIB.fits.gz,0'])

# prob is the unidimensional array containing the probability in each pixel

prob, header = hp.read_map( FILENAME, h=True)

# Map properties

# Total number of pixels
npix = len(prob)
print_and_write('Total number of pixels in map: ' + str( npix ), outfile)

# Area per pixel, el NSIDE es un parametro relacionado con la resolucion pero no tengo claro que es
area_per_pix, nside = SKY / npix, hp.npix2nside(npix)
print_and_write('Area per pixel (before eventual resolution decrease): '+str( area_per_pix )+' deg^2',outfile)

'''
# Decrease resolution and repeat
if DECREASE_RESOLUTION and ( nside > DECREASE_NSIDE ): 
  # si el NSIDE es mayor que el que queremos
  
  # Upgrade or degrade resolution of a map (or list of maps).
  # in degrading the resolution, ud_grade sets the value of the superpixel as 
  # the mean of the children pixels. power = -2 keeps the sum of the map 
  # invariant, which is useful given that it is a probabilistic map
  
  prob = hp.pixelfunc.ud_grade( prob, nside_out = DECREASE_NSIDE, power=-2 )
# The next lines are not calculated from skymap but "deterministically" calculated from given reduced resolution.
  nside = DECREASE_NSIDE
  npix  = hp.nside2npix(nside)
  area_per_pix = SKY / npix
  print_and_write('Area per pixel (after resolution decrease): '+str( area_per_pix )+' deg^2',outfile)
'''

# =============================================================================
#  Confidence Region of sky localization corresponding to CL
# =============================================================================

# First obtain the minimum value of prob array above the chosen CL
ConfReg_minimum_prob = get_ConfReg_minimum_prob(prob,CL)

# Using the get_ConfReg_minimun_prob this way, we obtain the maximun probability
# point value

Maximun = get_ConfReg_minimum_prob(prob, 0)

npix_ConfReg = np.sum( prob > ConfReg_minimum_prob )

print_and_write('Number of pixels in the confidence region: ' + str( npix_ConfReg ), outfile)
print_and_write('Solid angle of the confidence region:      ' + str( area_per_pix*npix_ConfReg )+' deg^2', outfile)
 
# Obtain all thetas, phis of the pixels in all the map
theta, phi = hp.pix2ang(nside, np.arange(npix))

# Convert thetas, phis of all the pixels to RA, Dec.
radecs = astropy.coordinates.SkyCoord( ra=phi*u.rad, dec= (0.5*np.pi - theta) * u.rad )

# get the declination of the maximun value (for the limit flux)

prob_dec = [(probi, (0.5*np.pi - thetai)*180/np.pi) for probi,thetai in zip(prob,theta)]

declination = [prob_deci[1] for prob_deci in prob_dec if prob_deci[0] == Maximun][0]


# =============================================================================
# Our observation and observatory properties (time, coordinates, event)
# =============================================================================

# bucle sobre el cabecero del fichero FITS

# Find the time (MJD = Modified Julian Date: A continuous measure in days since 
#midnight at the start of 17 November 1858. Based on UTC) in the header

time_index = -1
for i in range(len(header)):
  try:
    header[i].index('MJD-OBS') 
    time_index = i
    
  except:
      pass
  ''' In order to extract the name from the header (some fits doesn´t include it)
  try:
      header[i].index('OBJECT') 
      object_index = i
  except:
      pass
  '''
    
# Este bucle simplemente encuentra la posicion en la cabecera en la que está 
# MJD-OBS y OBJECT y la asigna a time_index y object_index, respectivamente

if time_index > -1 : 
    GWtime = astropy.time.Time( header[time_index][1], format='mjd' )
    
    #print(GWtime.datetime)
    
    # We construct the GWname from the date:
    
    GWname = "GW" + str(GWtime.datetime)[2:4] + str(GWtime.datetime)[5:7] + str(GWtime.datetime)[8:10] + "_" + str(GWtime.datetime)[11:13]+str(GWtime.datetime)[14:16]
    
    if GCN_mode == 'yes':
        
        GWname = GCN_ID
    
    #print(GWname)
else:               
    print('time_index invalid:',time_index,'could not find the MJD from the header file')
    sys.exit()
 
print("Executing CL_Coverage for", GWname)


'''   In order to extract the name from the header (some fits doesn´t include it)
if object_index > -1 : 
    GWname = header[object_index][1]
else:               
    print('object_index invalid:',object_index,'could not find the Object type from the header file')
    sys.exit()
'''



# =============================================================================
# Loop through the MJD during a day
# =============================================================================

# The number of iterations of the first loop will be taken into account
mjdindex = 0

# arrays containing the probability on the FOV and the fraction of the 90% CL
# region into the FOV of Auger
fov_prob = []
fov_frac = []


# Since the angular resolution is SCAN_RESOLUTION in degrees, the number of scans
# that we make is 360/SCAN_RESOLUTION + 1 (taking into account the 0 and the 360)


NSCANPOINTS = 360/SCAN_RESOLUTION+1 - (SCAN_RESOLUTION == 360) 
NSCANPOINTS = int(NSCANPOINTS)

# This last expression in the parentheses makes the argument of the following
# for-loop get only one element if the SCAN_RESOLUTION is 360, i.e. only the 
# actual time of the GW event is used

ifirst = 1

# ----- Start loop in MJD for 1 day

# The linspace makes equidistant time points between GWtime.mjd and the 
# following day GWtime.mjd+1

for mjd in np.linspace( GWtime.mjd, GWtime.mjd+1, NSCANPOINTS ):
    
    currtime = astropy.time.Time( mjd, format='mjd' )
    print_and_write( '\n\n'+str(currtime.datetime)+'\n' , outfile)
    
    # Observatory frame (unfortunately has to be redone for each MJD in the loop)
    
    frame = astropy.coordinates.AltAz( obstime=currtime, location=observatory )
    
    # radecs was the phis and thetas in RA and Dec from de FITS, we change
    # the coordinates system to the AltAz from the observatory with
    # (This is probably the longest taking step :/)
    
    altaz = radecs.transform_to(frame)
    
    #~ FoV definition - set of conditions fulfiled by each FoV
    # Arrays containing bool values telling if the points belong to
    # each FoV region
    
    fov_DGL = ( altaz.alt <= 30*u.deg ) & ( altaz.alt > 15*u.deg )
    fov_DGH = ( altaz.alt <= 15*u.deg ) & ( altaz.alt >  0*u.deg )
    fov_ES  = ( altaz.alt <=  0*u.deg ) & ( altaz.alt > -5*u.deg )
    
    all_fovs = [fov_DGL, fov_DGH, fov_ES]
    
    # sum_all_fovs = [np.sum(fovi) for fovi in all_fovs]
  
    # Do the following only once : it is intended to print things, not neccessary
    if ifirst == 1:

       #print('Fraction of sky in fov ES  [%]: ',sum_all_fovs[2]*100./len(fov_ES)) GCNListening
       #print('Fraction of sky in fov DGH [%]: ',sum_all_fovs[1]*100./len(fov_DGH))  GCNListening
       #print('Fraction of sky in fov DGL [%]: ',sum_all_fovs[0]*100./len(fov_DGL))  GCNListening
       #print(' ')  GCNListening
       
       ifirst = -1 
       
    
    '''
    Example to understand the following line
    
    a = np.array([1,2,3,4,5])
    b = np.array([0,4])
    
    then a[b] will be np.array([1,5])
    
    '''
    
    #~ Create 3 arrays containing the values of prob in the pixels that are in each FoV (ES, DGH, DGL) 
    # we are using here a mask to just select the wanted points
    
    all_prob_fovs = [prob[sub_fov] for sub_fov in all_fovs]
    
    #~ Loop through the FoV and get:
    # ---------------
    # Total probability (sum of the probabilities assigned to each pixel)
    # Total probability (sum of area of pixels divided by total area, i.e. 
    # all pixels have equal probability)
    # ---------------
    
    #~ Total probability of CL contour inside each FoV
    # We append empty lists to the lists defined out of the loop, and then we will
    # fill them
    fov_prob.append([])
    #~ Fraction of CL contour inside each FoV
    
    fov_frac.append([])
    
    # for loop trough the three FOV´s
    
    for i in range(len(all_prob_fovs)):
        
        # First we access the empty list and add the sum of probabilities
        # in the fov of each 
        
        fov_prob[mjdindex].append(np.sum(all_prob_fovs[i]))
        print_and_write('Total probability (sum of prob. pixels) for the event to be in the '+str(regions[i])+' FoV: '+str( fov_prob[mjdindex][i] ), outfile)
        
        # We do the same for the overlapping region
        
        fov_frac[mjdindex].append( float(np.sum(all_prob_fovs[i] > ConfReg_minimum_prob)) / npix_ConfReg )
        print_and_write('Fraction of 90 % confidence region in the '+str(regions[i])+' FoV:    '+str( fov_frac[mjdindex][i] ), outfile)
        print_and_write('', outfile)
        
    #~ Summary
    # We add the three probabilites to be in the DGH DGL and ES regions
    print_and_write('Total probability for the event to be in the total FoV: '+str( np.sum(fov_prob[mjdindex]) ), outfile)
    print_and_write('Fraction of 90 % confidence region in the total FoV:    '+str( np.sum(fov_frac[mjdindex]) ), outfile)
    
    mjdindex += 1

# Preserve the first overlap and probability

fov_prob_first = np.sum(fov_prob[0]) 
#print(fov_prob_first)
fov_frac_first = np.sum(fov_frac[0])

#print('------------------------------------------------------------------------------------------------------') GCN LISTENING
print_and_write('\n\nMaximum probability (sum of prob. pixels) for the event to be in the total FoV:    '+str( max(np.sum(fov_prob,axis=1)) )+' '+str( (np.sum(fov_prob,axis=1).argmax()/360.*SCAN_RESOLUTION*24)%24 )+' hours after the GW event.', outfile)
#print('------------------------------------------------------------------------------------------------------') GCN LISTENING
#print('Initial fraction of 90 % confidence region in the total FoV: ',fov_frac_first) GCN LISTENING

print_and_write('Maximum fraction of 90 % confidence region in the total FoV: '+str( max(np.sum(fov_frac,axis=1)) )+' '+str( (np.sum(fov_frac,axis=1).argmax()/360.*SCAN_RESOLUTION*24)%24 )+' hours after the GW event.\n', outfile)

outfile.close()


# =============================================================================
# Output files in pickle format
# =============================================================================

'''
for i in range(len(regions)):
  pcl.dump( np.array([np.linspace( 0, 24, NSCANPOINTS ), np.array(fov_prob)[:,i]]), open(PICKLE_FILENAME[:-len(SUFFIX)]+'_prob_vs_time_'+regions[i]+GWname+'.pickle','wb') )
  pcl.dump( np.array([np.linspace( 0, 24, NSCANPOINTS ), np.array(fov_frac)[:,i]]), open(PICKLE_FILENAME[:-len(SUFFIX)]+'_frac_vs_time_'+regions[i]+GWname+'.pickle'+GWname,'wb') )
'''

# =============================================================================
# Plots of:
# Sum of probabilities of localization map inside each field-of-view (ES, DGH, DGL) for each MJD time
# Fraction of localization map inside each field-of-view (ES, DGH, DGL) for each MJD time
# =============================================================================

mpl.rc('lines', linewidth=2.0)
    
xs = np.linspace( 0, 24, NSCANPOINTS ) # Per hour
fs1, fs2 = 16, 12  # Font sizes

# Array of two subplots
#You can specify nrows and ncols. sharex allows to share the x axis
fig, axarr = plt.subplots( 2, sharex=True ,figsize=(8,5) )

#---------------------------- Top subplot -----------------------------------

# Sum of probabilities of localization map inside each field-of-view (ES, DGH, DGL) 
# for each MJD time and for the sum of the 3 FoV

axarr[0].plot( xs  , np.sum(fov_prob,axis=1), color='black', label="sum" )
for i in range(len(regions)):
  axarr[0].plot( xs, np.array(fov_prob)[:,i], label=regions[i], color=colors[i])
  

axarr[0].set_ylabel( 'fov probability', fontsize=fs2 )
# axarr[0].set_ylims(0,1)
axarr[0].set_title( str(GWtime.datetime.year)+'-'+str(GWtime.datetime.month).zfill(2)+'-'+str(GWtime.datetime.day).zfill(2), fontsize=fs1, loc='right')
axarr[0].grid(color='gray', linestyle=':', linewidth=1)

plt.xlim( 0, 24 )
plt.xticks(np.arange(0, 25, 1.0))
plt.ylim( 0, 1 )

# Try to place legend so that it is not on top of lines: 
   # total prob at t=t_0 is not too high -> draw legend in upper left
if ( np.sum(fov_prob,axis=1)[0] < .5*np.max(np.sum(fov_prob,axis=1)) ): axarr[0].legend( loc = 'upper left' ) 
# first guess: move it a bit to the right
else:   axarr[0].legend( bbox_to_anchor=(0.3, .82, .1, .2) ) 

# -------------------------- Bottom subplot ----------------------------------
# Fraction of localization map inside each field-of-view (ES, DGH, DGL) for each MJD time
# and for the sum of the 3 FoV
axarr[1].plot( xs  , np.sum( fov_frac,axis=1), color='black', label="sum" )
for i in range(len(regions)):
  axarr[1].plot( xs, np.array(fov_frac)[:,i],label=regions[i],color=colors[i] )

axarr[1].set_ylabel( '90 % conf. region in fov',fontsize=fs2 )
# ~ axarr[1].set_xlabel('MJD')
axarr[1].set_xlabel( 'hours since GW event', fontsize=fs2 )
# axarr[1].set_ylims(0,1)
axarr[1].grid(color='gray', linestyle=':', linewidth=1)

plt.xlim( 0, 24 )
plt.xticks(np.arange(0, 25, 1.0))
maxValue = max( np.max(np.sum(fov_prob,axis=1)), np.max(np.sum(fov_frac,axis=1)) )
# plt.ylim( 0, min( maxValue-(maxValue%.1)+0.1, 1.) )
plt.ylim( 0, 1)

#print('----------------------FINISH--------------------------------')  GCN LISTENING
#print('approx. run time of the script', time.time()-scriptStartTime, 's') GCN LISTENING

plt.savefig("Plot_CL_Coverage/GW_confidence_region_coverage" + GWname + ".png")
# plt.savefig("CL_Coverage/GW_confidence_region_coverage" + nameGW + ".pdf")


# =============================================================================
# maximun times of coincidence calculation
# =============================================================================
    

dt = SCAN_RESOLUTION*60/15

# The added probability for the three channels

prob = np.sum(fov_prob,axis=1)

# We take the maximun position

maximo = np.max(prob)
posiciones = np.where(prob == maximo)[0]



# Once we have the maximun time, then we start two while loops moving to the
# right and to the left until find the time where we have decreased 70% of the
# maximun

# To the right
i = posiciones[0]
while prob[i] >= 0.95*maximo :
    i = i + 1
    if i == (len(prob) - 1):
        i = 0

t_fin = i*dt/60 # This is the right limit of the maximun region
i = posiciones[0]    

#To the left
while prob[i] >= 0.95*maximo :
    i = i - 1
    if i == 0:
        i = len(prob) - 1
    
    
t_in = i*dt/60 # This is the left limit of the maximun region


delta  = t_fin - t_in
# more than 1 hour and the duration information is included in the template
threshold = 1 

#print(t_in, '\n', t_fin)
#print(prob)
    
if delta > threshold:
    t_max = round(t_in, 1)
    maximum = str(round(maximo*100)) + '% for approximately '+ str(round(delta,2))+' hr'
else:
    t_max = round(posiciones[0]/60*dt,1)
    maximum = str(round(maximo*100)) + '%'
    
'''

The probability of the event to be coincident with the Field-of-View (FoV) where 
the SD of Auger is sensitive to UHE neutrinos (corresponding to inclined directions with 
respect to the vertical relative to the ground) was ${fov_prob_first} %, 
at the time T0 of the merger alert. The LVK 90% localization region
maximally overlapped with the Auger fov at ~ T0+${t_in} hours (i.e.
${t_in} hr after the merger) and remained until ~ T0+${t_fin} hours. 
'''
    
# =============================================================================
# We run the calculation of the flux limit or the email if neccessary
# =============================================================================

attachments = []

 
if "FOV_GW" in Scripts:
      attachments.append("Plot_FOV_Contours/fov_Auger_{0}_mollweide.png".format(GWname))
if "CL_Coverage" in Scripts:
      attachments.append("Plot_CL_Coverage/GW_confidence_region_coverage" + GWname + ".png")
if "Limit_Flux" in Scripts:
      attachments.append("Plots_Flux_limit/" + GWname + ".pdf")

if GCN_mode == "yes":
    
    # If we are in GCN mode we want to run the Limit flux and also capture the
    # fluence limit from the output in order to include it in the Circular
    
    print('Executing Limit_Flux.py for', GWname)
    Flux_limit = subprocess.run([ 'python3', 'Limit_Flux.py', 'yes', FILENAME, GCN_ID], capture_output=True, text=True)
    
    out_text = Flux_limit.stdout.split()
    
    try:
        limit_max   = round(float(out_text[0]),2)
    except ValueError:
        limit_max   = 'no'
        
    try:
        limit_av = round(float(out_text[1]),2)
    except ValueError:
        limit_av   = 'no'
        
    if limit_max != 'no'  and limit_av != 'no':

        limit = str(limit_av) + 'GeV·cm^(-2) averaged over declinations, and it is \n' + str(limit_max) + ' GeV·cm^(-2) at the most probable declination.'
    
    elif limit_max == 'no'  and limit_av != 'no':
        
        limit = str(limit_av) + 'GeV·cm^(-2) averaged over declinations.'
    
    elif limit_max != 'no'  and limit_av == 'no':
    
        limit = str(limit_max) + ' GeV·cm^(-2) at the most probable declination.'
        
    else:
        
        limit = str('is not provided for this GRB as its location is mostly outside the Auger FoV.')
    
    print("Maximun probability declination: ",declination, " the flux limit: ", limit)
    
    # Now we open the template as a str, including all the information that 
    # of the event
    
    with open('MT_Circular.txt', 'r', encoding='utf-8') as template_file:
            template_file_content = template_file.read()
            
            template_file = Template(template_file_content)
            
            MT_Circular = template_file.substitute(GWtime = str(GWtime.datetime), GWname = GWname, fov_prob_first = 100*round(fov_prob_first,2), name_trigger = GWname, 
                                                   limit = limit, GCN_ID = GCN_ID, t_max = t_max, maximum = maximum, time = time.strftime("%a, %d %b %Y %I:%M:%S %p", time.gmtime()) )
            
    
    # We finally send the email
    print("email")
        
    SDM.SendTheMail(SDM.sender, SDM.sendermail, SDM.receivers, SDM.receivermails, MT_Circular, attachments)

else:
    
    if "Limit_Flux" in Scripts:
        
        print('Executing Limit_Flux.py for', GWname)
        Flux_limit = subprocess.run([ 'python3', 'Limit_Flux.py', 'no', FILENAME], capture_output=True, text=True)
        
    if len(Scripts.split()) == 4: 
            
            # Now if we have launched everything in non GCN mode, a Circular-email
            # like is sent, again capturing the output of the Limit flux
            
            out_text = Flux_limit.stdout.split()
            print(Flux_limit.stderr)
            try:
                limit_max   = round(float(out_text[0]),2)
            except ValueError:
                limit_max   = 'no'
                
            try:
                limit_av = round(float(out_text[1]),2)
            except ValueError:
                limit_av   = 'no'
                
            if limit_max != 'yes'  and limit_av != 'yes ':

                limit = str(limit_av) + ' GeV·cm^(-2) averaged over declinations, and it is \n' + str(limit_max) + ' GeV·cm^(-2) at the most probable declination.'
            
            elif limit_max == 'yes'  and limit_av != 'yes ':
                
                limit = str(limit_av) + ' GeV·cm^(-2) averaged over declinations.'
            
            elif limit_max != 'yes'  and limit_av == 'yes ':
            
                limit = str(limit_max) + ' GeV·cm^(-2) at the most probable declination.'
                
            else:
                
                limit = str('is not provided for this GRB as its location is mostly outside the Auger FoV.')
             
            
            
            
            #print("declinacion del maximo: ",declination, " la cota es:", limit)
            
            with open('MT_Circular.txt', 'r', encoding='utf-8') as template_file:
                    template_file_content = template_file.read()
                    
                    template_file = Template(template_file_content)
                    
                    MT_Circular = template_file.substitute(GWtime = str(GWtime.datetime), GWname = GWname, fov_prob_first = 100*round(fov_prob_first,2), name_trigger = GWname, 
                                                           limit = limit, GCN_ID = GWname, t_max = t_max, maximum = maximum, time = time.strftime("%a, %d %b %Y %I:%M:%S %p", time.gmtime()) )
                    
                
            print("email")
                
            SDM.SendTheMail(SDM.sender, SDM.sendermail, SDM.receivers, SDM.receivermails, MT_Circular, attachments)
        
    elif 'email' in Scripts.split():
        
            # Otherwise the email is just send with the plots.
            
            MT = 'Resultados'
            
            print("email")
                
            SDM.SendTheMail(SDM.sender, SDM.sendermail, SDM.receivers, SDM.receivermails, MT, attachments)
            
    else:
        # If we don't want the email, we have finished 
            sys.exit()
