# -*- coding: utf-8 -*-
"""
Created on Sat Oct  8 11:43:14 2022

@author: yagol
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

import matplotlib.pyplot as plt
import matplotlib as mpl
import _pickle as pcl 

import time

# =============================================================================
# Manual parameters to the program
# =============================================================================

scriptStartTime = time.time()

PLOTvsTIME = True

#FILENAME = "GW170814_skymap.fits.gz"
FILENAME  = "bayestar.fits.gz,0"

DECREASE_RESOLUTION = False
DECREASE_NSIDE      = 512 # corresponds to 0.0131139632064 deg^2 per pixel ~ (0.115 deg)^2

#SCAN_RESOLUTION = 1./2.  # deg     # Scan every 2 minutes
#SCAN_RESOLUTION  = 1.     # deg     # Scan every 4 minutes
SCAN_RESOLUTION = 15.    # deg     # Scan every hour
#SCAN_RESOLUTION = 360.  # deg     # 1 point only at the time of the event



#PREFIX = 'Auger/alerts/'
SUFFIX = '.fits.gz'
PREFIX = "FITS_files/"
PICKLE_PREFIX = "Pickle/"


PICKLE_FILENAME = PICKLE_PREFIX + FILENAME


if PREFIX not in FILENAME:
  FILENAME = PREFIX + FILENAME
  
if SUFFIX not in FILENAME:
  SUFFIX='.fits' # assuming you used the unpacked fits file
  
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
    '''
    outputText : str to print and save
    outfile : str with the name of the file where outputText is saved
    
    This function prints outputText and then saves it in outputfile
    '''
    print(outputText)
    outfile.write(outputText+'\n')
    
# =============================================================================
# Constant definitions
# =============================================================================

SKY = 4 * 180**2 / np.pi
AUGERLAT, AUGERLONG, AUGERHEIGHT = -35.20666735*u.deg, -69.315833*u.deg, 1400*u.m
CL  = 0.9
regions = [ 'DGL', 'DGH', 'ES' ]
colors  =  [ 'g',   'b',   'r' ]

# =============================================================================
# Start here
# ============================================================================

# Name of output file (text format) coge el nombre original, le elimina el sufijo 
# y añade un nuevo sufijo _out.txt y lo abre para introducir informacion


outfile = open(FILENAME[:-len(SUFFIX)]+'_out.txt','w')


# =============================================================================
# Read Sky map
# =============================================================================

# Download sky map. Not implemented now
#~ subprocess.check_call([    'curl', '-O', '--netrc',    'https://gracedb.ligo.org/apiweb/events/G298389/files/LIB.fits.gz,0'])


'''
Como su nombre indica hp.read_map es una funcion que nos permite leer un 
healpix map de un archivo fits. Partial-sky files, if properly identified, 
are expanded to full size and filled with UNSEEN.

FILENAME el archivo fits

h = True se habilita para que devuelva también la cabecera

verbose deprecated ?? no hace ningun efecto segun la documentacion

'''

# prob is the unidimensional array containing the probability in each pixel

print(FILENAME)

prob, header = hp.read_map( FILENAME, h=True, verbose=False )

# Map properties

# Total number of pixels
npix = len(prob)
print_and_write('Total number of pixels in map: ' + str( npix ), outfile)

# Area per pixel, el NSIDE es un parametro relacionado con la resolucion pero no tengo claro que es
area_per_pix, nside = SKY / npix, hp.npix2nside(npix)
print_and_write('Area per pixel (before eventual resolution decrease): '+str( area_per_pix )+' deg^2',outfile)


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


print("The sum of the probabilities is " + str(sum(prob)))

# =============================================================================
#  Confidence Region of sky localization corresponding to CL
# =============================================================================

# First obtain the minimum value of prob array above the chosen CL
ConfReg_minimum_prob = get_ConfReg_minimum_prob(prob,CL)
# Count how many pixels are above the given minimum prob corresponding to the chosen CL

'''
Contar el numero de elementos que satisfacen una condicion se puede hacer con
np.sum tal como se hace aquí. si tenemos a = np.array([0,1,2,3,4])

np.sum( a < 2 ) dara como resultado 2, pues solo hay dos elementos en el
array tales que sean menores que 2.

Luego yo veo que npix_ConfReg es el numero de pixeles cuya probabilidad es
mayor que ConfReg_minimum_prob

'''
npix_ConfReg = np.sum( prob > ConfReg_minimum_prob )


print_and_write('Number of pixels in the confidence region: ' + str( npix_ConfReg ), outfile)
print_and_write('Solid angle of the confidence region:      ' + str( area_per_pix*npix_ConfReg )+' deg^2', outfile)
 
# Obtain all thetas, phis of the pixels in all the map
theta, phi = hp.pix2ang(nside, np.arange(npix))

# Convert thetas, phis of all the pixels to RA, Dec.
radecs = astropy.coordinates.SkyCoord( ra=phi*u.rad, dec= (0.5*np.pi - theta) * u.rad )


# =============================================================================
# Our observation and observatory properties (time, coordinates, event)
# =============================================================================

# bucle sobre el cabecero del fichero FITS

# Find the time (MJD = Modified Julian Date: A continuous measure in days since 
#midnight at the start of 17 November 1858. Based on UTC) in the header

timeindex = -1
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
    
    # We construct the GWname from the date:
    
    GWname = "GW" + str(GWtime.datetime)[2:4] + str(GWtime.datetime)[5:7] + str(GWtime.datetime)[8:10]
else:               
    print('time_index invalid:',time_index,'could not find the MJD from the header file')
    sys.exit()
    


'''   In order to extract the name from the header (some fits doesn´t include it)
if object_index > -1 : 
    GWname = header[object_index][1]
else:               
    print('object_index invalid:',object_index,'could not find the Object type from the header file')
    sys.exit()
'''

observatory = astropy.coordinates.EarthLocation( lat=AUGERLAT, lon=AUGERLONG, height=AUGERHEIGHT )

# =============================================================================
# Loop through the MJD during a day
# =============================================================================

# The number of iterations of the first loop will be taken into account
mjdindex = 0
fov_prob = []
fov_frac = []


'''
Since the angular resolution is SCAN_RESOLUTION in degrees, the number of scans
that we make is 360/SCAN_RESOLUTION + 1 (taking into account the 0 and the 360)
'''
NSCANPOINTS = 360/SCAN_RESOLUTION+1 - (SCAN_RESOLUTION == 360) 
NSCANPOINTS = int(NSCANPOINTS)
# This last expression in the parentheses makes argument of the following for-loop get only one element 
# if the SCAN_RESOLUTION is 360, i.e. only the actual time of the GW event is used

ifirst = 1

# ----- Start loop in MJD for 1 day

# The linspace makes equidistant time points between GWtime.mjd and the 
# following day GWtime.mjd+1

for mjd in np.linspace( GWtime.mjd, GWtime.mjd+1, NSCANPOINTS ):
    currtime = astropy.time.Time( mjd, format='mjd' )
    print_and_write( '\n\n'+str(currtime.datetime)+'\n' , outfile)
    
    #~ Observatory frame (unfortunately has to be redone for each MJD in the loop)
    
    frame = astropy.coordinates.AltAz( obstime=currtime, location=observatory )
    
    # radecs was the phis and thetas in RA and Dec from de fits, we change
    # the coordinates system to the AltAz from the observatory with
    # (This is probably the longest taking step :/)
    
    altaz = radecs.transform_to(frame)
    
    #~ FoV definition - set of conditions fulfiled by each FoV
    fov_DGL = ( altaz.alt <= 30*u.deg ) & ( altaz.alt > 15*u.deg )
    fov_DGH = ( altaz.alt <= 15*u.deg ) & ( altaz.alt >  0*u.deg )
    fov_ES  = ( altaz.alt <=  0*u.deg ) & ( altaz.alt > -5*u.deg )
    
    all_fovs = [fov_DGL,fov_DGH,fov_ES]
    
    sum_all_fovs = [np.sum(fovi) for fovi in all_fovs]
  
    # Do the following only once
    if ifirst == 1:

       print('Fraction of sky in fov ES  [%]: ',sum_all_fovs[2]*100./len(fov_ES))
       print('Fraction of sky in fov DGH [%]: ',sum_all_fovs[1]*100./len(fov_DGH))
       print('Fraction of sky in fov DGL [%]: ',sum_all_fovs[0]*100./len(fov_DGL))
       print(' ')
       ifirst = -1 
       
    
    '''
    Example to understand the following line
    
    a = np.array([1,2,3,4,5])
    b = np.array([0,4])
    
    then a[b] will be np.array([1,5])
    
    '''
    
    #~ Create 3 arrays containing the values of prob in the pixels that are in each FoV (ES, DGH, DGL) 
    all_prob_fovs = [prob[sub_fov] for sub_fov in all_fovs]
    
    #~ Loop through the FoV and get:
    # ---------------
    # Total probability (sum of the probabilities assigned to each pixel)
    # Total probability (sum of area of pixels divided by total area, i.e. 
    # all pixels have equal probability)
    # ---------------
    
    #~ Total probability of ##CL contour## inside each FoV
    # We append empty lists
    fov_prob.append([])
    #~ Fraction of CL contour inside each FoV
    fov_frac.append([])
    
    # for loop trough the three FOV´s
    
    for i in range(len(all_prob_fovs)):
        
        # First we access the empty list and add the sum of probabilities
        # in the fov of each 
        fov_prob[mjdindex].append(np.sum(all_prob_fovs[i]))
        print_and_write('Total probability (sum of prob. pixels) for the event to be in the '+str(regions[i])+' FoV: '+str( fov_prob[mjdindex][i] ), outfile)
        
        
        fov_frac[mjdindex].append( float(np.sum(all_prob_fovs[i] > ConfReg_minimum_prob)) / npix_ConfReg )
        print_and_write('Fraction of 90 % confidence region in the '+str(regions[i])+' FoV:    '+str( fov_frac[mjdindex][i] ), outfile)
        print_and_write('', outfile)
        
    #~ Summary
    # We add the three probabilites to be in the DGH DGL and ES regions
    print_and_write('Total probability for the event to be in the total FoV: '+str( np.sum(fov_prob[mjdindex]) ), outfile)
    print_and_write('Fraction of 90 % confidence region in the total FoV:    '+str( np.sum(fov_frac[mjdindex]) ), outfile)
    
    mjdindex += 1
    
fov_prob_first = np.sum(fov_prob[0])
fov_frac_first = np.sum(fov_frac[0])

print('------------------------------------------------------------------------------------------------------')
print_and_write('\n\nMaximum probability (sum of prob. pixels) for the event to be in the total FoV:    '+str( max(np.sum(fov_prob,axis=1)) )+' '+str( (np.sum(fov_prob,axis=1).argmax()/360.*SCAN_RESOLUTION*24)%24 )+' hours after the GW event.', outfile)
print('------------------------------------------------------------------------------------------------------')
print('Initial fraction of 90 % confidence region in the total FoV: ',fov_frac_first)

print_and_write('Maximum fraction of 90 % confidence region in the total FoV: '+str( max(np.sum(fov_frac,axis=1)) )+' '+str( (np.sum(fov_frac,axis=1).argmax()/360.*SCAN_RESOLUTION*24)%24 )+' hours after the GW event.\n', outfile)

outfile.close()


# =============================================================================
# Output files in pickle format
# =============================================================================

for i in range(len(regions)):
  pcl.dump( np.array([np.linspace( 0, 24, NSCANPOINTS ), np.array(fov_prob)[:,i]]), open(PICKLE_FILENAME[:-len(SUFFIX)]+'_prob_vs_time_'+regions[i]+GWname+'.pickle','wb') )
  pcl.dump( np.array([np.linspace( 0, 24, NSCANPOINTS ), np.array(fov_frac)[:,i]]), open(PICKLE_FILENAME[:-len(SUFFIX)]+'_frac_vs_time_'+regions[i]+GWname+'.pickle'+GWname,'wb') )


# =============================================================================
# Plots of:
# Sum of probabilities of localization map inside each field-of-view (ES, DGH, DGL) for each MJD time
# Fraction of localization map inside each field-of-view (ES, DGH, DGL) for each MJD time
# =============================================================================

mpl.rc('lines', linewidth=2.0)

if PLOTvsTIME:
    
  xs = np.linspace( 0, 24, NSCANPOINTS ) # Per hour
  fs1, fs2 = 16, 12  # Font sizes
  
  # Array of two subplots
  #You can specify nrows and ncols. sharex allows to share the x axis
  fig, axarr = plt.subplots( 2, sharex=True ,figsize=(8,5) )

# ---------------------------- Top subplot -----------------------------------

# Sum of probabilities of localization map inside each field-of-view (ES, DGH, DGL) for each MJD time
# and for the sum of the 3 FoV

  axarr[0].plot( xs  , np.sum(fov_prob,axis=1), color='black', label="sum" )
  for i in range(len(regions)):
    axarr[0].plot( xs, np.array(fov_prob)[:,i], label=regions[i], color=colors[i] )

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

print('----------------------FINISH--------------------------------')
print('approx. run time of the script', time.time()-scriptStartTime, 's')

plt.savefig("CL_Coverage/GW_confidence_region_coverage" + GWname + ".png")
# plt.savefig("CL_Coverage/GW_confidence_region_coverage" + nameGW + ".pdf")


