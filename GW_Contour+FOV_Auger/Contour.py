# -*- coding: utf-8 -*-
"""
Created on Mon Oct 10 20:33:46 2022

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
from astropy.io import fits as pyfits

import subprocess
import matplotlib.pyplot as plt
import matplotlib as mpl
import _pickle as pcl 

import time

scriptStartTime = time.time()

#FILENAME = "GW170814_skymap.fits.gz"
#FILENAME  = "bayestar.fits.gz,0"
FILENAME  = "LALInference_2.fits.gz"


DECREASE_RESOLUTION = False
DECREASE_NSIDE      = 512 # corresponds to 0.0131139632064 deg^2 per pixel ~ (0.115 deg)^2


#PREFIX = 'Auger/alerts/'
PREFIX = "FITS_files/"
SUFFIX = '.fits.gz'

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
    
# =============================================================================
# Constant definitions
# =============================================================================

SKY = 4 * 180**2 / np.pi
AUGERLAT, AUGERLONG, AUGERHEIGHT = -35.20666735*u.deg, -69.315833*u.deg, 1400*u.m
CL  = 0.9
regions = [ 'DGL', 'DGH', 'ES' ]
colors =  [ 'g',   'b',   'r' ]

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

prob, header = hp.read_map( FILENAME, h=True, verbose=False )

# Map properties

# Total number of pixels
npix = len(prob)

# Area per pixel, el NSIDE es un parametro relacionado con la resolucion pero no tengo claro que es
area_per_pix, nside = SKY / npix, hp.npix2nside(npix)

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
  npix = hp.nside2npix(nside)
  area_per_pix = SKY / npix

print(sum(prob))


# =============================================================================
# Our observation and observatory properties (time, coordinates)
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

    
# Este bucle simplemente encuentra la posicion en la cabecera en la que está 
# MJD-OBS  y la asigna a time_index

if time_index > -1 : 
    GWtime = astropy.time.Time( header[time_index][1], format='mjd' )
    
    # We construct the GWname from the date:
    
    GWname = "GW" + str(GWtime.datetime)[2:4] + str(GWtime.datetime)[5:7] + str(GWtime.datetime)[8:10]
else:               
    print('time_index invalid:',time_index,'could not find the MJD from the header file')
    sys.exit()

print(GWtime.datetime)


  
# Initial prob plot

plt.figure(1)
hp.mollview(
    prob,
    coord=["C"],
    title="Histogram equalized Galactic",
    unit="Probability",
    norm=None,
)
hp.graticule()
plt.savefig("Plots_Prob/Prob_plot"+GWname+".png")

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


# Obtain all thetas, phis of the pixels in all the map
theta, phi = hp.pix2ang(nside, np.arange(npix))

# Convert thetas, phis of all the pixels to RA, Dec.
radecs = astropy.coordinates.SkyCoord( ra=phi*u.rad, dec= (0.5*np.pi - theta) * u.rad )


# =============================================================================
# Obtaining the Contour
# =============================================================================

# First we obtain the points of the CL region

prob_above_min = np.array([0 if probi <= ConfReg_minimum_prob else probi for probi in prob])

print(np.sum(prob_above_min > 0 ))
print(npix_ConfReg)

# Region of the sky with nonzero is the CL region in the following plot

'''

plt.figure(2)
hp.mollview(
    prob_above_min,
    coord=["C"],
    title="Histogram equalized Galactic",
    unit="Probability",
    norm=None,
)
hp.graticule()
plt.savefig("CL_region.png")
'''

'''
Now that we have all the points in the 90% CL region, we can construct the
contour points as the points that change from zero to nonzero values
'''

contour = np.zeros(npix)
deroiz  = np.zeros(npix)

for i in range(npix-1):
    if (prob_above_min[i] == 0) and (prob_above_min[i+1] != 0):
        contour[i+1] = 1
        deroiz[i]
    elif (prob_above_min[i] != 0) and (prob_above_min[i+1] == 0):
        contour[i] = 1
    else:
        contour[i] == 0
        
print(np.sum(contour > 0 ))
        
''' Plot with hp.mollview that does not represent the contour well
plt.figure(3)
hp.mollview(
    contour,
    coord=["C"],
    title="Histogram equalized Galactic",
    unit="Probability",
    norm=None,
)
hp.graticule()
plt.savefig("3.png")
'''

dec_contour = np.array([ 0.5*np.pi - thetai for thetai, contouri in zip(theta,contour) if contouri != 0])
ar_contour   = np.array([ phii for phii, contouri in zip(phi,contour) if contouri != 0]) 


# =============================================================================
# Preparing the contour for the plot
# =============================================================================

def Ordering(theta,phi):
    
    '''
    theta : declination array of the contour
    
    phi : AR of the contour
    
    This function takes the declination (theta) and the AR(phi) positions of 
    the contours, that are not placed in an orderly manner, and them orders it 
    around the entire contour  by taking as the first element the one with the 
    greatest declination
    '''
    
    def dis(p1,p2):
        return np.sqrt((p1[0]-p2[0])**2+(p1[1]-p2[1])**2)
    
    n = len(theta)
    
    ord_th = np.zeros(n) 
    ord_ph = np.zeros(n) 
    
    ord_th[0] = theta[0]
    ord_ph[0] = phi[0]
    
    theta_0 = theta[0]
    phi_0   = phi[0]
    
    theta = np.delete(theta,0)
    phi   = np.delete(phi,0)
    
    count = 1
    
    while len(theta) > 1:
        
        distances = np.array([dis([phi_0,theta_0],[phii,thetai]) for phii,thetai in zip(phi,theta)])
        
        minimun = np.min(distances)
        pos_min = np.where(distances == minimun)[0][0]
            
        # cambiamos las posiciones de theta0 y phi0 por el minimo
        
        aux1 = theta[0]
        aux2 = phi[0]
        
        theta[0] = theta[pos_min]
        phi[0]   = phi[pos_min]
        
        theta[pos_min] = aux1
        phi[pos_min]   = aux2
    
        ord_th[count] = theta[0]
        ord_ph[count] = phi[0]
    
        theta_0 = theta[0]
        phi_0   = phi[0]
        
        theta = np.delete(theta,0)
        phi   = np.delete(phi,0)  
        
        count += 1
    
    ord_th[-1] = theta[0]
    ord_ph[-1] = phi[0]  
    
    
    return ord_th, ord_ph
    
'''
After ordering there still exists the posibility for the contour to be separated
separated into lobes, so we need to build a function that separate this lobes
'''

def Separate_Lobes(theta, phi, n_average):
    
    '''
    theta : array with the dec points of the Contours
    phi : array with the ar points of the Contours
    n_average : number of times the average of the separation between consecutive
    points of a Contour are needed to finish a contour and start the following new
    
    
    This function takes the ordered points of the Contour and separates it in
    different Lobes in the case that there exists
    
    returns lobes_theta and lobes_Phi which are lists of lists that contain the
    dec and ar points of every lobe.
    '''
    
    def dis(p1,p2):
        return np.sqrt((p1[0]-p2[0])**2+(p1[1]-p2[1])**2)
    
    lobes_theta = []
    lobes_phi   = []
    
    list_theta  = [theta[0],theta[1]] 
    list_phi    = [phi[0],phi[1]] 
    
    theta = np.delete(theta,0) 
    phi   = np.delete(phi,0)   
    
    suma = dis([list_theta[0],list_phi[0]], [list_theta[1],list_phi[1]])
    
    n = 1
    
    lobe = 0
    
    while len(theta) > 1 :
        
        average = suma/n
        
        new_distance = dis([theta[0],phi[0]], [theta[1],phi[1]])
        
        '''
        We calculate the distance of the new point and compare it with
        the average of distances multiplied by n_average, only if 
        new_distance >> average, a new lobe is created
        '''
        
        
        if new_distance < n_average*average :
                                    
            list_theta.append(theta[1])
            list_phi.append(phi[1])
            
            theta = np.delete(theta,0) 
            phi   = np.delete(phi,0)
            
            n = n + 1
            suma = suma + new_distance
            
            '''
            In the case that len(theta) == 1 after the previus delete. We calculate the distance
            of the last point and the point in theta/phi and decide if it is necessary 
            to create a new lobe with one point (that after all will not be taken into account)
            or include it in the last lobe
            '''
            
            if len(theta) == 1 :
                 new_distance = dis([list_theta[-1],list_phi[-1]], [theta[0],phi[0]])
                 
                 if new_distance < n_average*average : 
                     
                     
                     lobe = lobe + 1 ; print("Lobe", lobe)
                     
                     list_theta.append(theta[0])
                     list_phi.append(phi[0])
                     
                     lobes_theta.append(list_theta)
                     lobes_phi.append(list_phi)
                 else : 
                     
                     lobe = lobe + 1 ; print("Lobe", lobe)
                     
                     lobes_theta.append([theta[0]])
                     lobes_phi.append([phi[0]])
            
        else:
            
            lobe = lobe + 1 ; print("Lobe", lobe)
            
            theta = np.delete(theta,0) 
            phi   = np.delete(phi,0)
            
            lobes_theta.append(list_theta)
            lobes_phi.append(list_phi)
            
            '''
            It exists the posibility to entry here with len(theta) = 2, then
            after the two deletes from above we arrive to a situation with 
            len(theta) = 1 which is treated separately
            '''
            
            if len(theta) > 1: 
                
                new_distance = dis([theta[0],phi[0]], [theta[1],phi[1]])
                
                '''
                This while tries to take into account posible points not contained
                in lobes that could exist between lobes after the ordering procedure
                '''
                
                while new_distance >= n_average * average :
                    
                    lobe = lobe + 1 ; print("Lobe", lobe)
                    
                   
                    lobes_theta.append([theta[0]])
                    lobes_phi.append([phi[0]])
                    
                    theta = np.delete(theta,0) 
                    phi   = np.delete(phi,0) 
                    
                    new_distance = dis([theta[0],phi[0]], [theta[1],phi[1]])
                    
                    
                
                if new_distance < n_average*average : 
            
                    list_theta = [theta[0], theta[1]]
                    list_phi   = [phi[0], phi[1]]
                    
                    theta = np.delete(theta,0) 
                    phi   = np.delete(phi,0)   
                    
                    suma = dis([list_theta[0],list_phi[0]], [list_theta[1],list_phi[1]])
    
                    n = 1
                    
            elif len(theta) == 1:
                
                lobe = lobe + 1 ; print("Lobe", lobe)
                
                
                lobes_theta.append([theta[0]])
                lobes_phi.append([phi[0]])
        
        
    " We eliminate possible 1 length lobes (I finally consider this unnecessary) "
    #lobes_theta = [lobesi for lobesi in lobes_theta if len(lobesi) > 1]
    #lobes_phi = [lobesi for lobesi in lobes_phi if len(lobesi) > 1]
        
                
    return(lobes_theta, lobes_phi)

dec_contour_new, ar_contour_new = Ordering(dec_contour, ar_contour)
lobes_dec, lobes_ar = Separate_Lobes(dec_contour_new, ar_contour_new, 8)

# =============================================================================
# Plane plot
# =============================================================================

mpl.rc('lines', linewidth=2.0)
'''

#plt.figure(5)


#plt.plot(ar_contour - np.pi, dec_contour,"b." )

#plt.savefig("Contour_probe.png")

# In order to prove that the Contour is well traversed

plt.figure(5)
for i in range(len(ar_contour_new)):
    
    if i%100 < 50:
    
        plt.plot(ar_contour_new[i], dec_contour_new[i],"r." , markersize = 0.1, alpha = 1)    
    else:
        plt.plot(ar_contour_new[i], dec_contour_new[i],"b." , markersize = 0.1, alpha = 1)
    if i%50 == 0 : 
        plt.savefig("Contour_probe%s.png"%i)

plt.plot(ar_contour_new - np.pi, dec_contour_new,"-b",  label = "CL contour")
plt.plot((ar_contour_new - np.pi)[-2], dec_contour_new[-2],"*k")

plt.legend()

plt.savefig("Contour_probe.png")

'''
# =============================================================================
# Plot of the contour mollweide
# =============================================================================

plt.clf()

fig = plt.figure(4,figsize=(10,6))

ax = fig.add_subplot((111), projection="mollweide")
ax.set_xticklabels(['2h','4h','6h','8h','10h','12h','14h','16h','18h','20h','22h'])
ax.grid(True)

# We substract pi in order to have the 0 h at the left of the image
plt.plot(ar_contour_new - np.pi, dec_contour_new,"-b",  label = "CL contour")
n = 0
plt.savefig("Plots_Contour/Contour"+GWname+".png")

plt.clf()

ax = fig.add_subplot((111), projection="mollweide")
ax.set_xticklabels(['2h','4h','6h','8h','10h','12h','14h','16h','18h','20h','22h'])
ax.grid(True)


cmap = plt.get_cmap('gnuplot')
c = [cmap(i) for i in np.linspace(0, 1, len(lobes_ar))]
for i,j,ci in zip(lobes_ar, lobes_dec, c): 

    plt.plot(np.array(i) - np.pi, np.array(j), "-", label = "CL contour %s" %n, color = ci)
    n = n + 1
    


plt.legend(bbox_to_anchor=(1.12, 1.18), loc='upper right', borderaxespad=1)

#plt.scatter(ar_contour - np.pi, dec_contour, s = 0.005)

plt.savefig("Plots_Contour/ContourLobes"+GWname+".png")

# =============================================================================
# Now we save the data
# =============================================================================

# Delete old data

file_dec = open("Contour_Data/decContour"+GWname+".txt", "w")
file_ar  = open("Contour_Data/arContour"+GWname+".txt", "w")
file_dec.close()
file_ar.close()
    
# Introduce new data

file_dec = open("Contour_Data/decContour"+GWname+".txt", "a")
file_ar  = open("Contour_Data/arContour"+GWname+".txt", "a")

def write(a, file):
    
    for i in a:
        
        file.write(str(i))
        file.write(" ")
        

for lobe_dec, lobe_ar in zip(lobes_dec,lobes_ar):
    
    write(lobe_dec, file_dec)
    file_dec.write("\n")
    
    write(lobe_ar, file_ar)
    file_ar.write("\n")
    
file_dec.close()
file_ar.close()
    
#np.savetxt("Pruebas/arContourGW"+GWname+".txt", dec_contour_new)
#np.savetxt("Pruebas/decContourGW"+GWname+".txt", ar_contour_new )


