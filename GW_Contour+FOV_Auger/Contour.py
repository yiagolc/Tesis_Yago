# -*- coding: utf-8 -*-
"""
This program is intended to obtain the Contours of 90% confidence level from
a GW event. To do that it only needs the FITS_file of the event. The script 
reads it, calculates the contour points and then sorts them following the contour.
Finally separates in different "lobes", because the contour can be splitted in 
different regions.

It may happen that, for optimal results, the parameter_lobes  parameter should be 
different for each event, but even with suboptimal performance, the contour plot
with all lobes in black colour looks approximately correct. In the event of
bad performance, trick with parameter_lobes in the Runner.py.

The Separate_Lobes function compares the average distance between points in a lobe
and the distance of a new point with respect to the last one, when following the
Ordered points array. So parameter_lobes is just the number of times greater than 
the distance to a new point should be with respect to the average in order to split
and create a new lobe.

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
import os
import subprocess

from numba import njit

#import time

#print(os.getpid())

currentdir = os.getcwd()

#scriptStartTime = time.time()

# This is just in the case we want to decrease the resolution of the healpy map

DECREASE_RESOLUTION = False
DECREASE_NSIDE      = 512 # corresponds to 0.0131139632064 deg^2 per pixel ~ (0.115 deg)^2

#GCN_Mode = 'yes'


# =============================================================================
# Functions
# =============================================================================


def get_ConfReg_minimum_prob(prob,CL):
    '''    
    Parameters
    ----------
    prob : 1D np.ndarray 
        This array contains the probability as a function of declination and
        right ascencion, but using healpix pixelization schemes.
    CL : float
        Confidence level that we want to consider, tipically for the GCN automation
        we take CL = 0.9
    Returns
    -------
    i : int
        This function takes prob and sums all the elements from the largest to 
        the lowest, until the sum is greater than CL. Then returns the prob element
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
    
    
# functions for constructing the contours

@njit
def Ordering(theta,phi):
    '''
    

    Parameters
    ----------
    theta : 1D np.ndarray
        Declination coordinates of all points along the contour.
    phi : 1D np.ndarray
        Right Ascencion coordinates of all points along the contour

    Returns
    -------
    ord_th, ord_ph 1D np.ndarray
        This function takes the declination (theta) and the AR(phi) positions of 
        the contours, that are not placed in an orderly manner, and them orders it 
        along the entire contour clockwise or anticlockwise by taking as the 
        first element the one with the greatest declination. If the contour
        it's splitted in several lobes, it just jumps from one to the following
        one.

    '''   
    
    def dis(p1,p2):
        # Angular distance between two points
        return np.sqrt((p1[0]-p2[0])**2+(p1[1]-p2[1])**2)
    
    def zipp(a1,a2):
        # We define the zip here in order to use it along with the njit decorator
        n = len(a1)
        a = []
        
        for i in range(n):
            
            a.append([a1[i],a2[i]])
            
        return a
        
        
    
    n = len(theta)
    
    # We start creating the ordered arrays, adding the first points, and creating
    # auxiliary variables, that we call theta_0, phi_0. Recall that in healpy 
    # the first point is the northernmost point
    
    ord_th = np.zeros(n) 
    ord_ph = np.zeros(n) 
    
    ord_th[0] = theta[0]
    ord_ph[0] = phi[0]
    
    theta_0 = theta[0]
    phi_0   = phi[0]
    
    # We will be deleting the points from the original arrays once thay are passed
    # two the ordered arrays
    
    theta = np.delete(theta,0)
    phi   = np.delete(phi,0)
    
    count = 1
    
    while len(theta) > 1:
        
        # We compute the distance of every point to the last one in the ordered
        # list. We take the closer one to be the following in the ordering algorythm
        
        distances = np.array([dis([phi_0,theta_0],[phii,thetai]) for phii,thetai in zipp(phi,theta)])
        
        minimun = np.min(distances)
        pos_min = np.where(distances == minimun)[0][0]
            
        # update the auxiliary variables
        
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
        
        # delete the new point from the old list
        
        theta = np.delete(theta,0)
        phi   = np.delete(phi,0)  
        
        count += 1
        
    # We have to include also the last point in the ordered array
    
    ord_th[-1] = theta[0]
    ord_ph[-1] = phi[0]  
    
    
    return ord_th, ord_ph
    
'''
After ordering there still exists the posibility for the contour to be separated
separated into lobes, so we need to build a function that separate this lobes
'''

@njit
def Separate_Lobes(theta, phi, n_average):
    '''
    Parameters
    ----------
    theta : 1D np.ndarray
        Declination coordinates of ordered points along the contour.
    phi : 1D np.ndarray
        Right ascencion coordinates of ordered points along the contour.
    n_average : float
        number of times that the distance to a new point respect to the latest
        has to be with respect to the average distance of the last lobe
        in order to create a new lobe, tipically 8

    Returns
    -------
    lobes_theta, lobes_phi  : list of lists
    
        This function takes the ordered points of the Contour and separates it in
        different Lobes in the case that there exists. returns lobes_theta and 
        lobes_Phi which are lists of lists that contain the dec and ar points
        of every lobe.
    '''
   

    def dis(p1,p2):
        # Angular distance between two points
        return np.sqrt((p1[0]-p2[0])**2+(p1[1]-p2[1])**2)
    
    def zipp(a1,a2):
        # We define the zip here in order to use it along with the njit decorator
        n = len(a1)
        a = []
        
        for i in range(n):
            
            a.append([a1[i],a2[i]])
            
        return a
    
    # We start by creating the lobes lists
    
    lobes_theta = []
    lobes_phi   = []
    
    # list_theta and list_phi will be auxiliary lists, corresponding to the points
    # of a new lobe. We include the first point, and we delet it from the original
    # array
    
    list_theta  = [theta[0],theta[1]] 
    list_phi    = [phi[0],phi[1]] 
    
    theta = np.delete(theta,0) 
    phi   = np.delete(phi,0)   
    
    # we are going to create a variable called suma, which then we will divide
    # by the number of points considered minus one, which will be the number 
    # of distances averaged
    
    suma = dis([list_theta[0],list_phi[0]], [list_theta[1],list_phi[1]])
    
    # Number of points in the lobe and number of lobes
    
    n = 1
    lobe = 0
    
    while len(theta) > 1 :
        
        average = suma/n
        
        new_distance = dis([theta[0],phi[0]], [theta[1],phi[1]])
        
        # We calculate the distance of the new point and compare it with
        # the average of distances multiplied by n_average, only if 
        # new_distance >> average, a new lobe is created
        
        if new_distance < n_average*average :
            
            # Here a new lobe is not created, so we add the point to the lobe
            # in which we are, and delete from the original list
                                    
            list_theta.append(theta[1])
            list_phi.append(phi[1])
            
            theta = np.delete(theta,0) 
            phi   = np.delete(phi,0)
            
            n = n + 1
            suma = suma + new_distance
            
            # If after this lines we have just one point, the while loop stops. 
            # Then lets introduce here that case. It may happen that the ordering
            # algorithm does not work perfectly. For instace if we have something like
            # the following:
            #              .
            #           .  .
            #        .     . 
            #     .        .
            # If we are ordering this along the contour using the nearest point
            # and clokwise direction, then the corner point would be skipped.
            # In that case this point might end up in the end of ordered arrays.
            # Taking this into account we check whether this point is far from the 
            # last point of the contour or not, just like we have been doing

            
            if len(theta) == 1 :
                 new_distance = dis([list_theta[-1],list_phi[-1]], [theta[0],phi[0]])
                 
                 if new_distance < n_average*average : 
                     
                     # If the last point is near, then we include it in the lobe
                     # and we have finished
                                         
                     lobe = lobe + 1 ; #print("Lobe", lobe)
                     
                     list_theta.append(theta[0])
                     list_phi.append(phi[0])
                     
                     lobes_theta.append(list_theta)
                     lobes_phi.append(list_phi)
                     
                 else : 
                     
                     # Otherwise we create a new lobe with just this point.
                     # This of course would be bad if we were trying to do 
                     # a perfect splitting. But in order to plot the contours
                     # that's enough, we won't make any distinction between
                     # lobes, all will be plotted in black
                     
                     lobe = lobe + 1 ; #print("Lobe", lobe)
                     
                     lobes_theta.append([theta[0]])
                     lobes_phi.append([phi[0]])
            
        else:
            
            # In the case that the new distance is greater, then we first consider
            # the last point of the finished lobe to be the last included, which is still
            # in the original arrays in order to compute the distance to the following.
            # Then it should be eliminated from these lists. Besides, we include 
            # this finished lobes in the lobes lists
            
            lobe = lobe + 1 ; #print("Lobe", lobe)
            
            theta = np.delete(theta,0) 
            phi   = np.delete(phi,0)
            
            lobes_theta.append(list_theta)
            lobes_phi.append(list_phi)
            
            # Now again, like we did previously we take care of the case of len(theta) = 1
            # after this lines, in the event that it happens

            
            if len(theta) > 1: 
                
                new_distance = dis([theta[0],phi[0]], [theta[1],phi[1]])
                
                # The following while loop stands, for the possibility that
                # the situation of bad ordering performance previously described
                # happens more than just once. In that case it would be neccessary
                # to creat new lobes, each for each independent  bad-ordered point.
                # Even if the bad reconstructed point are between lobes and not 
                # just in the end this loop will take into account that
            
                
                while new_distance >= n_average * average :
                    
                    # The new independent points lobes, will be created while the
                    # distance to the new point from the last continues to be greater
                    # to the previous average in a lobe.
                    
                    lobe = lobe + 1 ; #print("Lobe", lobe)
                   
                    lobes_theta.append([theta[0]])
                    lobes_phi.append([phi[0]])
                    
                    theta = np.delete(theta,0) 
                    phi   = np.delete(phi,0) 
                    
                    # We could reach here the case of finishing the points, so
                    # we treat again the len(theta) = 1 case
                    
                    if len(theta) > 1:
                        new_distance = dis([theta[0],phi[0]], [theta[1],phi[1]])
                        
                    else:
                       lobe = lobe + 1 ; #print("Lobe", lobe) 
                       
                       lobes_theta.append([theta[0]])
                       lobes_phi.append([phi[0]])
                       
                       return(lobes_theta, lobes_phi)
                    
                
                if new_distance < n_average*average : 
                    
                    # In this case we can start new lobes. So we restart the list
                    # theta and the list_phi arrays with the new point.
            
                    list_theta = [theta[0], theta[1]]
                    list_phi   = [phi[0], phi[1]]
                    
                    theta = np.delete(theta,0) 
                    phi   = np.delete(phi,0)   
                    
                    suma = dis([list_theta[0],list_phi[0]], [list_theta[1],list_phi[1]])
    
                    n = 1
                    
            elif len(theta) == 1:
                
                lobe = lobe + 1 ; #print("Lobe", lobe)
                
                
                lobes_theta.append([theta[0]])
                lobes_phi.append([phi[0]])
                
                return(lobes_theta, lobes_phi)        
                
    return(lobes_theta, lobes_phi)


# =============================================================================
#%% Input information, from Listening_kafka.py or Runner.py
# =============================================================================


if(len(sys.argv) < 2) :
    print('You have not include a fits file in the input')
    sys.exit()

else:
    
    # Two situations here first, are we in GCN_mode?
    GCN_mode = sys.argv[1]
    
    if GCN_mode == 'yes':
        
        # If we are in GCN mode then the parameters are fixed accordingly to what 
        # it is stipulated. We read from the nput just the file name and the ID
        
        list_path_file = [sys.argv[2]]
        GCN_ID = sys.argv[3]
        
        Scripts = "FOV_GW CL_Coverage Limit_Flux email" # In GCN mode we want all
        parameters = '1. no yes' # all other parameters in GCN mode
        parameter_lobes = 8  # This has to be fixed in this script
        
    else:
        
        # If we are not in GCN mode, then all the parameters and the Scripts
        # that we have to run are included in the information of the input,
        # as it was established in the Runner.py Script

        Scripts         = sys.argv[2]
        parameters      = sys.argv[3]
        
        # We keep this two parameters, but for this Script we only need the 
        # second one
        SCAN_RESOLUTION = float(parameters.split()[1])
        parameter_lobes = float(parameters.split()[0])
        
        # We eliminate parameter_lobes from the parameter list, because it 
        # won't be needed anymore
        parameters_list = parameters.split()[1:]
        parameters = parameters_list[0]
        
        for i in range(len(parameters_list) - 1):
            i = i + 1
            
            parameters = parameters + ' ' + parameters_list[i]
            
        # When we are not in the GCN mode, we will run all the established 
        # scripts and parameters over all the FITS files in the /Data_FITS
        # folder. We build here a list with all the files and paths
        
        path = os.path.join(currentdir, "Data_FITS")
        
        text_in_file = ".fits"
        #text_in_file = "GW170814_skymap.fits.gz"
        #text_in_file = ".gz"
        
        list_path_file = []   
        list_file      = [] 
        
        for dirpath, dirnames, filenames in os.walk(path):
            for filename in [f for f in filenames if text_in_file in f]:
                
                path_with_filename_e=os.path.join(dirpath, filename)
                list_path_file.append(path_with_filename_e)
                list_file.append(filename)
            list_path_file.sort() 
            list_file.sort()  

#list_path_file = [os.path.join(script_dir, "GW170608_skymap.fits.gz")]

# In the case that we don't want the FOV+Contours plots, then we just don't run
# the rest of the script, we just start running the CL_coverage for every FITS file 


if "FOV_GW" not in Scripts.split():
    
    for FILENAME in list_path_file:
        
         process = subprocess.run([ 'python3', 'CL_coverage.py','no', FILENAME, Scripts, str(SCAN_RESOLUTION)])
        
    sys.exit()  # if FOV_GW is not in the Scripts str, we don´t want to run Contour or fov_Auger_builder
                # we stop here
      
# =============================================================================
# Constant definitions
# =============================================================================

SKY = 4 * 180**2 / np.pi
AUGERLAT, AUGERLONG, AUGERHEIGHT = -35.20666735*u.deg, -69.315833*u.deg, 1400*u.m
CL  = 0.9
regions = [ 'DGL', 'DGH', 'ES' ]
colors =  [ 'g',   'b',   'r' ]



for FILENAME in list_path_file:

    # =========================================================================
    # Read Sky map
    # =========================================================================
    
    # Download sky map. Not implemented now
    #~ subprocess.check_call([    'curl', '-O', '--netrc',    'https://gracedb.ligo.org/apiweb/events/G298389/files/LIB.fits.gz,0'])
    
    

    # As the name suggests hp.read_map is a function that allows as reading a healpix
    # map from a FITS file.Partial-sky files, if properly identified, 
    # are expanded to full size and filled with UNSEEN. With h = True we store
    # the header too
    
    
    # prob is the unidimensional array containing the probability in each pixel

    prob, header = hp.read_map( FILENAME, h=True)
    
    ######### Map properties #########
    
    # Total number of pixels
    npix = len(prob)
    
    # Area per pixel,  NSIDE is a parameter related to the resolution 
    area_per_pix, nside = SKY / npix, hp.npix2nside(npix)
    
    # The following lines are in the case that we consider a Decrease in resolution
    # which is not interesting for us in this Automation
    
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
    
    
    # =========================================================================
    # Our observation and observatory properties (time, coordinates)
    # =========================================================================
    
    # loop over the header of the FITS file
    
    # Find the time (MJD = Modified Julian Date: A continuous measure in days since 
    # midnight at the start of 17 November 1858. Based on UTC) in the header
    
    # First we find the index of the MJD in case it is provided
    
    timeindex = -1
    for i in range(len(header)):
      try:
        header[i].index('MJD-OBS') 
        time_index = i
        
      except:
          pass
    
    
    if time_index > -1 : 
        
        # Start using astropy to define dates
        
        GWtime = astropy.time.Time( header[time_index][1], format='mjd' )
        
        # GWname its GCN_ID if GCN_mode == yes, otherwise we construct the name from the time
    
        if GCN_mode == 'yes' : 
            
            GWname = GCN_ID
        else:
            
            GWname = "GW" + str(GWtime.datetime)[2:4] + str(GWtime.datetime)[5:7] + str(GWtime.datetime)[8:10] + "_" + str(GWtime.datetime)[11:13]+str(GWtime.datetime)[14:16]
    else:               
        print('time_index invalid:',time_index,'could not find the MJD from the header file')
        sys.exit()
        
    print('Calculation Contour for', GWname)
      
    # This code is just to plot the probability in equatorial coordinates
    
    plt.figure(1)
    
    hp.projview(
    prob,
    coord=["C"],
    rot= (-180,0,0),
    graticule=True,
    graticule_labels=True,
    phi_convention = 'clockwise',
    flip = 'geo',
    unit="probability",
    xlabel="RA",
    ylabel="declination",
    cb_orientation="vertical",
    latitude_grid_spacing=15,
    longitude_grid_spacing=45,
    projection_type="mollweide",
    title="Mollweide projection",
    #flip="geo",
    #phi_convention="clockwise",
    xtick_label_color='white'
    )
    #plt.subplot(projection="mollweide")
    plt.tight_layout()
    try:
        plt.savefig("Plots_Prob/Prob_plot"+GWname+".png")
    
    except:
        os.makedirs("Plots_Prob")
        plt.savefig("Plots_Prob/Prob_plot"+GWname+".png")
    
    
    
    # =========================================================================
    #  Confidence Region of sky localization corresponding to CL
    # =========================================================================
    
    # First obtain the minimum value of prob array above the chosen CL
    ConfReg_minimum_prob = get_ConfReg_minimum_prob(prob,CL)
    
    '''
    Contar el numero de elementos que satisfacen una condicion se puede hacer con
    np.sum tal como se hace aquí:
    
    si tenemos a = np.array([0,1,2,3,4])
    
    np.sum( a < 2 ) dara como resultado 2, pues solo hay dos elementos en el
    array tales que sean menores que 2.
    
    Luego yo veo que npix_ConfReg es el numero de pixeles cuya probabilidad es
    mayor que ConfReg_minimum_prob
    
    '''
    
    # Count how many pixels are above the given minimum prob corresponding to the chosen CL
    npix_ConfReg = np.sum( prob > ConfReg_minimum_prob )
    
    # Obtain all thetas, phis of the pixels in all the map
    theta, phi = hp.pix2ang(nside, np.arange(npix))
    
    # Convert thetas, phis of all the pixels to RA, Dec.
    radecs = astropy.coordinates.SkyCoord( ra=phi*u.rad, dec= (0.5*np.pi - theta) * u.rad )
    
    
    # =========================================================================
    # Obtaining the Contour
    # =========================================================================
    
    # First we obtain the points of the CL region
    
    prob_above_min = np.array([0 if probi <= ConfReg_minimum_prob else probi for probi in prob])
    
    # Now that we have all the points in the 90% CL region, we can construct the
    # contour points as the points that change from zero to nonzero values,
    # We remenber that in the healpix pixelization Schemes, we go for each 
    # declination, every Right ascention, and both increasing.
    
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
    
    # Now we create the contour arrays with just the declinations and right 
    # ascention of non void points
    
    dec_contour = np.array([ 0.5*np.pi - thetai for thetai, contouri in zip(theta,contour) if contouri != 0])
    ar_contour   = np.array([ phii for phii, contouri in zip(phi,contour) if contouri != 0]) 
    
    
    # =========================================================================
    # Preparing the contour for the plot
    # =========================================================================
    
    # Now we apply the functions, in order to obtain the lobes
    
    dec_contour_new, ar_contour_new = Ordering(dec_contour, ar_contour)
    lobes_dec, lobes_ar = Separate_Lobes(dec_contour_new, ar_contour_new, parameter_lobes)
    
    
    # =========================================================================
    # Plane contours plot (without the Mollweide projection), here we plot 
    # all the contour
    # =========================================================================
    '''
    mpl.rc('lines', linewidth=2.0)
    
    
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
    # Plot of the contour mollweide (different colours for the contours)
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
    ax.grid(True)
    ax.set_xticklabels(['2h','4h','6h','8h','10h','12h','14h','16h','18h','20h','22h'])
    
    cmap = plt.get_cmap('gnuplot')
    c = [cmap(i) for i in np.linspace(0, 1, len(lobes_ar))]
    for i,j,ci in zip(lobes_ar, lobes_dec, c): 
    
        plt.plot(np.array(i) - np.pi, np.array(j), "-", label = "CL contour %s" %n, color = ci)
        n = n + 1
        
    
    
    plt.legend(bbox_to_anchor=(1.12, 1.18), loc='upper right', borderaxespad=1)
    
    #plt.scatter(ar_contour - np.pi, dec_contour, s = 0.005)
    
    
    try:
        plt.savefig("Plots_Contour/ContourLobes"+GWname+".png")
    
    except:
        os.makedirs("Plots_Contour")
        plt.savefig("Plots_Contour/ContourLobes"+GWname+".png")
    
    
    # =============================================================================
    # Now we save the data
    # =============================================================================
    
    # Delete old data 
    
    file_dec = open("Data_contour/decContour"+GWname+".txt", "w")
    file_ar  = open("Data_contour/arContour"+GWname+".txt", "w")
    file_dec.close()
    file_ar.close()
    
    #os.remove("Data_contour/decContour"+GWname+".txt")
    #os.remove("Data_contour/arContour"+GWname+".txt")
        
    # Introduce new data
    
    file_dec = open("Data_contour/decContour"+GWname+".txt", "a")
    file_ar  = open("Data_contour/arContour"+GWname+".txt", "a")
    
    # We define this litle function in order to do it easier
    
    def write(a, file):
        
        for i in a:
            
            file.write(str(i))
            file.write(" ")
            
    
    for lobe_dec, lobe_ar in zip(lobes_dec,lobes_ar):
        
        # We save the declination and the Right ascencion for each lobe
        
        write(lobe_dec, file_dec)
        file_dec.write("\n")
        
        write(lobe_ar, file_ar)
        file_ar.write("\n")
        
    file_dec.close()
    file_ar.close()
        
    #np.savetxt("Pruebas/arContourGW"+GWname+".txt", dec_contour_new)
    #np.savetxt("Pruebas/decContourGW"+GWname+".txt", ar_contour_new )
    
    # Now if the script have reached this point it means that we want the GW-contour plot
    # so we launch fov_Auger_GW_builder.py in each case
    
    if GCN_mode == "yes":
        process = subprocess.run([ 'python3', 'fov_Auger_GW_builder.py', 'yes', FILENAME, GCN_ID, parameters])
    
    else:
        process = subprocess.run([ 'python3', 'fov_Auger_GW_builder.py', 'no', FILENAME, Scripts, parameters])

#Scriptfinaltime = time.time()

#print(Scriptfinaltime-scriptStartTime)