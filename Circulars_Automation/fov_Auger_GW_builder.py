
'''
This program has the aim to obtain the plot of the contours at 90% Confidence
Level of the GW event, together with the Field of view of the Pierre Auger 
cosmic rays observatory, on a mollweide sky diagram, at the time of detection.

the inputs of the program are two. First the fits file of the event, that it is 
only used to obtain the exact time of detection. With that time, we calculate
the Auger FOV for neutrinos ES (-5º, 0º in altitude), DGH (0º, 15º), DGL
(15º, 30º) and we include the FOV to photons to (30º, 60º), but the program
is able to calculate this last FOV changing the last value between 30º and 90º
as we wish.
Then the second file is the result of Contour.py,a txt file which contains the 
data of calculated 90% CL contours.
'''


# Python standard library imports
import sys

# Third-party imports
import numpy as np
import astropy.coordinates
import healpy as hp
import astropy.time
import astropy.units as u
import os
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

import subprocess
import time

#print(os.getpid())

#T = 'yes'
#Mollweide = 'yes'
#GWname = 'GW150914_0950'
#GWtime = astropy.time.Time('2015-09-14 09:50:45.413937')
# =============================================================================
# Plot_fov function
# =============================================================================


def plot_fov(datos0, datos1, datos2, datos3, datos4, alt):
    '''
    This function plots the FOV of the Pierre Auger observatory, in equatorial 
    coordinates, for the three channels, DGH, DGL and ES, and in a region between 
    30º altitude and the alt parameter, given arrays with the position of
    the limits for these four regions. Usually alt = 60º, corresponding
    to the FOV for photons for completeness. 

    Parameters
    ----------
    datos0 : numpy.ndarray (3 columns)
        Array containing the data of a curve in equatorial coordinates, in this 
        case, it contains data from altitude 'alt', which can be chosen from 30º to
        90º as desired and it is passed as the fifth argument of the function. 
        The first column should be the declination in degrees
        The second column should be the RA in hours 
    datos1 : numpy.ndarray (3 columns)
        Array containing the data of a curve in equatorial coordinates, in this 
        case, it contains data from altitude 30º.
        The first column should be the declination in degrees
        The second column should be the RA in hours
    datos2 : numpy.ndarray (3 columns)
        Array containing the data of a curve in equatorial coordinates, in this 
        case, it contains data from altitude 15º.
        The first column should be the declination in degrees
        The second column should be the RA in hours
    datos3 : numpy.ndarray (3 columns)
        Array containing the data of a curve in equatorial coordinates, in this 
        case, it contains data from altitude 0º.
        The first column should be the declination in degrees
        The second column should be the RA in hours
    datos4 : numpy.ndarray (3 columns)
        Array containing the data of a curve in equatorial coordinates, in this 
        case, it contains data from altitude -5º.
        The first column should be the declination in degrees
        The second column should be the RA in hours
    alt : float
        The elevation of the upper limit to the completeness band to be plotted.
        Apart from DGH, DGL and ES, this function plots an extra band between
        30º and this parameter, in degrees

    Returns
    -------
    None
    '''
    
    # The latitude from Auger corresponds approximately with the value of the 
    # altitude 'alt' (in absolute value) which starts making data0 a closed curve.
    
    lat_Auger=-35.235937 
    
    # We start by reading RA and dec, from the data arrays, then we interpolate
    # the data.
    
    print('Plotting Auger fov at the moment of merger '+GWname)

    dec1 = datos1[:,0]
    RA1  = datos1[:,1]
    
    '''
    interp1d(x,y) devuelve una funcion f, tal que f(x) = y en los puntos indicados
    y continúa analíticamente entre para cualquier x en el intervalo de interpolacion
    uniendo con rectas.
    '''
    g1   = interp1d(RA1,dec1)
    ##

    dec2 = datos2[:,0]
    RA2  = datos2[:,1]
    g2   = interp1d(RA2,dec2)
    ##
    
    dec3 = datos3[:,0]
    RA3  = datos3[:,1]
    g3   = interp1d(RA3,dec3)
    ##

    dec4 = datos4[:,0]
    RA4  = datos4[:,1]
    g4   = interp1d(RA4,dec4)
    ##
    
    # We don't interpolate the alt curve, because it might be a closed curve
    # in fact it is for alt = 60º
    
    dec0 = datos0[:,0]
    RA0  = datos0[:,1]
    
    
    hr2angle = 15.         # degrees per hour
    deg2rad = np.pi/180.   # radians per degree

    aux_label = GWname
    
    # Recall that the RA was stored in hours and we had also sorted it
    # from smallest to largest, setting the first value to 0 and the last to 24.
    # A1, A2, A3 y A4 are not equal but are approximately equal.
    # We create labels, and we use fill_between to plot the FOV for the three
    # easy regions, DGH, DGL, ES.
    
    label_ES  = r'Auger fov: $\theta\in[90^{\rm o}, 95^{\rm o}]$ at UTC of '+aux_label    
    plt.fill_between((RA1*hr2angle-180.)*deg2rad,g3(RA1)*deg2rad,g4(RA1)*deg2rad,
    facecolor = "red",alpha=0.4,label=label_ES)
    
    label_DGH = r'Auger fov: $\theta\in[75^{\rm o}, 90^{\rm o}]$'
    plt.fill_between((RA1*hr2angle-180.)*deg2rad,g2(RA1)*deg2rad,g3(RA1)*deg2rad,
    facecolor = "blue",alpha=0.4,label=label_DGH)
#    plt.fill_between(RA1,g2(RA1),g3(RA1),facecolor="blue",alpha=0.4,label=label_DGH)

    label_DGL = r'Auger fov: $\theta\in[60^{\rm o}, 75^{\rm o}]$'
    plt.fill_between((RA1*hr2angle-180.)*deg2rad,g1(RA1)*deg2rad,g2(RA1)*deg2rad,
    facecolor = "green",alpha=0.4,label=label_DGL)
#    plt.fill_between(RA1,g3(RA1),g4(RA1),facecolor="green",alpha=0.4,label=label_DGL)

    
    label_DG = r'Auger fov: $\theta\in[30^{\rm o}, 60^{\rm o}]$'
    
    
    # When the altitude is greater than the latitude of Auger, then the contour
    # is closed and the plot is a bit harder to do
    
    if alt < abs(lat_Auger):
        
        # In this case, we can plot in the same way as we did before
        
        g0   = interp1d(RA0,dec0)
        
        
        label_DGL = r'Auger fov: $\theta\in[60^{\rm o}, 75^{\rm o}]$'
        plt.fill_between((RA1*hr2angle-180.)*deg2rad,g0(RA1)*deg2rad,g1(RA1)*deg2rad,
        facecolor = "gray",alpha=0.4,label=label_DG)
        
    
    else:
        
        # In order to treat this case well, we first define two functions
        
        def plot_out(x, label = 'no'):
            '''
            This function plots the FOV between the 30º curve and the 
            south celestial pole. It is useful to plot out of the closed curve
            at declination decl
            
            Parameters
            ----------
            x : np.ndarray
                RA values of the points to plot in radians
            label : str
                Choose any value to not label and 'yes' to do so. 
                The default is 'no'.
            Returns
            -------
            None.
            '''
                        
            if label == 'yes':
            
                plt.fill_between(x , -np.pi/2*np.ones(len(x)), g1(((x/deg2rad)+180.)/hr2angle)*deg2rad,
                             facecolor = "gray",alpha=0.4,label=label_DG)
            else: 
                plt.fill_between(x , -np.pi/2*np.ones(len(x)), g1(((x/deg2rad)+180.)/hr2angle)*deg2rad,
                             facecolor = "gray",alpha=0.4)
            
        def plot_in(x, y, separator):
            '''
            This function firstly plots the FOV between the 30º and the superior
            part of the closed curve for alt altitude, and then plots the FOV
            below the inferior part of the closed curve for alt altitude and the
            south celestial pole

            Parameters
            ----------
            x : np.ndarray
                RA values of the points to plot in radians
            y : np.ndarray
                dec values of the points to plot in radians
            separator : float
                RA value of the point used to separate the upper and lower regions
                of the curve

            Returns
            -------
            None.

            '''
        
            y_separator = y[np.where(x == separator)[0][0]]
            
            x_in_up = np.array([xi for xi,yi in zip(x,y) if yi >= y_separator])
            y_in_up = np.array([yi for xi,yi in zip(x,y) if yi >= y_separator])
            
            x_in_down = np.array([xi for xi,yi in zip(x,y) if yi <= y_separator])
            y_in_down = np.array([yi for xi,yi in zip(x,y) if yi <= y_separator])
            
            # we need to sort this arrays from the minor x to the highest
            
            xy_in_up = [(xi,yi) for xi,yi in zip(x_in_up,y_in_up)]
            xy_in_down = [(xi,yi) for xi,yi in zip(x_in_down,y_in_down)]
            
            xy_up_s   = sorted(xy_in_up, key = lambda x : x[0])
            xy_down_s = sorted(xy_in_down, key = lambda x : x[0])
            
            x_in_up   = np.array([xy[0] for xy in xy_up_s ])
            x_in_down = np.array([xy[0] for xy in xy_down_s])
            
            y_in_up   = np.array([xy[1] for xy in xy_up_s ])
            y_in_down = np.array([xy[1] for xy in xy_down_s ])
            
            gup     = interp1d(x_in_up,y_in_up)
            gdown   = interp1d(x_in_down,y_in_down)
            
            plt.fill_between(x_in_up , gup(x_in_up), g1(((x_in_up/deg2rad)+180.)/hr2angle)*deg2rad,
                             facecolor = "gray",alpha=0.4)
            plt.fill_between(x_in_down , -np.pi/2*np.ones(len(x_in_down)), gdown(x_in_down),
                             facecolor = "gray",alpha=0.4)
            
        # We change to radians, because we have to plot in radians for Mollweide
        # projection in matplotlib
        
        x = (RA0*hr2angle-180.)*deg2rad ; y = dec0*deg2rad
        
        
        # We find the minimun and the maximun in the alt curve
        
        minimo = np.min(x) ; maximo = np.max(x)
        
        # One thing that might happen is that the closed contour is splitted 
        # in two parts. That could happen because it basically starts on, let's
        # say in 22 hours, but it finishes in 3 hours in RA.
        #
        # to differentiate one situation from the other, we look to the borders
        # in right ascencion. If we have points, very close to the border, let's say
        #  1/100*2*np.pi, from the border. We consider that the curve is splitted
        
        
        if abs(maximo - np.pi) > 1/100*2*np.pi or abs(minimo + np.pi) > 1/100*2*np.pi :  
            
            # In this case, the curve is not splitted in two, so minimo and maximo
            # are basically the extremes in RA of the curve. 
            
            x_out = (RA1*hr2angle-180.)*deg2rad
            
            # if the points have RA higher than maximo, or lower than minimo, 
            # then we are don't have curve
            
            x_out_1 = x_out[x_out < minimo]  
            x_out_2 = x_out[x_out > maximo] 
            
            # We plot with the label only once, otherwise we would have two labels.
            
            plot_out(x_out_1, label = "yes") ; plot_out(x_out_2)
            
            # The point where the upper and lower curves are splitted can be chosen 
            # either maximo or minimo
            
            plot_in(x,y, minimo)
        
        else:
            
            # In this case, the curve is splitted in two regions, minimo and 
            # maximo doesn't have any useful information
            
            x_out = (RA1*hr2angle-180.)*deg2rad
            
            # We want to find the borders of the curves, the equivalent points
            # to the maximo and minimo in the one region case. To do so we first
            # sort the RA coordinates
            
            x_sorted = np.sort(x)
            
            # We compute the distance between consecutive points and we average
            # that distance
            
            distances = np.array([x_sorted[i+1]- x_sorted[i] for i in range(len(x) - 1)])
            a = np.average(distances)
            
            # We expect that only the separation between the splitted regions 
            # would be higher than this average, so separation should contain
            # just one value, and we compute the index
            
            separation = distances[distances > a][0] 
            separation_index = np.where(distances == separation)[0][0]
            
            # Now the equivalents to maximo and minimo would be
            
            maximo = x_sorted[separation_index]
            minimo = x_sorted[separation_index + 1]
            
            # The region between maximo and minimo wouldn't have decl curve
            
            x_out = x_out[x_out > maximo ] ; x_out = x_out[x_out < minimo]
                    
            plot_out(x_out, label = 'yes')
            
            # Now we use plot_in to both splitted regions
            
            xy_1 = np.array([(xi,yi) for xi, yi in zip(x, y) if xi <= maximo])
            xy_2 = np.array([(xi,yi) for xi, yi in zip(x, y) if xi >= minimo])
            
            x_1 = np.array([i[0] for i in xy_1])
            y_1 = np.array([i[1] for i in xy_1])
            
            x_2 = np.array([i[0] for i in xy_2])
            y_2 = np.array([i[1] for i in xy_2])
            
            plot_in(x_1,y_1, maximo)
            plot_in(x_2,y_2, minimo)
        
        
    
        
    
# =============================================================================
#%% Auger Coordinates
# =============================================================================

#~ # Geodetic coordinates of observatory (here: Pierre Auger Observatory)
alt_Auger = 1400.
#--# Center of the array - see GAP 2001-038
lat_Auger = -35.235937 
lon_Auger = -69.249965
#--# Max. Northing point (North)
#lat_Auger=-34.91607 
#lon_Auger=-69.055999
#--# Min. Northing point (South) 
#lat_Auger=-35.478049 
#lon_Auger=-69.205341
#--# Max. Easting point (East)
#lat_Auger=-35.126829 
#lon_Auger=-68.90806
#--# Min. Easting point (West)
#lat_Auger=-35.218812 
#lon_Auger=-69.64961

#print('---------------------------------------------------------------------------')
#print('Pierre Auger SD coordinates')
#print('lat (deg), long (deg), height (m)')
#print(lat_Auger,lon_Auger,alt_Auger)
observatory = astropy.coordinates.EarthLocation(lat=lat_Auger*u.deg, lon=lon_Auger*u.deg, height=alt_Auger*u.m)



# Dictionary that contains the files that will be analyzed (old code)

#dicGW = {}

# =============================================================================
#%% Directory (only useful if we want to run over /Data_FITS)
# =============================================================================

'''

################################################################################
#Library of GW events:
#~ # GW number and Time of the gravitational-wave event...

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
    




for i in list_path_file:
    
    prob, header = hp.read_map( i, h=True, verbose=False )
    
    timeindex = -1
    for i in range(len(header)):
        try:
            header[i].index('MJD-OBS') 
            time_index = i
    
        except:
            pass
    
    if time_index > -1 : 
        GWtime = astropy.time.Time( header[time_index][1], format='mjd' )
        
        
    
        # We construct the GWname from the date:
    
        GWname = "GW" + str(GWtime.datetime)[2:4] + str(GWtime.datetime)[5:7] + str(GWtime.datetime)[8:10]
    else:               
        print('time_index invalid:',time_index,'could not find the MJD from the header file')
        sys.exit()
        
    dicGW[GWname] = GWtime    

'''




# =============================================================================
#%% Input information
# =============================================================================

if(len(sys.argv) < 2) :
    print('You have not include a fits file in the input')
    sys.exit()

else:
    GCN_mode = sys.argv[1]
    file = sys.argv[2] 
    parameters = sys.argv[4]
    
    SCAN_RESOLUTION = float(parameters.split()[0])
    
    T = parameters.split()[1]
    Mollweide = parameters.split()[2]
    
    
    if GCN_mode == "yes":
        
        GCN_ID  = sys.argv[3]
        Scripts = "FOV_GW CL_Coverage Limit_Flux email"
    
    else:
        
        Scripts = sys.argv[3]
        

# We read the FITS file in order to extract the date and the name
        
        
prob, header = hp.read_map( file, h=True)

timeindex = -1
for i in range(len(header)):
        try:
            header[i].index('MJD-OBS') 
            time_index = i
    
        except:
            pass
    
if time_index > -1 : 
        GWtime = astropy.time.Time( header[time_index][1], format='mjd' )
        
        # We construct the GWname from the date if we are not in GCN_mode
    
        GWname = "GW" + str(GWtime.datetime)[2:4] + str(GWtime.datetime)[5:7] + str(GWtime.datetime)[8:10] + "_" + str(GWtime.datetime)[11:13]+str(GWtime.datetime)[14:16]

        if GCN_mode == 'yes':
            
            GWname = GCN_ID
else:               
        print('time_index invalid:',time_index,'could not find the MJD from the header file')
        sys.exit()

#dicGW[GWname] = GWtime


# =============================================================================
#%% Calculation of data of FOV_Auger
# =============================================================================

#dicGW['GW190517'] = astropy.time.Time(58619.348211976736, format='mjd')

# alt_deg_array contains the four altitudes considered for the four curves
# az_range contains azhimuts for each altitude

# Note: theta = 90 - alt_deg
alt_deg_array=(-5,0,15,30,60) # the splitting limit is between 35,27 y 35,28 
az_range = np.arange(0,360+0.5,0.1)

#    print(az_range[-1]) , with this definition this is 360.

print('Calculation fov Auger for ',GWname)

    
#t_GS = time_GW.sidereal_time('apparent', 'greenwich')
    
#time_GW = astropy.time.Time( time_GW.mjd - t_GS.hour/24, format='mjd' )
    
    
#~ # Alt/az reference frame at observatory, now
Auger = astropy.coordinates.AltAz(obstime= GWtime, location=observatory)
    
# We start here the actual computation
    
result = []
    
for alt_deg in alt_deg_array:
        
        # altaz contains all the points with alt_deg, this means, that contains
        # the information of one curve
        
        altaz = astropy.coordinates.SkyCoord(alt=alt_deg*u.degree, az=az_range*u.degree, frame=Auger)
        
        # We change to equatorial coordinates
        
        radec = altaz.transform_to('icrs')
            
        if T == 'yes':
            
            # T was a parameter contained in the input. Here we change the x
            # coordinate to be RA - t_Sidereal_source
            
            t_GS = GWtime.sidereal_time('apparent', 'greenwich')
            hr_angle = [radeci.ra.hour - t_GS.hour for radeci in radec]
            
            # Of course we have to correct in the case that the RA are lower than 0
            
            hr_angle = [24.0 + hr_anglei if hr_anglei < 0 else hr_anglei for hr_anglei in hr_angle ]
            #if (hr_angle < 0): hr_angle = 24.0 + hr_angle

        
            # we store the declination, the RA and altitude   
            out = [[radeci.dec.deg,hr_anglei,90.-altazi.alt.deg] for radeci, altazi,hr_anglei in zip(radec, altaz,hr_angle)]
        
        else:
            # we store the declination, the RA and altitude 
            out = [[radeci.dec.deg,radeci.ra.hour,90.-altazi.alt.deg] for radeci, altazi in zip(radec, altaz)]
        
        out = np.array(out)
        
        
        '''
        np.argsort devuelve las posiciones que ordenan el array
        a = np.array([2,4,1,3]) ; b = np.argsort(a)
        devuelve np.array([2,0,3,1]) tal que a[b] = np.array([1,2,3,4])
        '''
        
        # We sort out with the RA points just in the case that the curve is not
        # splitted, otherwise doesn't matter
        
        if alt_deg <= abs(lat_Auger):
            
            out  = out[np.argsort(out[:,1])]
            
            # In the old Script this limits were fixed too, is not something critical
            
            out[0,1] = 0.
            out[len(az_range) - 1,1] = 24. # está bien pues el out[721] es para acimut 360,5 grados
        
        
        '''
        La sentencia exec esencialmente ejecuta el código escrito dentro
        '''
        
        # And finally we save the data for each curve
        
        if alt_deg == -5:
           # exec('np.savetxt("GW%s/GW%s_95deg.dat",out)'%(nameGW,nameGW))  # Put output in corresponding GW dir
            exec('np.savetxt("Data_FOV/%s_95deg.dat",out)'%(GWname))               # Put output in file in current dir
        if alt_deg == 0:
           # exec('np.savetxt("GW%s/GW%s_90deg.dat",out)'%(nameGW,nameGW))
            exec('np.savetxt("Data_FOV/%s_90deg.dat",out)'%(GWname))
        if alt_deg == 15:
           # exec('np.savetxt("GW%s/GW%s_75deg.dat",out)'%(nameGW,nameGW))
            exec('np.savetxt("Data_FOV/%s_75deg.dat",out)'%(GWname))
        if alt_deg == 30:
           # exec('np.savetxt("GW%s/GW%s_60deg.dat",out)'%(nameGW,nameGW))
            exec('np.savetxt("Data_FOV/%s_60deg.dat",out)'%(GWname))
        if alt_deg == 60:
           # exec('np.savetxt("GW%s/GW%s_60deg.dat",out)'%(nameGW,nameGW))
            exec('np.savetxt("Data_FOV/%s_30deg.dat",out)'%(GWname))
            
        result.append(out)
        
  
    
# =============================================================================
#%% We read the contour data
# =============================================================================

# read the files and store data

file_dec = open("Data_contour/decContour"+GWname+".txt", "r")
file_ar  = open("Data_contour/arContour"+GWname+".txt", "r")


lines_file_dec = file_dec.readlines()
lines_file_ar  = file_ar.readlines()

file_dec.close()
file_ar.close()

lobes_dec = []
lobes_ar  = []

for line_dec, line_ar in zip(lines_file_dec, lines_file_ar):
    
    line_dec = line_dec.split()
    line_ar  = line_ar.split()
    
    if len(line_ar) == 0 : continue
    if len(line_dec) == 0 : continue
    
    if line_dec[0] == "\n" : continue
    if line_ar[0] == "\n" : continue


    else:
        #store data in lobes lists
        lobes_dec.append(np.array(line_dec, dtype = float))
        lobes_ar.append(np.array(line_ar, dtype = float))

# =============================================================================
#%% Now we start the plot
# =============================================================================

plt.clf()

plt.rc("font",family="sans-serif",size=10)

fig = plt.figure(1, figsize = [9,8])

if Mollweide == 'yes':
    ax = fig.add_subplot((111), projection="mollweide")
    
    bbox_to_anchor = 1.4
    
else:
    ax = fig.add_subplot(111)
    
    ax.set_ylim(-np.pi/2, np.pi/2)
    ax.set_yticks(np.array([-90, -75, -60, -45, -30, -15, 0, 15, 30, 45, 60, 75, 90])*np.pi/180)
    ax.set_yticklabels([-90, -75, -60, -45, -30, -15, 0, 15, 30, 45, 60, 75, 90])
    
    ax.set_xlim(-np.pi,np.pi)
    ax.set_xticks(np.array([2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22])*np.pi/(180)*15 - np.pi)
    
    bbox_to_anchor = 1.25
    
# the plot will be done with data in radians, but this can be done in the x axis    

ax.set_xticklabels(['2h','4h','6h','8h','10h','12h','14h','16h','18h','20h','22h'])

plt.text(0.1, 1.2, GWname, horizontalalignment='center',
         verticalalignment='center', transform=plt.gca().transAxes, fontsize=14)

# We plot the GW 90% CL Contours 

if T == 'yes':
    
   
    # In this case it may happen that by moving the points, a lobe that was previously 
    # on the left becomes part on the right, so I separate the contour in two regions, 
    # to avoid the horizontal lines that join both parts.

    lobes_ar_right  = []
    lobes_dec_right = []
    
    lobes_ar_left   = []
    lobes_dec_left  = []
    
    lobes_ar_shifted = []
    
    for i in range(len(lobes_ar)):
        
        # we first shift the lobes
        
        lobes_ar_shifted.append(np.array([ra_j*180/np.pi - t_GS.hour*15 for ra_j in lobes_ar[i]]))
        
        # Now we separate the lobes into right and left parts depending whether
        # it is neccessary to add 360 or not

        lobes_ar_right.append(np.array([360 + ra_j  for ra_j in lobes_ar_shifted[i] if ra_j < 0])*np.pi/180)
        lobes_dec_right.append(np.array([dec_j  for ra_j , dec_j in zip(lobes_ar_shifted[i],lobes_dec[i]) if ra_j < 0]))
        
        lobes_dec_left.append(np.array([dec_j for ra_j , dec_j in zip(lobes_ar_shifted[i],lobes_dec[i]) if ra_j >= 0]))
        lobes_ar_left.append(np.array([ra_j for ra_j in lobes_ar_shifted[i] if ra_j >= 0])*np.pi/180)
        
    # Now we plot each lobe
    
    for i in range(len(lobes_dec_right)): 
        if len(lobes_dec_right[i]) == 0 : continue
        plt.plot(np.array(lobes_ar_right[i]) - np.pi, np.array(lobes_dec_right[i]), "-k")
    
    # This j variable is intended to plot just one lobe with the label
    # it ensures that the arrays chosen are not empty
    
    j = 0
    
    for i in range(len(lobes_dec_left)): 
        if i == j : 
            if len(lobes_dec_left[i]) != 0 :
                plt.plot(lobes_ar_left[j]-np.pi, lobes_dec_left[j], "k-", 
                         label = "$\\sim 90\\%$ CL contour "+GWname+" (LIGO + VIRGO)")
            else:
                j += 1
                
        if len(lobes_dec_left[i]) == 0 : continue
        plt.plot(np.array(lobes_ar_left[i]) - np.pi, np.array(lobes_dec_left[i]), "-k")
          
        
else:       
    plt.plot(lobes_ar[0]-np.pi, lobes_dec[0], "k-", 
             label = "$\\sim 90\\%$ CL contour "+GWname+" (LIGO + VIRGO)")
    
    for i in range(len(lobes_dec)): 
        if i == 0 : continue
        plt.plot(np.array(lobes_ar[i]) - np.pi, np.array(lobes_dec[i]), "-k")
        
#plt.show()
#sys.exit()
    
# We plot the Auger fov using the function defined to do so

x = plot_fov(result[4],result[3],result[2],result[1], result[0], alt_deg)
    
handles, labels = plt.gca().get_legend_handles_labels()

    
plt.gca().grid(linestyle='-.', linewidth=1)


plt.legend(handles,labels,loc='upper right', bbox_to_anchor=(1,bbox_to_anchor), ncol=1, 
               fancybox=True,shadow=False,numpoints=1)
fig.tight_layout()

# save the plot


plt.savefig("Plot_FOV_Contours/fov_Auger_{0}_mollweide.png".format(GWname))

if GCN_mode == "yes":
    
    
    #print('Executing CL_coverage.py for', GWname)
    subprocess.run(['python3', 'CL_coverage.py', 'yes', file, GCN_ID, str(SCAN_RESOLUTION)])

else:
    #print('Executing CL_coverage.py for', GWname)
    subprocess.run([ 'python3', 'CL_coverage.py','no', file, Scripts, str(SCAN_RESOLUTION)])




        
        
        


