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

'''
Program with the aim of calculating the field of view of Auger for a given
time of an event, written in the dictionary dicGW

'''

def plot_fov(datos1, datos2, datos3, datos4):
    
    # Now we will plot the part of the sky that we can see at the moment of the GW

    print('Plotting Auger fov at the moment of merger '+nameGW)

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
    
    hr2angle = 15.         # grados por hora
    deg2rad = np.pi/180.   # radianes por grado

    aux_label = nameGW
    
    # Recordamos que la AR fue guardada en horas y además lo habíamos ordenado
    # de menor a mayor, fijando el primer valor a 0 y el ultimo a 24. 
    # A1, A2, A3 y A4 no son exactamente iguales pero sí aproximadamente
    
    label_ES  = r'Auger fov: $\theta\in[90^{\rm o}, 95^{\rm o}]$ at UTC of '+aux_label    
    plt.fill_between((RA1*hr2angle-180.)*deg2rad,g1(RA1)*deg2rad,g2(RA1)*deg2rad,
    facecolor = "red",alpha=0.4,label=label_ES)
    
    label_DGH = r'Auger fov: $\theta\in[75^{\rm o}, 90^{\rm o}]$'
    plt.fill_between((RA1*hr2angle-180.)*deg2rad,g2(RA1)*deg2rad,g3(RA1)*deg2rad,
    facecolor = "blue",alpha=0.4,label=label_DGH)
#    plt.fill_between(RA1,g2(RA1),g3(RA1),facecolor="blue",alpha=0.4,label=label_DGH)

    label_DGL = r'Auger fov: $\theta\in[60^{\rm o}, 75^{\rm o}]$'
    plt.fill_between((RA1*hr2angle-180.)*deg2rad,g3(RA1)*deg2rad,g4(RA1)*deg2rad,
    facecolor = "green",alpha=0.4,label=label_DGL)
#    plt.fill_between(RA1,g3(RA1),g4(RA1),facecolor="green",alpha=0.4,label=label_DGL)


################################################################################
#~ # Geodetic coordinates of observatory (here: Pierre Auger Observatory)
alt_Auger=1400.
#--# Center of the array - see GAP 2001-038
lat_Auger=-35.235937 
lon_Auger=-69.249965
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

print('---------------------------------------------------------------------------')
print('Pierre Auger SD coordinates')
print('lat (deg), long (deg), height (m)')
print(lat_Auger,lon_Auger,alt_Auger)
observatory = astropy.coordinates.EarthLocation(lat=lat_Auger*u.deg, lon=lon_Auger*u.deg, height=alt_Auger*u.m)

#%%

################################################################################
#Library of GW events:
#~ # GW number and Time of the gravitational-wave event...

script_dir = os.getcwd()

path = os.path.join(script_dir, "FITS_files")

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
    
dicGW = {}

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


# Note: theta = 90 - alt_deg
alt_deg_array=(-5,0,15,30)
az_range = np.arange(0,360+0.5,0.5)

#    print(az_range[-1]) , with this definition this is 360.

# Loop over GW events
for j in dicGW:
    
    nameGW=j
    print('Calculation fov Auger for ',nameGW)
    time_GW = dicGW[nameGW]
    
    #~ # Alt/az reference frame at observatory, now
    Auger = astropy.coordinates.AltAz(obstime=time_GW, location=observatory)
    
    ################################################################################
    # Calculate field-of-view
    # Loop over azimuth for a fixed altitude 
    #(for instance 0 deg. corresponding to 90 deg. zenith angle)  
    
    result = []
    
    for alt_deg in alt_deg_array:

        out=[]
    
        for az_deg in az_range: 
            altaz = astropy.coordinates.SkyCoord(alt=alt_deg*u.degree, az=az_deg*u.degree, frame=Auger)
            radec = altaz.transform_to('icrs')
            
            # we store the declination, the AR and 
            
            out += [[radec.dec.deg,radec.ra.hour,90.-altaz.alt.deg]] # append
        
        out = np.array(out)
        
        # Ordenar salida por la segunda columna (RA)
        
        '''
        np.argsort devuelve las posiciones que ordenan el array
        a = np.array([2,4,1,3]) ; b = np.argsort(a)
        devuelve np.array([2,0,3,1]) tal que a[b] = np.array([1,2,3,4])
        
        '''
        
        out  = out[np.argsort(out[:,1])]
        
        #Redondeamos a 0 y a 24 el primer y ultimo ra
        out[0,1] = 0.
        out[len(az_range) - 1,1] = 24. # está bien pues el out[721] es para acimut 360,5 grados
        
        
        '''
        La sentencia exec esencialmente ejecuta el código escrito dentro
        '''
        
        if alt_deg == -5:
           # exec('np.savetxt("GW%s/GW%s_95deg.dat",out)'%(nameGW,nameGW))  # Put output in corresponding GW dir
            exec('np.savetxt("GW_Auger_fov/%s_95deg.dat",out)'%(nameGW))               # Put output in file in current dir
        if alt_deg == 0:
           # exec('np.savetxt("GW%s/GW%s_90deg.dat",out)'%(nameGW,nameGW))
            exec('np.savetxt("GW_Auger_fov/%s_90deg.dat",out)'%(nameGW))
        if alt_deg == 15:
           # exec('np.savetxt("GW%s/GW%s_75deg.dat",out)'%(nameGW,nameGW))
            exec('np.savetxt("GW_Auger_fov/%s_75deg.dat",out)'%(nameGW))
        if alt_deg == 30:
           # exec('np.savetxt("GW%s/GW%s_60deg.dat",out)'%(nameGW,nameGW))
            exec('np.savetxt("GW_Auger_fov/%s_60deg.dat",out)'%(nameGW))
            
        result.append(out)
    
# =============================================================================
# We read the contour data
# =============================================================================

    file_dec = open("Contour_Data/decContour"+nameGW+".txt", "r")
    file_ar  = open("Contour_Data/arContour"+nameGW+".txt", "r")
    
    
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
            lobes_dec.append(np.array(line_dec, dtype = float))
            lobes_ar.append(np.array(line_ar, dtype = float))

# =============================================================================
# Now we start the plot
# =============================================================================
    
    plt.clf()
    
    plt.rc("font",family="sans-serif",size=10)
    
    fig = plt.figure(1, figsize = [9,8])
        
    ax = fig.add_subplot(111, projection="mollweide")
    ax.set_xticklabels(['2h','4h','6h','8h','10h','12h','14h','16h','18h','20h','22h'])
    
    plt.text(0.1, 1.2, nameGW, horizontalalignment='center',
             verticalalignment='center', transform=plt.gca().transAxes, fontsize=14)
    
    # We plot the GW 90% CL Contours 
    
    plt.plot(lobes_ar[0]-np.pi, lobes_dec[0], "k-", 
             label = "$\\sim 90\\%$ CL contour "+nameGW+" (LIGO + VIRGO)")
    
    for i in range(len(lobes_dec)): 
        if i == 0 : continue
        plt.plot(np.array(lobes_ar[i]) - np.pi, np.array(lobes_dec[i]), "-k")
        
    # We plot the Auger fov
    
    plot_fov(result[0],result[1],result[2],result[3])
        
    handles, labels = plt.gca().get_legend_handles_labels()

        
    plt.gca().grid(linestyle='-.', linewidth=1)
    
    
    plt.legend(handles,labels,loc='upper right', bbox_to_anchor=(1,1.4), ncol=1, 
                   fancybox=True,shadow=False,numpoints=1)
    fig.tight_layout()

    plt.savefig("Fov_Contours/fov_Auger_{0}_mollweide.png".format(nameGW))


        
        
        


