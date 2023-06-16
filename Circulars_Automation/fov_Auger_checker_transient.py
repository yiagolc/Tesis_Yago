import sys
import numpy as np
import matplotlib.pyplot as plt

from scipy.interpolate import interp1d
import astropy.coordinates
import astropy.time
import astropy.units as u

import os

currentdir = os.getcwd()

#from astropy.utils.iers import iers_auto_url_mirror
#print('###### ',iers.iers_auto_url_mirror)

#from astroplan import download_IERS_A
#download_IERS_A()

################################################################################
def plot_position_source():

#-------------------------------------------------------------------------------
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
#--# Coordinates M. Schimp 
    #lat_Auger=-35.20666735
    #lon_Auger=-69.315833
#--# Original coordinates
    #lat_Auger=-35.15
    #lon_Auger=-69.2
   print('---------------------------------------------------------------------------')
   print('Pierre Auger SD coordinates')
   print('lat (deg), long (deg), height (m)')
   print(lat_Auger,lon_Auger,alt_Auger)

# Define location of Pierre Auger Observatory
   observatory = astropy.coordinates.EarthLocation(lat=lat_Auger*u.deg, lon=lon_Auger*u.deg, height=alt_Auger*u.m)

#-------------------------------------------------------------------------------
#   datos=open("GW170817.inp")
#   datos=open("GRB190114C.inp")
#   datos=open("TXS_0506+056.inp")
#   datos=open("GW170818_center.inp")
#   datos=open("IceCube_190331A.inp")
#   datos=open("IceCube_200926A.inp")
#   datos=open("GW190408_center_approx.inp")
#   datos=open("GW190412_center1_approx.inp")
#   datos=open("GW190412_center2_approx.inp")


   if len(sys.argv) > 1:
       FILENAME = sys.argv[1]
   else:
       print('USAGE EXAMPLE: python fov_Auger_checker_transient.py GW170817.inp')
       FILENAME = currentdir + "/Data_transient/GW170817.inp"

   datos = open(FILENAME)
   print('Reading...',FILENAME)

   for line in datos:           
        f=line.split()      # splits line in columns and puts values in array c
        if(f[0]=="NAME"):    name=f[1]
        if(f[0]=="RA_hr"):   ra_hr=float(f[1])
        if(f[0]=="RA_min"):  ra_min=float(f[1])
        if(f[0]=="RA_sec"):  ra_sec=float(f[1])
        if(f[0]=="DEC_deg"): dec_deg=float(f[1])
        if(f[0]=="DEC_min"): dec_min=float(f[1])
        if(f[0]=="DEC_sec"): dec_sec=float(f[1])
        if(f[0]=="UTC"):     utc=f[1]+' '+f[2]

    #print(name,ra_hr,ra_min,ra_sec,dec_deg,dec_min,dec_sec,utc)

   ra_degrees = ra_hr*15.+ra_min*15./60.+ra_sec*15./3600. 
   dec_degrees = dec_deg+dec_min/60.+dec_sec/3600. 

   print('---------------------------------------------------------------------------')
   print('Source name:')
   print(name)
   print('Source coordinates')
   print('ra (deg)  dec (deg)')
   print(ra_degrees,dec_degrees)

#-------------------------------------------------------------------------------
#~ # Time of the bursting event
   time_source = astropy.time.Time(utc, format='iso', scale='utc', location=observatory)
   
   frame = astropy.coordinates.AltAz(obstime=time_source, location=observatory)
   radec = astropy.coordinates.SkyCoord(ra=ra_degrees*u.degree, dec=dec_degrees*u.degree)
   altaz = radec.transform_to(frame)
   
   print('UTC source (hr,min,sec) ')
   print(time_source)
   print('---------------------------------------------------------------------------')
   print('alt (deg),zenith,azimuth (deg) as seen from Auger SD')
   print(altaz.alt.deg,90.-altaz.alt.deg,altaz.az.deg)
   print('NOTE: Azimuth is measured clockwise from North')
   #print('---------------------------------------------------------------------------')

#-------------------------------------------------------------------------------
# Convert UTC to Greenwich sidereal time

   '''
   En todos los plots que hemos realizado hasta ahora conforme el tiempo aumenta
   la Tierra rota y el FOV de Auger se desplaza hacia mayores valores de AR.
   
   Con el objetivo de mantener el FOV de Auger constante en todos los plots 
   que hagamos de aquí en adelante, tomaremos la resta de la AR y el tiempo
   sidéreo medido desde Greenwich
   '''

   t_GS = time_source.sidereal_time('apparent', 'greenwich')
   print('---------------------------------------------------------------------------')
   print('Greenwich sidereal time')
   print(t_GS)   

   print('Hour angle - Greenwich (hours)')
   hr_angle = ra_degrees/15. - t_GS.hour 
   if (hr_angle < 0): hr_angle = 24.0 + hr_angle
   print(hr_angle) 

# Vertical line at alpha-t_GS of source 
   x_source=hr_angle
   plt.axvline(x=x_source, ls='-.', color='black', linewidth=1.5)

# Horizontal line at declination of GW170817
   y_source=dec_degrees
   plt.axhline(y=y_source, ls='-.', color='black', linewidth=1.0)

# Arrow with text
   x=x_source
   y=y_source 
   
   print(x, y)
# Position label depending on size of the plot where the transient source is
   if (x <= 12.) and (y <= 0):
           plt.annotate(name, xy=(x, y), xytext=(x + 1.0, y + 10.0), color='black',
           arrowprops=dict(color='black', shrink=0.05), fontsize=18)
            #arrowprops=dict(facecolor='black', arrowstyle='-|>'),
            
   elif (x > 12.) and (y <= 0):
           print(x, y)
           plt.annotate(name, xy=(x, y), xytext=(x - 4.0, y + 10.0), color='black',
           arrowprops=dict(color='black', shrink=0.05), fontsize=18)
            #arrowprops=dict(facecolor='black', arrowstyle='-|>'),
   elif (x <= 12.) and (y > 0):
           plt.annotate(name, xy=(x, y), xytext=(x + 1.0, y - 10.0), color='black',
           arrowprops=dict(color='black', shrink=0.05), fontsize=18)
            #arrowprops=dict(facecolor='black', arrowstyle='-|>'),
   else:
           plt.annotate(name, xy=(x, y), xytext=(x - 4.0, y - 20.0), 
           color='black', arrowprops=dict(color='black',shrink=0.05), fontsize=18)
            #arrowprops=dict(facecolor='black', arrowstyle='-|>'),


#####################################################################
# Auger fov - Plot boundaries of fov corresponding to each channel 
# at t_GS = 0 (see JCAP 2019 paper for details)
# and fill between boundaries with color
# 
# fov files needed obtained with:
# /Users/jaime/Auger/Neutrinos/Multi-messenger/Sky_tools/gfortran/fov_Auger_UTC_GST_zero.make
#####################################################################
def plot_fov():
    
    datos4=np.loadtxt(currentdir + "/Data_FOV_transient/fov_95deg.dat")   

    datos3=np.loadtxt(currentdir + "/Data_FOV_transient/fov_90deg.dat")   

    datos2=np.loadtxt(currentdir + "/Data_FOV_transient/fov_75deg.dat")   

    datos1=np.loadtxt(currentdir + "/Data_FOV_transient/fov_60deg.dat")   
    
    datos0=np.loadtxt(currentdir + "/Data_FOV_transient/fov_30deg.dat")  

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
    dec0 = datos0[:,0]
    RA0  = datos0[:,1]

    
    # Recordamos que la AR fue guardada en horas y además lo habíamos ordenado
    # de menor a mayor, fijando el primer valor a 0 y el ultimo a 24. 
    # A1, A2, A3 y A4 no son exactamente iguales pero sí aproximadamente
    
    label_ES  = r'Auger fov: $\theta\in[90^{\rm o}, 95^{\rm o}]$'    
    plt.fill_between(RA1,g3(RA1),g4(RA1),
    facecolor = "red",alpha=0.4,label=label_ES)
    
    label_DGH = r'Auger fov: $\theta\in[75^{\rm o}, 90^{\rm o}]$'
    plt.fill_between(RA1,g2(RA1),g3(RA1),
    facecolor = "blue",alpha=0.4,label=label_DGH)
#    plt.fill_between(RA1,g2(RA1),g3(RA1),facecolor="blue",alpha=0.4,label=label_DGH)

    label_DGL = r'Auger fov: $\theta\in[60^{\rm o}, 75^{\rm o}]$'
    plt.fill_between(RA1,g1(RA1),g2(RA1),
    facecolor = "green",alpha=0.4,label=label_DGL)
#    plt.fill_between(RA1,g3(RA1),g4(RA1),facecolor="green",alpha=0.4,label=label_DGL)

    
    label_DG = r'Auger fov: $\theta\in[30^{\rm o}, 60^{\rm o}]$'
    
    
    # When the altitude is greater than the latitude of Auger, then the contour
    # is closed and the plot is a bit harder to do
    

    
    def plot_out(x, label = 'no'):
        
        '''Plot when there is no new contour'''
        
        if label == 'yes':
        
            plt.fill_between(x , -90*np.ones(len(x)), g1(x),
                         facecolor = "gray",alpha=0.4,label=label_DG)
        else: 
            plt.fill_between(x , -90*np.ones(len(x)), g1(x),
                         facecolor = "gray",alpha=0.4)
        
    def plot_in(x, y, separator):
        
        '''
        Plots FOV above and below a contour, and below de 60º contour
        
        x: contour , 
        separator: x point to separate the upper and the lower points
        of the contour
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
        
        plt.fill_between(x_in_up , gup(x_in_up), g1(x_in_up),
                         facecolor = "gray",alpha=0.4)
        plt.fill_between(x_in_down , -90*np.ones(len(x_in_down)), gdown(x_in_down),
                         facecolor = "gray",alpha=0.4)
        
    
    
    minimo = np.min(RA0) ; maximo = np.max(RA0)
    
    
            
    # Also the new closed contour can be separated into two parts
    # we decide this if the minimum and the maxima are in the border of
    # the plots.
    
     
        
        
    x_out = RA1
    
    x_out_1 = x_out[x_out < minimo]
    x_out_2 = x_out[x_out > maximo] 
        
    plot_out(x_out_1, label = "yes") ; plot_out(x_out_2)
        
    plot_in(RA0,dec0, minimo)
    
    

#####################################################################
#####################################################################
if __name__ == '__main__':

    plt.rc("font",family="sans-serif",size=14)
    fig = plt.figure(figsize=(12,8))


# Labels
#    plt.text(0.85, 0.05,'PRELIMINARY', horizontalalignment='center',
#             verticalalignment='center', transform=plt.gca().transAxes, fontsize=20)

#----------------------------------------------------------------------------------------------
    plot_fov()
    plot_position_source()
#----------------------------------------------------------------------------------------------

# Labels
    plt.xlabel(r"$\alpha - t_{\rm GS}$ (hr)",fontsize=18)
    plt.ylabel(r"Declination $\delta$ (deg)",fontsize=18)
    plt.tick_params(length=8, width=1, top=True, right=True, labelsize=18)

# Legends
    handles, labels = plt.gca().get_legend_handles_labels()
    plt.legend(handles,labels,loc='upper center', bbox_to_anchor=(0.5,0.995), ncol=2, fancybox=True, shadow=False,numpoints=1,framealpha=1)

# Axes: ticks, scale, lims,...
    xmin=0.
    xmax=24.
    plt.gca().set_xlim(xmin, xmax)
    plt.xticks(np.arange(xmin, xmax+1, 1))

    dec_min=-90.
    dec_max=90.
    plt.gca().set_ylim(dec_min, dec_max)
    plt.gca().grid(linestyle=':', linewidth=0.5)
    plt.yticks(np.arange(dec_min, dec_max+1, 15))


    fig.tight_layout()
    
    plt.savefig(currentdir + "/Plots_transient/fov_Auger_checker_transient.pdf")

#    plt.show()

