# -*- coding: utf-8 -*-
"""


The Script computes the integral of Exposures*E^(-2) in energy for the all channels
channels (DGL, DGH and ES and CC CN for differnts neutrino types) and then calculates 
the limits and combines them.

"""
import matplotlib.pyplot as plt
import numpy as np
import os 
import sys
from scipy.interpolate import interp1d
from scipy.integrate import quad
import healpy as hp
import astropy.time
from numba import njit

#print(os.pid)

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


@njit
def CLlim_Prob_dec(prob,dec,CL,prob_above_min):
    '''
    This function is just part of the code optimized with njit, takes the 
    probabilities, the declinations associated to that probabilities, the CL
    the prob_above_min (which has 0 value out of the CL region, and prob inside it)
    
    returns extremos_CL_dec, which is an array containing the declination extremes
    of the CL region (the declinations where the region starts or finishes. Returns 
    aux_prob, which is an array that adds the probability for each declination 
    in RA and aux_dec is an array which contains all the different  declinations

    '''
    
    npix = len(prob)
    
    # We are going to add all probabilities for each declination, 

    aux_dec         = [] 
    aux_prob        = []
    
    
    aux_dec.append(dec[0]) 
    aux      = 0
    
    for i in range(npix):
        
        if dec[i] != aux_dec[-1]: 
            
            # Once this happen, we have to change declination, so 
            # we add the results to the lists
            
            aux_dec.append(dec[i])
            aux_prob.append(aux)
            
            aux = 0
            
        else : 
            
            # In aux we store the sum of probabilities
            
            aux = aux + prob_above_min[i]
    
    aux_prob.append(aux) # In order to add the last sum
    
    
    # Here we store the Limits in declination for the 90% CL region
      
    extremos_CL_dec = []
    
    # If the first declination (north pole) starts with non-zero value, then we store 
    # its declination
    
    if aux_prob[0] != 0 : 
        extremos_CL_dec.append(aux_dec[0])
        
    for i in range(len(aux_dec)-1):
        
        # Basically if we change from 0 to 1 or vice versa, we are in a declination limit
        
        if (aux_prob[i] == 0) and (aux_prob[i+1] != 0):
            extremos_CL_dec.append(aux_dec[i+1])
        elif (aux_prob[i] != 0) and (aux_prob[i+1] == 0):
            extremos_CL_dec.append(aux_dec[i])
    
    # If the last declination (south pole) starts with non-zero value, then we store 
    # its declination
    
    if aux_prob[-1] != 0 : 
        extremos_CL_dec.append(aux_dec[-1])
        
    return extremos_CL_dec, aux_prob, aux_dec

# =============================================================================
# input GCN
# =============================================================================



if(len(sys.argv) < 2) :
    print('You have not include a fits file in the input')
    sys.exit()

else:
    FILENAME = sys.argv[2]
    GCN_mode = sys.argv[1]
    
    if GCN_mode == 'yes':
        GCN_ID   = sys.argv[3]
    

'''

if(len(sys.argv) < 2) :
    print('You have not include declination')
    sys.exit()

else:
    declination = float(sys.argv[1])
'''

# =============================================================================
# input files
# =============================================================================
'''
script_dir = os.getcwd()

path = os.path.join(script_dir, "Data_FITS")

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
'''

# =============================================================================
# Now read the Fluence data
# =============================================================================

cwd = os.path.dirname(os.path.realpath(__file__))

path_ES  = cwd + "/Data_exposures/Fluence_ES.txt"
path_DGH = cwd + "/Data_exposures/Fluence_DGH.txt"
path_DGL = cwd + "/Data_exposures/Fluence_DGL.txt"
path     = cwd + "/Data_exposures/Fluence_total.txt"


with open(path_ES) as file:
    
    filas_ES  = file.readlines()
    
with open(path_DGH) as file:    
    filas_DGH = file.readlines()

with open(path_DGL) as file:    
    filas_DGL = file.readlines()
    
with open(path) as file:    
    filas     = file.readlines()
    
k_ES  = []
k_DGH = []
k_DGL = []
k     = []

    
for fila in filas_ES:
    
    fila = fila.split()
    
    k_ES.append(fila[1])
    
for fila in filas_DGH:
    
    fila = fila.split()
    
    k_DGH.append(fila[1])

for fila in filas_DGL:
    
    fila = fila.split()
    
    k_DGL.append(fila[1])
    
for fila in filas:
    
    fila = fila.split()
    
    k.append(fila[1])

# Arrays with declinations

angulos_ES  = np.arange(-55,63) # the last number is not included in arange
angulos_DGH = np.arange(-69,56)
angulos_DGL = np.arange(-86,42)
angulos     = np.arange(-86,63)
        
k_ES = np.array(k_ES) ; k_DGH = np.array(k_DGH) ; k_DGL = np.array(k_DGL)

# limits in k interpolated

k_function  = interp1d(angulos,k)
k_function_ES  = interp1d(angulos_ES,k_ES)
k_function_DGH = interp1d(angulos_DGH,k_DGH)
k_function_DGL = interp1d(angulos_DGL,k_DGL)



def k_fun(x):
    if x >= -86 and x <= 62:
        return k_function(x)
    else:
        return 10**5 # in order to define limits in every declination

# =============================================================================
#%% Apart from this we have to compute the probability as a function of declination
# =============================================================================

angulos_x_ES  = np.linspace(-55,62, 1000)
angulos_x_DGH = np.linspace(-69,55, 1000)
angulos_x_DGL = np.linspace(-86,41, 1000)
angulos_x     = np.linspace(-86,62, 1000)

SKY = 4 * 180**2 / np.pi
CL = 0.9


    
    
prob, header = hp.read_map( FILENAME, h=True)

npix = len(prob)

area_per_pix, nside = SKY / npix, hp.npix2nside(npix)

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
    
    GWname = "GW" + str(GWtime.datetime)[2:4] + str(GWtime.datetime)[5:7] + str(GWtime.datetime)[8:10] + "_" + str(GWtime.datetime)[11:13]+str(GWtime.datetime)[14:16]

    if GCN_mode == 'yes' : 
        GWname = GCN_ID
else:               
    print('time_index invalid:',time_index,'could not find the MJD from the header file')
    sys.exit()

# =============================================================================
# Calculation of the CL extremes and probability function
# =============================================================================

# Arrays containing right ascencion and declination

theta, phi = hp.pix2ang(nside, np.arange(npix))    
dec = (0.5*np.pi - theta)*180/np.pi

# get the declination of the maximun value (for the limit flux)
Maximun = get_ConfReg_minimum_prob(prob, 0)
prob_dec = [(probi, (0.5*np.pi - thetai)*180/np.pi) for probi,thetai in zip(prob,theta)]
declination = [prob_deci[1] for prob_deci in prob_dec if prob_deci[0] == Maximun][0]


# obtain the 90% CL region
ConfReg_minimum_prob = get_ConfReg_minimum_prob(prob,CL)
prob_above_min = np.array([0 if probi <= ConfReg_minimum_prob else probi for probi in prob])


#print("Calculation for "+GWname)

prob = prob.astype('float64') # We do this for the @njit to work.

# Part of the code optimized with njit

extremos_CL_dec, probabilities, declinations = CLlim_Prob_dec(prob,dec,CL,prob_above_min)
    
#print(extremos_CL_dec)


# I'm going to change the maximun declination to be the point that adds more 
# declination probability and not the maximun probability point

declination = declinations[np.argmax(probabilities)]


# Now we compute the average limit

ks = []

for dec in declinations:
    ks.append(k_fun(np.array(dec)))

ks = np.array(ks)

inv_k_av = ks**(-1.)@np.array(probabilities)/sum(probabilities)
k_av = 1/inv_k_av


# We print the average and maximum limit values to pass it to CL_Coverage,
# but in the case of the maximum probability declination of the GW event being 
# out of the FOV of Auger, we don't include this information. The same for 
# the average value. If it is greater than 10**3, we don't include it

if declination >= min(angulos_x) and declination <= max(angulos_x):
    print(k_fun(declination)) 
    print(k_av)
    out_ind = 'no'
    out_av = 'no'
    
else:
    out_ind = 'yes'
    #print(k_fun(declination)) 
    print(out_ind)  
    
    if k_av > 10**3:
        out_av = 'yes'
        print(out_av)
        
    else:
        out_av = 'no'
        print(k_av)
        
    
# =============================================================================
#  Plot
# =============================================================================

top_offset = .04
left_offset = .10
right_offset = .23
bottom_offset = .10
hgap = .05
ax_width = 1-left_offset - right_offset
top_ax_height = (1-top_offset - bottom_offset - hgap)*3/4
bottom_ax_height = (1-top_offset - bottom_offset - hgap)*1/4



fig = plt.figure(FILENAME, figsize = [7,5])

ax = fig.add_axes([left_offset, bottom_offset + bottom_ax_height + hgap, ax_width, top_ax_height])

plt.semilogy(angulos_x_ES,k_function_ES(angulos_x_ES),"r-",label = "ES $\\left( 95^{o} > \\theta > 90^{o} \\right)$")
plt.semilogy(angulos_x_DGH,k_function_DGH(angulos_x_DGH),"b-", label = "DGH $\\left( 90^{o} > \\theta > 75^{o} \\right)$")
plt.semilogy(angulos_x_DGL,k_function_DGL(angulos_x_DGL),"g-", label = "DGL $\\left( 75^{o} > \\theta > 60^{o} \\right)$") 
plt.semilogy(angulos_x,k_function(angulos_x),"k-", label = "DGL+DGH+ES")
#plt.plot(angulos_ES,k_ES,"b.")
#plt.plot(angulos_DGH,k_DGH,"b.")
plt.grid()

plt.ylabel(" Fluence $\\left[GeV \cdot cm^{-2} \\right]$")

plt.ylim(1,10**3.5)
plt.xlim(-90,90)

#print(declination, float(k_fun(declination)))

# Resulta que los limites en y se definen entre 0 y 1 siendo 0 el minimo y 1 el máximo del plot,
# divido entre 3.5 porque son los ordenes de magnitud del plot en log

plt.axvline(declination, 0, np.log10(float(k_fun(declination)))/(3.5), color="darkcyan", linestyle="--" )

if out_ind == 'no' and out_av == 'no':

    plt.plot(declination, float(k_fun(declination)), color="darkcyan", marker= 'o', label = '$\\mathcal{F}_{GW}^{90} = $' + str(round(float(k_fun(declination)), 2)) + ' $GeV \cdot cm^{-2} $ \n$\\mathcal{F}_{GW}^{90, av} = $' + str(round(float(k_av), 2)) + ' $GeV \cdot cm^{-2} $')

elif out_ind == 'yes' and out_av == 'no':
    
    plt.plot(declination, float(k_fun(declination)), color="darkcyan", marker= 'o', label = '$\\mathcal{F}_{GW}^{90, av} = $' + str(round(float(k_av), 2)) + ' $GeV \cdot cm^{-2} $')

elif out_ind == 'no' and out_av == 'yes':

    plt.plot(declination, float(k_fun(declination)), color="darkcyan", marker= 'o', label = '$\\mathcal{F}_{GW}^{90} = $' + str(round(float(k_fun(declination)), 2)) + ' $GeV \cdot cm^{-2} $ ')
    
else:
    
    plt.plot(declination, float(k_fun(declination)), color="darkcyan", marker= 'o')
             
plt.xticks(np.arange(-90,105,15))
#ax.spines['bottom'].set_color('white')
ax.tick_params(axis='x', colors='white')

for i in range(0,len(extremos_CL_dec)-1,2):
    if i == 0:
        ax.axvspan(extremos_CL_dec[i], extremos_CL_dec[i+1], alpha=0.2, color='gray', label = "90$\\%$ CL region")
    else:
         ax.axvspan(extremos_CL_dec[i], extremos_CL_dec[i+1], alpha=0.2, color='gray')
         

plt.legend(bbox_to_anchor=(0.83, 1.1), loc='upper left', borderaxespad=1,fancybox=True,shadow=True)

ax_sub = fig.add_axes([left_offset, bottom_offset, ax_width, bottom_ax_height])

plt.xticks(np.arange(-90,105,15))
#plt.yticks([0.005,0.010,0.015])
#plt.ylim(0,0.017)
plt.xlim(-90,90)

plt.yticks([0])

plt.plot(declinations, np.array(probabilities), "k-")
plt.grid()

plt.ylabel("relative \n probability")
plt.xlabel("$\\delta$ $\\left[{}^{o}\\right]$")
plt.title(GWname)

for i in range(0,len(extremos_CL_dec)-1,2):
    ax_sub.axvspan(extremos_CL_dec[i], extremos_CL_dec[i+1], alpha=0.2, color='gray')
    

plt.axvline(declination, max(probabilities)/(max(probabilities) - min(probabilities)) , 1, color="darkcyan", linestyle="--" )
plt.plot(declination, max(probabilities), color="darkcyan", marker= '.')

plt.tight_layout()


#plt.show()
plt.savefig("Plots_Flux_limit/"+GWname+".pdf", format = 'pdf')


