# -*- coding: utf-8 -*-
"""
Script that computes the one day flux limit for a puntual source, assuming a normal 
operating status for the SD of the Pierre Auger Observatory. It takes as input the
1 day average exposures of the detector (passed as a file), and the declination 
of the source (in the sys.argv).

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


def get_ConfReg_minimum_prob(prob,CL):
    
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
    
    returns extremos_CL_dec, which is an array containing the declination extremes of
    the CL region. Returns aux_prob, which is an array that adds the probability for
    each declination and aux_dec is an array which contains all the different declinations
    '''
    
    npix = len(prob)
    
    extremos_CL_dec = []
    aux_dec         = [] 
    aux_prob        = []
    
    
    aux_dec.append(dec[0]) 
    aux      = 0
    
    for i in range(npix):
        
        if dec[i] != aux_dec[-1]: 
            aux_dec.append(dec[i])
            aux_prob.append(aux)
            
            aux = 0
        else : 
            aux = aux + prob_above_min[i]
    
    aux_prob.append(aux) # hay que añadir la última suma
        
    
    if aux_prob[0] != 0 : 
        extremos_CL_dec.append(aux_dec[0])
        
    for i in range(len(aux_dec)-1):
        if (aux_prob[i] == 0) and (aux_prob[i+1] != 0):
            extremos_CL_dec.append(aux_dec[i+1])
        elif (aux_prob[i] != 0) and (aux_prob[i+1] == 0):
            extremos_CL_dec.append(aux_dec[i])
    
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
    list_path_file = [sys.argv[1]]

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

# =============================================================================
# Now compute the flux limit function
# =============================================================================

cwd = os.path.dirname(os.path.realpath(__file__))

path_ES  = cwd + "/Data_exposures/exposure_average1day_ES.dat"
path_DGH = cwd + "/Data_exposures/exposure_average1day_DGH.dat"


with open(path_ES) as file:
    
    filas_ES  = file.readlines()
    
with open(path_DGH) as file:    
    filas_DGH = file.readlines()
    
k_ES  = []
k_DGH = []
k_DGL = []

exposure = []
energy   = []

    
for fila in filas_ES:
    
    fila = fila.split()
    
    if (float(fila[1]) == 0 and float(fila[0]) == 0) == False:
    
        energy.append(float(fila[1]))
        exposure.append(float(fila[2]))
        
        angulo_act = float(fila[0])
    
    
    if float(fila[1]) == 0 and float(fila[0]) == 0:  # to only innclude separations
        
        energy = np.array(energy) ; exposure = np.array(exposure)

        f = interp1d(energy,exposure)
    
        def funcion(x):
            return (x**(-2))*f(x)
    
        Integration = quad(funcion, min(energy), max(energy),points = energy)
    
        if Integration[0] != 0:
    
            ki = 2.39/Integration[0]*24*3600
    
            k_ES.append(ki)
        else:
            k_ES.append(10**(5))
        
        if angulo_act == 55:
            
            #print("Energies")
            #print(" ")
            #print(energy)
            #print(" ")
            #print("Exposures for declination = ", angulo_act)
            #print(" ")
            #print(exposure)
            #print(" ")
            #print("flux limit for declination = ", angulo_act)
            #print(" ")
            #print(ki)
            
            x = np.linspace(min(energy),max(energy),100000)
        '''    
            plt.figure(2)
            plt.loglog(x*10**9,f(x),"b-")
            plt.loglog(energy*10**9, exposure,"r.")
            
            plt.xlabel("$E_{\\nu}$  $\\left[ eV \\right]$")
            plt.ylabel("$\\epsilon$  $\\left[ cm^2 \\cdot s \\right]$")
            
            plt.title("$\\delta ="+str(angulo_act)+"^{o}$")
            
            plt.figure(4)
            plt.loglog(x,f(x)*x**(-2),"b-")
            plt.loglog(energy, exposure*energy**(-2),"r.")
            
            plt.xlabel("$E_{\\nu}$  $\\left[ GeV \\right]$")
            plt.ylabel("$\\epsilon \\cdot E_{\\nu}^{-2}$  $\\left[ cm^2 \\cdot s \\cdot GeV^{-2}\\right]$")
            
            plt.title("$\\delta ="+str(angulo_act)+"^{o}$")
        '''
            
        exposure = []
        energy   = []
        
#%%

for fila in filas_DGH:
    
    fila = fila.split()
    
    if (float(fila[1]) == 0 and float(fila[0]) == 0) == False:
    
        energy.append(float(fila[1]))
        exposure.append([float(fila[2]),float(fila[3]),float(fila[4]),float(fila[5])])
        
        angulo_act = float(fila[0])
    
    
    if float(fila[1]) == 0 and float(fila[0]) == 0:  # to only include separations
        
        energy = np.array(energy) ; exposure = np.array(exposure)

        f_1 = interp1d(energy,exposure[:,0])
        f_2 = interp1d(energy,exposure[:,1])
        f_3 = interp1d(energy,exposure[:,2])
        f_4 = interp1d(energy,exposure[:,3])
        
    
        def funcion_1(x):
            return (x**(-2))*f_1(x)
        def funcion_2(x):
            return (x**(-2))*f_2(x)
        def funcion_3(x):
            return (x**(-2))*f_3(x)
        def funcion_4(x):
            return (x**(-2))*f_4(x)
    
        Integration_1 = quad(funcion_1, min(energy), max(energy),points = energy)
        Integration_2 = quad(funcion_2, min(energy), max(energy),points = energy)
        Integration_3 = quad(funcion_3, min(energy), max(energy),points = energy)
        Integration_4 = quad(funcion_4, min(energy), max(energy),points = energy)
        
        Integration = np.array([Integration_1[0],Integration_1[0],Integration_1[0],Integration_1[0]])
    
        if np.all(Integration != 0):
    
            k1 = 2.39/Integration_1[0]*24*3600
            k2 = 2.39/Integration_2[0]*24*3600
            k3 = 2.39/Integration_3[0]*24*3600
            k4 = 2.39/Integration_4[0]*24*3600
            
            ki = 1/(1/k1+1/k2+1/k3+1/k4)
            
            k_DGH.append(ki)
        else:
            k_DGH.append(10**(3))
        
        if angulo_act == 45:
            
            x = np.linspace(min(energy),max(energy),100000)
        '''    
            plt.figure(3)
            plt.loglog(x*10**9,f_1(x),"b-")
            plt.loglog(energy*10**9, exposure[:,0],"r.")
            
            plt.xlabel("$E_{\\nu}$  $\\left[ eV \\right]$")
            plt.ylabel("$\\epsilon$  $\\left[ cm^2 \\cdot s \\right]$")
            
            plt.title("$\\delta = "+str(angulo_act)+"^{o}$")
        '''
            
        exposure = []
        energy   = []


angulos_ES  = np.arange(-55,63)
angulos_DGH = np.arange(-69,56)
angulos     = np.arange(-69,63)
        
k_ES = np.array(k_ES) ; k_DGH = np.array(k_DGH)

kES = np.zeros(63) ; kDGH = np.zeros(63)

dif_ES  = 69 - 55 
dif_DGH = 63 - 56 


kES = np.array([10**5 if i < dif_ES else k_ES[i - dif_ES] for i in range(len(k_ES)+ dif_ES)])
kDGH = np.array([10**5 if i >= len(k_DGH) else k_DGH[i] for i in range(len(k_DGH)+ dif_DGH)])



k = 1/(1/kES+1/kDGH)

k_function  = interp1d(angulos,k)
k_function_ES  = interp1d(angulos_ES,k_ES)
k_function_DGH = interp1d(angulos_DGH,k_DGH)


#print(k_function(declination))


def k_fun(x):
    if x >= -55 and x <= 69:
        return k_function(x)
    else:
        return 10**3

# =============================================================================
#%% Now that we have the function we start the iteration over the files
# =============================================================================

angulos_x_ES  = np.linspace(-55,62, 1000)
angulos_x_DGH = np.linspace(-69,55, 1000)
angulos_x     = np.linspace(-69,62, 1000)

SKY = 4 * 180**2 / np.pi
CL = 0.9

for FILENAME in list_path_file:
    
    
    prob, header = hp.read_map( FILENAME, h=True, verbose=False )
    
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
    else:               
        print('time_index invalid:',time_index,'could not find the MJD from the header file')
        sys.exit()
    
# =============================================================================
# Calculation of the CL extremes and probability function
# =============================================================================
    
    theta, phi = hp.pix2ang(nside, np.arange(npix))    
    
    dec = (0.5*np.pi - theta)*180/np.pi
    
    # get the declination of the maximun value (for the limit flux)
    Maximun = get_ConfReg_minimum_prob(prob, 0)
    prob_dec = [(probi, (0.5*np.pi - thetai)*180/np.pi) for probi,thetai in zip(prob,theta)]
    declination = [prob_deci[1] for prob_deci in prob_dec if prob_deci[0] == Maximun][0]
    
    print(k_fun(declination)) #this value is passed to CL_Coverage in GCN_mode
    
    ConfReg_minimum_prob = get_ConfReg_minimum_prob(prob,CL)
    prob_above_min = np.array([0 if probi <= ConfReg_minimum_prob else probi for probi in prob])
    
    
    #print("Calculation for "+GWname)
    
    prob = prob.astype('float64') # We do this for the @njit to work.
    
    extremos_CL_dec, probabilities, declinations = CLlim_Prob_dec(prob,dec,CL,prob_above_min)
        
    #print(extremos_CL_dec)
    
# =============================================================================
#  Plot
# =============================================================================
    
    fig = plt.figure(FILENAME)
    
    ax = fig.add_axes([.15,0.35,.8,.6])
    
    plt.semilogy(angulos_x_ES,k_function_ES(angulos_x_ES),"r-",label = "ES")
    plt.semilogy(angulos_x_DGH,k_function_DGH(angulos_x_DGH),"g-", label = "DGH")
    plt.semilogy(angulos_x,k_function(angulos_x),"k-", label = "adding the three")
    #plt.plot(angulos_ES,k_ES,"b.")
    #plt.plot(angulos_DGH,k_DGH,"b.")
    plt.grid()
    
    plt.ylabel("$E_{\\nu}^2 \cdot \phi$ $\\left[GeV \cdot cm^{-2} \\cdot s^{-1}\\right]$")
    
    plt.ylim(1,10**3)
    plt.xlim(-90,90)
    
    plt.xticks(np.arange(-90,105,15))
    #ax.spines['bottom'].set_color('white')
    ax.tick_params(axis='x', colors='white')
    
    for i in range(0,len(extremos_CL_dec)-1,2):
        if i == 0:
            ax.axvspan(extremos_CL_dec[i], extremos_CL_dec[i+1], alpha=0.2, color='red', label = "90$\\%$ CL region")
        else:
             ax.axvspan(extremos_CL_dec[i], extremos_CL_dec[i+1], alpha=0.2, color='red')

    plt.legend()
    
    
    ax_sub = fig.add_axes([.15,0.1,.8,.2])
    
    plt.xticks(np.arange(-90,105,15))
    #plt.yticks([0.005,0.010,0.015])
    #plt.ylim(0,0.017)
    plt.xlim(-90,90)
    
    plt.yticks([])
    
    plt.plot(declinations, np.array(probabilities), "b-")
    plt.grid()
    
    plt.ylabel("probability")
    plt.xlabel("$\\delta$ $\\left[{}^{o}\\right]$")
    
    for i in range(0,len(extremos_CL_dec)-1,2):
        ax_sub.axvspan(extremos_CL_dec[i], extremos_CL_dec[i+1], alpha=0.2, color='red')
        
    plt.tight_layout()
    plt.savefig("Plots_Flux_limit/"+GWname+".png")
    
    
