# -*- coding: utf-8 -*-
"""
Script that computes the one day flux limit for a puntual source, assuming a normal 
operating status for the SD of the Pierre Auger Observatory. It takes as input the
1 day average exposures of the detector (passed as a file). And saves the limit as
a function of declination
"""

import matplotlib.pyplot as plt
import numpy as np
import os 
import sys
from scipy.interpolate import interp1d
from scipy.integrate import quad

# =============================================================================
# Start reading the exposures
# =============================================================================

cwd = os.path.dirname(os.path.realpath(__file__))

path_ES  = "exposure_average1day_ES.dat"
path_DGH = "exposure_average1day_DGH.dat"
path_DGL = "DeclExposureDGLAvg.dat"


with open(path_ES) as file:
    
    filas_ES  = file.readlines()
    
with open(path_DGH) as file:    
    filas_DGH = file.readlines()
    
with open(path_DGL) as file:    
    filas_DGL = file.readlines()
    
# =============================================================================
# Fluence limit calculation fro ES channel
# =============================================================================

# Introduce a k list for each channel

k_ES  = []
k_DGH = []
k_DGL = []

# Energy and exposure list

exposure = []
energy   = []

# start reading and computing

#%%
    
for fila in filas_ES:
    
    fila = fila.split()
    
    if not (float(fila[1]) == 0 and float(fila[0]) == 0):
        
        # Here we store data just whenever fila[1] != 0 and fila[0] != 0 in
        # order to skip separation rows
        #  0. 0. 0. 
        # We include the second element in order to take into account the 
        # 0 degrees declination
    
        energy.append(float(fila[1]))
        exposure.append(float(fila[2]))
        
        # Current declination angle
        
        angulo_act = float(fila[0])
    
    
    if float(fila[1]) == 0 and float(fila[0]) == 0: 
        
        # Once we have read everything we start the computation
        
        energy = np.array(energy) ; exposure = np.array(exposure)
        
        # Interpolate exposure as a function of energy

        f = interp1d(energy,exposure)
    
        def funcion(x):
            # Integrand of the denominator of the calculation of the limit flux
            return (x**(-2))*f(x)
        
        # Perform the integration using the energies as breaking points
        
        Integration = quad(funcion, min(energy), max(energy),points = energy)
    
        if Integration[0] != 0:
            
            # If we have exposure different from zero we compute the limit.
            # The fluence limit is obtained by multiplying the flux by E^2*T
            # So multiplying by 1 day we obtain de Fluence limit.

            ki = 2.39/Integration[0]*24*3600
            k_ES.append(ki)
            
        else:
            
            # In the case that we don't have exposure, we assing an arbitrarily
            # large value
            
            k_ES.append(10**(5))
        
        # Code to print stuff and plot exposure against declination
        
        #if angulo_act == 55: 
            
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
            
            #x = np.linspace(min(energy),max(energy),100000)
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
        
        # After storing the k computed we initialize again the exposure and energy lists
        
        exposure = []
        energy   = []
        
#%%

for fila in filas_DGH:
    
    fila = fila.split()
    
    if not (float(fila[1]) == 0 and float(fila[0]) == 0):
        
        # In DGH is almost the same, but we have 4 channels, CC_e, CC_mu, CC_tau
        # and NC
        
        energy.append(float(fila[1]))
        exposure.append([float(fila[2]),float(fila[3]),float(fila[4]),float(fila[5])])
        
        angulo_act = float(fila[0])
    
    
    if float(fila[1]) == 0 and float(fila[0]) == 0: 
        
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
            
            # In order to obtain a combined limit we do the following
            
            ki = 1/(1/k1+1/k2+1/k3+1/k4)
            
            k_DGH.append(ki)
        else:
            k_DGH.append(10**(5))
        
        #if angulo_act == 45:
            
           # x = np.linspace(min(energy),max(energy),100000)
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
        
#%%      

energy   = []
exposure = []

for fila in filas_DGL:
    
    line = fila.split()
    
    if line == []:
        continue
    
    if line[0] != 'declination':
        
        energia = line[0].split('^')
        energia = int(energia[0])**float(energia[1])/10**9
        
        exposure_val = float(line[-1])*(10**10)
        
        energy.append(energia)
        exposure.append(exposure_val)
        
    else:
        if energy == []:
            continue
        
        energy = np.array(energy) ; exposure = np.array(exposure)
        
        # Interpolate exposure as a function of energy

        f = interp1d(energy,exposure)
    
        def funcion(x):
            # Integrand of the denominator of the calculation of the limit flux
            return (x**(-2))*f(x)
        
        # Perform the integration using the energies as breaking points
        
        Integration = quad(funcion, min(energy), max(energy),points = energy)
        
        if Integration[0] != 0:
            # If we have exposure different from zero we compute the limit.
            # The fluence limit is obtained by multiplying the flux by E^2*T
            # So multiplying by 1 day we obtain de Fluence limit.

            ki = 2.39/Integration[0]*24*3600
            k_DGL.append(ki)
            
        else:
            
            # In the case that we don't have exposure, we assing an arbitrarily
            # large value
            
            k_DGL.append(10**(5))
            
        
        energy   = []
        exposure = []

# Last one is not included, the file is different for DGL
        
energy = np.array(energy) ; exposure = np.array(exposure)
        
# Interpolate exposure as a function of energy

f = interp1d(energy,exposure)

def funcion(x):
    # Integrand of the denominator of the calculation of the limit flux
    return (x**(-2))*f(x)

# Perform the integration using the energies as breaking points

Integration = quad(funcion, min(energy), max(energy),points = energy)

if Integration[0] != 0:
    # If we have exposure different from zero we compute the limit.
    # The fluence limit is obtained by multiplying the flux by E^2*T
    # So multiplying by 1 day we obtain de Fluence limit.

    ki = 2.39/Integration[0]*24*3600
    k_DGL.append(ki)
    
else:
    
    # In the case that we don't have exposure, we assing an arbitrarily
    # large value
    
    k_DGL.append(10**(5))    

#%%

# Now the declination limits for each band are differents, so when merging 
# different chanels, for instance ES and DGH, we have to take into account 
# this, we construct arrays that contain all

angulos_ES  = np.arange(-55,63)
angulos_DGH = np.arange(-69,56)
angulos_DGL = np.arange(-86,42)
angulos     = np.arange(-86,63)
        
#k_ES = np.array(k_ES) ; k_DGH = np.array(k_DGH)

dif_ES  = 69 - 55 
dif_DGH = 63 - 56


# We include 10**5 values for each declination not included in one specific channel

kES  = []
kDGH = []
kDGL = []
k    = []


i_ES  = 0
i_DGH = 0
i_DGL = 0

for a in angulos:
    
    if a not in angulos_ES:
        kES.append(10**5)
    if a not in angulos_DGH:
        kDGH.append(10**5)
    if a not in angulos_DGL:
        kDGL.append(10**5)
        
    if a in angulos_ES and i_ES == 0:
        kES += k_ES
        i_ES = 1
        
    if a in angulos_DGH and i_DGH == 0:
        kDGH += k_DGH
        i_DGH = 1
    
    if a in angulos_DGL and i_DGL == 0:
        kDGL += k_DGL
        i_DGL = 1
    
    
    
    

'''
kES = np.array([10**5 if i < dif_ES else k_ES[i - dif_ES] for i in range(len(k_ES)+ dif_ES)])
kDGH = np.array([10**5 if i >= len(k_DGH) else k_DGH[i] for i in range(len(k_DGH)+ dif_DGH)])
'''
# This would be the Fluence limit values to be stored

k = 1/(1/np.array(kES)+1/np.array(kDGH)+1/np.array(kDGL))

# After this is just for the plot  

k_function  = interp1d(angulos,k)
k_function_ES  = interp1d(angulos_ES,k_ES)
k_function_DGH = interp1d(angulos_DGH,k_DGH)
k_function_DGL = interp1d(angulos_DGL,k_DGL)


# =============================================================================
#%% Plot
# =============================================================================


angulos_x_ES  = np.linspace(-55,62, 1000)
angulos_x_DGH = np.linspace(-69,55, 1000)
angulos_x_DGL = np.linspace(-86,41, 1000)
angulos_x     = np.linspace(-86,62, 1000)


fig = plt.figure(1)
    
plt.semilogy(angulos_x_ES,k_function_ES(angulos_x_ES),"r-",label = "ES")
plt.semilogy(angulos_x_DGH,k_function_DGH(angulos_x_DGH),"b-", label = "DGH")
plt.semilogy(angulos_x_DGL,k_function_DGL(angulos_x_DGL),"g-", label = "DGH")
plt.semilogy(angulos_x,k_function(angulos_x),"k-", label = "adding the three")
    #plt.plot(angulos_ES,k_ES,"b.")
    #plt.plot(angulos_DGH,k_DGH,"b.")
plt.grid()
    
plt.ylabel(" Fluence $\\left[GeV \cdot cm^{-2} \\right]$")
plt.xlabel("$\\delta$ $\\left[{}^{o}\\right]$")
    
plt.ylim(1,10**3.5)
plt.xlim(-90,90)
    
plt.xticks(np.arange(-90,105,15))

plt.tight_layout()
plt.savefig("Plot_fluences")

# =============================================================================
#%% Saving data 
# =============================================================================

with open("Fluence_ES.txt","w") as file:
    
    for dec,fluence in zip(angulos_ES,k_ES):
        
        file.write(str(dec) + '  ' + str(fluence)+'\n')
        
        
    
with open("Fluence_DGH.txt","w") as file:
        
    for dec,fluence in zip(angulos_DGH,k_DGH):
        
        file.write(str(dec) + '  ' + str(fluence)+'\n')
        
        
with open("Fluence_DGL.txt","w") as file:
        
    for dec,fluence in zip(angulos_DGL,k_DGL):
        
        file.write(str(dec) + '  ' + str(fluence)+'\n')
        
        
with open("Fluence_total.txt","w") as file:
        
    for dec,fluence in zip(angulos, k):
        
        file.write(str(dec) + '  ' + str(fluence)+'\n')


