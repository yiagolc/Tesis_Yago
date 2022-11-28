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


# =============================================================================
# input
# =============================================================================

if(len(sys.argv) < 2) :
    print('You have not include declination')
    sys.exit()

else:
    declination = float(sys.argv[1])


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
  
#%%

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


print(k_function(declination))

#%%
'''
angulos_x_ES  = np.linspace(-55,62, 1000)
angulos_x_DGH = np.linspace(-69,55, 1000)
angulos_x     = np.linspace(-69,62, 1000)

plt.figure(1)
plt.semilogy(angulos_x_ES,k_function_ES(angulos_x_ES),"r-",label = "ES")
plt.semilogy(angulos_x_DGH,k_function_DGH(angulos_x_DGH),"g-", label = "DGH")
plt.semilogy(angulos_x,k_function(angulos_x),"k-", label = "suma")
plt.plot(angulos_ES,k_ES,"b.")
plt.plot(angulos_DGH,k_DGH,"b.")
plt.grid()

plt.xticks(np.arange(-70,70,10))

plt.xlabel("$\\delta$ $\\left[{}^{o}\\right]$")
plt.ylabel("$E_{\\nu}^2 \cdot \phi$ $\\left[GeV \cdot cm^{-2} \\cdot s^{-1}\\right]$")

'''
