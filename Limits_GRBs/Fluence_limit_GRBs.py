# -*- coding: utf-8 -*-
"""
Created on Thu Feb 16 17:50:40 2023

@author: yagol
"""

import numpy as np
import sqlite3
import pandas
import requests
import astropy.time
import astropy.coordinates
import astropy.units as u
import matplotlib.pyplot as plt

import linecache
import paramiko
import os

from scipy.interpolate import interp1d
from scipy.integrate import quad


Fermi = False # if we want to run just with fermi or with summary

#%%

class GRB:
    
    def __init__(self, data, Fermi = False):
        
        self.data = data
        self.Fermi = Fermi
        
        if Fermi == False:
            self.name = data['GRB_name']
        else:
            self.name = data['GRB_name_Fermi']
        
    def UTCtime(self):
        'Returns the time as a astropy.time.Time object'
        
        dic = self.data
        Fermi = self.Fermi
        
        if Fermi == False:
        
            date = dic['mjd']
            
            time = astropy.time.Time(date, scale = 'utc', format = 'mjd')
            time.format = 'iso'
            
        else:
            
            date = dic['datum'] + ' '
            hours = dic['t_trigger'] // 3600
            minutes = (dic['t_trigger'] % 3600) // 60
            seconds = (dic['t_trigger'] % 3600) % 60
            
            if hours < 10:
                hours = '0' + str(int(hours))
            else:
                hours = str(int(hours))

            
            if minutes < 10:
                minutes = '0' + str(int(minutes))
            else:
                minutes = str(int(minutes))
            
            time = date + hours + ':' + minutes + ':' + str(seconds)
            time = astropy.time.Time(time, scale = 'utc', format = 'iso')
            
        return time
    
    def RA(self):
        
        'returns the RA in degrees'
        
        dic = self.data
        
        RA = dic['ra']
        
        return RA
    
    def Dec(self):
        
        'returns the Declination in degrees'
        
        dic = self.data
        
        Dec = dic['decl']
    
        return Dec   
    
    def T90(self):
        
        dic = self.data
        
        T90 = dic['T90']
        
        return T90
        
    def Alt(self, Observatory):
        
        GRB_ra = self.RA()
        GRB_dec= self.Dec()
        
        GRB_time = self.UTCtime()
        
        Auger = astropy.coordinates.AltAz(obstime= GRB_time, location=Observatory)
        ardec = astropy.coordinates.SkyCoord(ra=GRB_ra*u.degree, dec=GRB_dec*u.degree, frame='icrs')
        altaz = ardec.transform_to(Auger)
        
        return altaz.alt.deg
        
    def Az(self, Observatory):
        
        GRB_ra = self.RA()
        GRB_dec= self.Dec()
        
        GRB_time = self.UTCtime()
        
        Auger = astropy.coordinates.AltAz(obstime= GRB_time, location=Observatory)
        ardec = astropy.coordinates.SkyCoord(ra=GRB_ra*u.degree, dec=GRB_dec*u.degree, frame='icrs')
        altaz = ardec.transform_to(Auger)
        
        return altaz.az.deg

#%%    
    
r = requests.get("https://icecube.wisc.edu/~grbweb_public/GRBweb2.sqlite")
f = open('GRBweb2.sqlite', 'wb').write(r.content)


db = sqlite3.connect('GRBweb2.sqlite')

table_names = pandas.read_sql_query("SELECT * from sqlite_sequence", db)


Fermi_table = pandas.read_sql_query("SELECT * from Fermi_GBM", db)
Summary_table = pandas.read_sql_query("SELECT * from Summary", db)

 
#%%

First_Date = astropy.time.Time('2004-01-01 00:00:00')
#date_limit = astropy.time.Time('2018-09-03 00:00:00')
date_limit = astropy.time.Time('2022-01-03 00:00:00')

if Fermi == True:
    
    dic_Fermi   = {}
    
        
    for i in range(len(Fermi_table)):
        
        fila = Fermi_table.loc[i]
        name = fila['GRB_name_Fermi']
        
        dic_Fermi[name] = GRB(fila, Fermi = True)
        
        # We have to eliminate  GRB without T90, ar or dec, and GRB before 
        # First Date, and Date limit   
        
        if (dic_Fermi[name].UTCtime() <= First_Date) or (dic_Fermi[name].UTCtime() >= date_limit) : 
            
            del dic_Fermi[name]
            continue
            
        if np.isnan(dic_Fermi[name].T90()): 
        
            del dic_Fermi[name]
        
        
else:

    dic_Summary = {}
    
    for i in range(len(Summary_table)):
        
        fila = Summary_table.loc[i]
        name = fila['GRB_name']
        
        dic_Summary[name] = GRB(fila)
        
        if (dic_Summary[name].UTCtime() <= First_Date) or (dic_Summary[name].UTCtime() >= date_limit) : 
            
            del dic_Summary[name]
            continue
            
        if np.isnan(dic_Summary[name].T90()): 
        
            del dic_Summary[name]
    

# =============================================================================
#%% Now we compute the Auger sky for the location of each GRB in order to 
# see if they were into the field of view
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

observatory = astropy.coordinates.EarthLocation(lat=lat_Auger*u.deg, lon=lon_Auger*u.deg, height=alt_Auger*u.m)


GRB_DGH = {}
GRB_DGL = {}
GRB_ES  = {}

T90_DGH = []
T90_DGL = []
T90_ES  = []


discarted = 0

if Fermi == True:

    for key in dic_Fermi:
        
        GRB = dic_Fermi[key]
    
        altitude = GRB.Alt(observatory)
        azimuth  = GRB.Az(observatory)
        
        
        if (altitude > -5) and (altitude <= 0):
            GRB_ES[key]  = GRB
            T90_ES.append(GRB.T90())
        elif (altitude > 0) and (altitude <= 15):
            GRB_DGH[key]  = GRB
            T90_DGH.append(GRB.T90())
        elif (altitude > 15) and (altitude <= 30):
            GRB_DGL[key]  = GRB
            T90_DGL.append(GRB.T90())
        else:
            discarted += 1
            
else: 
    
    for key in dic_Summary:
        
        GRB = dic_Summary[key]
    
        altitude = GRB.Alt(observatory)
        azimuth  = GRB.Az(observatory)
        
        
        if (altitude > -5) and (altitude <= 0):
            GRB_ES[key]  = GRB
            T90_ES.append(GRB.T90())
        elif (altitude > 0) and (altitude <= 15):
            GRB_DGH[key]  = GRB
            T90_DGH.append(GRB.T90())
        elif (altitude > 15) and (altitude <= 30):
            GRB_DGL[key]  = GRB
            T90_DGL.append(GRB.T90())
        else:
            discarted += 1

    
# =============================================================================
#%% T90 Histogram
# =============================================================================

box = {'facecolor': 'none',
       'boxstyle': 'round'
      }


plt.figure(10, figsize = [6,6])

plt.subplot(3,1,1)
hist, bins = np.histogram(T90_DGL, bins=40)
logbins = np.logspace(np.log10(bins[0]),np.log10(bins[-1]),len(bins))
plt.hist(T90_DGL, bins=logbins, color = 'g')
plt.xscale('log')
plt.xticks(color = 'white')
plt.grid()
plt.ylim(0,54)
plt.xlim(10**-2,10**3)
plt.ylabel('Counts DGL')
plt.text(10**(-1),25, str(len(GRB_DGL)) +  ' DGL GRBs', bbox = box )


plt.subplot(3,1,2)
hist, bins = np.histogram(T90_DGH, bins=40)
logbins = np.logspace(np.log10(bins[0]),np.log10(bins[-1]),len(bins))
plt.hist(T90_DGH, bins=logbins, color = 'b')
plt.xscale('log')
plt.xticks(color = 'white')
plt.grid()
plt.ylim(0,54)
plt.xlim(10**-2,10**3)
plt.ylabel('Counts DGH')
plt.text(10**(-1),25, str(len(GRB_DGH)) +  ' DGH GRBs', bbox = box )

plt.subplot(3,1,3)
hist, bins = np.histogram(T90_ES, bins=40)
logbins = np.logspace(np.log10(bins[0]),np.log10(bins[-1]),len(bins))
plt.hist(T90_ES, bins=logbins, color = 'r')
plt.xscale('log')
plt.ylim(0,54) 
plt.xlim(10**-2,10**3)
plt.grid()
plt.ylabel('Counts ES')
plt.xlabel('T90 $\\left[ s \\right]$')
plt.text(10**(-1),25, str(len(GRB_ES)) +  ' ES GRBs', bbox = box )

plt.tight_layout()

plt.savefig('GRB_histogram_Fermi.png')

# =============================================================================
#%% All GRB plot in the FOV of Auger
# =============================================================================

# First we plot the FOV of Auger

x = np.linspace(0,360,1000)

fig, ax = plt.subplots(1, 1, figsize=(9, 7))

ax.fill_between(x, -5, 0, color = 'r', alpha = 0.5, label = 'ES band')
ax.fill_between(x, 0, 15, color = 'b', alpha = 0.5, label = 'DGH band')
ax.fill_between(x, 15, 30, color = 'g', alpha = 0.5,label = 'DGL band')

ax.set_ylim(-90,90)
ax.set_xlim(0,360)

Azimut = []
Altitud= []

n_GRBs = 0


if Fermi == True:
    for key in dic_Fermi:
        
        if dic_Fermi[key].UTCtime() < date_limit:
            
            n_GRBs += 1
        
            Azimut.append(dic_Fermi[key].Az(observatory))
            Altitud.append(dic_Fermi[key].Alt(observatory))

else:
    
    for key in dic_Summary:
        
        if dic_Summary[key].UTCtime() < date_limit:
            
            n_GRBs += 1
        
            Azimut.append(dic_Summary[key].Az(observatory))
            Altitud.append(dic_Summary[key].Alt(observatory))
    
plt.plot(Azimut,Altitud,'x', label = 'GRB')

plt.xlabel('Azhimut $\\left[ degrees \\right]$')
plt.ylabel('Altitude $\\left[ degrees \\right]$')

ax.legend(loc='upper right',bbox_to_anchor=(1.2,1.11),fancybox=True,shadow=True)
plt.tight_layout()
plt.savefig('GRB_sample.png')



# =============================================================================
#%% Fluence limit for each GRB calculation, the number of stations is obtained
# reading the t2 files from remote
# =============================================================================

def Num_Integral(exposure, energy):
    
    f = interp1d(energy,exposure)
    
    def funcion(x):
        return f(x)*x**(-4)
    
    Integration = quad(funcion, min(energy), max(energy),points = energy)
    
    return Integration[0]

def ntanks(Date):
    
    year  = Date.datetime.year ; year = str(year)
    month = Date.datetime.month
    
    if month < 10:
        month = '0' + str(month)
    else:
        month = str(month)
    day   = Date.datetime.day
    
    if day < 10:
        day = '0' + str(day)
    else:
        day = str(day)
    
    hours  = Date.datetime.hour
    minutes = Date.datetime.minute
    seconds = Date.datetime.second + Date.datetime.microsecond*10**(-6)
    
    GWtime = hours*3600 + minutes*60 + seconds
    
    filename = 't2_' + year + '_' + month + '_' + day + '_00h00.dat'

    remotepath = '/data2/auger/data/t2Files/' + year + '/' + month + '/' + filename
    localpath  = os.getcwd() + '/' + filename
    
    sftp.get(remotepath,localpath) # para recibir el archivo del ssh 
    
    time = 0
    tank = 0
    
    aux = 0
    
    if os.path.getsize(filename) == 0:
        os.remove(filename)
        return 0
    
    with open(filename, 'r') as f:    
            lines = f.readlines()   
            i = 0
            while time <= GWtime + 1 + T90:
                
                
                try:
                    j = lines[i]
                    
                    if j[0] == '\x00': #some files are 'damaged'
                        
                        aux = 1
                        break
                    
                    #In each line, we split it and we add the GPS time in the list time
                    line=j.split()
                    time = int(line[1])
                    
                    if i == 0:
                        initial_time = time
                        
                    time = time - initial_time
                    
                    if time > GWtime: break
                    
                    #If it is the first time, the number of tank previous is 0, in another case, is the last entry of ntank
                    if i == 0:
                        aux_tank = 0
                    else:
                        aux_tank = tank
                    
                    #If the first element of the line is '+' we add the length of the line minus two (the + and the time)
                    if line[0]=='+':
                        tank = aux_tank+len(line)-2
                    #If the first element of the line is '-' we substract the length of the line minus two (the - and the time) divided for 4
                    elif line[0]=='-':
                        tank = aux_tank-(len(line)-2)/4
                        
                except IndexError:
                    
                    # In the case that the t2File is not completed because 
                    # all stations stopped
                    
                    aux = 1
                    break
                    
                  
                i = i + 1
    
    #once everything is done, we delete the t2File
    
    os.remove(filename)
    
    if aux : # just when IndexError or damaged file
        return 0
    
    else:
        
        #print(time)
        
        
        return int(tank)
    
    
def list_ntanks(Date, T90, leap):
    
    year  = Date.datetime.year ; year = str(year)
    month = Date.datetime.month
    
    if month < 10:
        month = '0' + str(month)
    else:
        month = str(month)
    day   = Date.datetime.day
    
    if day < 10:
        day = '0' + str(day)
    else:
        day = str(day)
    
    # In GPS time the leap seconds are not considered, so the origin of the day for 
    # the t2 file would be UTC first second of the day - leap seconds
    
    day_GRB = astropy.time.Time(year+'-'+month+'-'+day+' 00:00:00', format = 'iso', scale = 'utc') - leap*u.second
    
    hours  = Date.datetime.hour
    minutes = Date.datetime.minute
    seconds = Date.datetime.second + Date.datetime.microsecond*10**(-6)
    
    GWtime = hours*3600 + minutes*60 + seconds
    
    filename = 't2_' + year + '_' + month + '_' + day + '_00h00.dat'
    
    remotepath = '/data2/auger/data/t2Files/' + year + '/' + month + '/' + filename
    localpath  = os.getcwd() + '/' + filename
    
    sftp.get(remotepath,localpath) # para recibir el archivo del ssh 
    
    if os.path.getsize(filename) == 0:
        os.remove(filename)
        return 0, True
    
    time  = []
    ntank = []
    
    
    with open(filename, 'r') as f:
        for j in f:
            #In each line, we split it and we add the GPS time in the list time
            line=j.split()
            
            if j[0] == '\x00': #some files are 'damaged'
                        f.close()
                        os.remove(filename) 
                        return  0, True 
                    
            time+=[int(line[1])] 
            
            #If it is the first time, the number of tank previous is 0, in another case, is the last entry of ntank
            if len(ntank)==0:
                aux_tank=0
            else:
                aux_tank=ntank[-1]
            
            #If the first element of the line is '+' we add the length of the line minus two (the + and the time)
            if line[0]=='+':
                ntank+=[aux_tank+len(line)-2]
            #If the first element of the line is '-' we substract the length of the line minus two (the - and the time) divided for 4
            elif line[0]=='-':
                ntank+=[aux_tank-(len(line)-2)/4]
                
    # The origin of times should be at the beginning of the day but it might not
    # for some t2-files which are damaged
    
    origin_time = astropy.time.Time(time[0], format = 'gps', scale = 'utc')
    final_time  = astropy.time.Time(time[-1], format = 'gps', scale = 'utc')
    GWtime = Date
    
    time = np.array(time) 
    ntank = np.array(ntank)
    
    # Sometimes we have both highs and lows in terms of the number of stations
    # in that case, the first data is from highs, the second for lows, all the information
    # for that case is contained in the last of the two 
    
    time_2 = [time[0]]
    ntank_2 = [ntank[0]]
    
    for i in range(1,len(time)-1):
        if time[i] == time[i + 1]:
            time_2.append(time[i+1])
            ntank_2.append(ntank[i+1])
        elif time[i] == time[i-1]:
            continue
        else:
            time_2.append(time[i])
            ntank_2.append(ntank[i])
    
    
    # For the last calculation we take the origin of the day to be 0
        
    time   =  np.array(time_2)
    time   =  np.array([timei - day_GRB.gps for timei in time])
    ntank  =  ntank_2
    GWtime =  GWtime.gps - day_GRB.gps
        
    time_list = []
    tank_list = []
    
    aux = origin_time.gps - day_GRB.gps
    end = final_time.gps  - day_GRB.gps
    i = 0
    
    #print(aux, end, GWtime, GWtime + T90)
    
    primera = True
    complete = True # It should be False in the case of having a GRB not completely into the t2
    
    
    # Caso 1: GW sucediendo antes que el data del T2, en este caso no se devuelve nada
    
    if aux > GWtime + T90:
        os.remove(filename)
        return  0, True
    
        
    # Caso 2: GW empieza fuera del T2 fil, pero termina dentro, entonces el complete es Falso
    # pero ponemos límites con lo que podemos aprovechar del GRB
        
    if aux > GWtime:
                complete = False
                
                
    # Caso 5 : GW empieza fuera del T2, despúes del tiempo recogido por el fichero.
                
    if end < GWtime:
            os.remove(filename)
            return  0, True
            
    
    while aux < GWtime + T90:        
        
            try:
                
                # Continuación del caso 2, caso 3 y caso 4
                
                # Casos 2 y 4 aprovechamos todo el tiempo que el GRB esté en el tiempo del t2
                # Caso 3 : el GRB está completamente contenido en el T2
                # Caso 4 : el GRB comienza en el t2 pero acaba fuera de él
                
        
            
                if aux >= GWtime :
                    
                    if primera == True:
                        if aux == GWtime: 
                            primera = False
                            
                        else:
                            if complete == False:
                                # In case 2 the first value will be superior to 
                                # the GWtime, so we should take the first value of 
                                # the times
                                i = i + 1
                            
                            time_list.append(time[i - 1])
                            tank_list.append(ntank[i - 1])
                            
                            primera = False
                    
                    time_list.append(aux)
                    tank_list.append(ntank[i])
                    
                i = i + 1
                aux = time[i]
            
            except IndexError:
                
                        # Caso 4 : En el caso de que el T2 no esté completo y no
                        # contenga el GRB entonces se aprovecha lo que se puede y 
                        # se pone el complete a 0
                        
                        # In the case that the t2File is not completed because 
                        # all stations stopped
                        #os.remove(filename)
                        #return 0
                        # Cuando esto sucede deberia poder aprovecharse cierta cantidad de GRBs
                        complete = False
                        break
        
    
    if len(time_list) == 0:
        time_list = T90
        tank_list = ntank[i - 1]
        
        os.remove(filename)
        
        
        return time_list*tank_list, complete
    
    else:
        
        time_list = time_list[1 :]
        time_list = [GWtime] + time_list + [GWtime + T90]
        time_list = np.array([time_list[i + 1] - time_list[i] for i in range(len(time_list) - 1)])
        tank_list = np.array(tank_list)
        
        
        factor = time_list@tank_list
        
        os.remove(filename)
        
        return factor, complete
        
        
        
    


def av_ntanks(Date):
    
    file = 'UseT2Files_010104_311221.txt'
    First_Date = astropy.time.Time('2004-01-01 00:00:00')
    
    # TimeDelta astropy object
    dt = Date - First_Date
    # absolute difference in days between dates
    dt = dt.value 
    # 3 days periods
    period = int(dt/3)  + 1 # since reading lines with linecache starts at 1.
    
    
    Line = linecache.getline(file, period, module_globals=None).split()
    
    ntank = int(Line[2])
    
    return int(ntank)
        
    
#GRB_DGH = {'GRB190420981' : GRB_DGH['GRB190420981']}

# It will be necessary to use this in order to read the times correctly. 
# At the beginning of the T2-files there is the GPS time without the leap 
# seconds, so we have to eliminate them somehow. Also, before 2004 the number
# of leap seconds is 13

leap_seconds = astropy.time.Time(['2005-12-31 23:59:60','2008-12-31 23:59:60','2016-12-31 23:59:60',
                                  '2012-06-30 23:59:60','2015-06-30 23:59:60'])
leap_seconds = np.array([i.mjd for i in leap_seconds])

First_Date = astropy.time.Time('2004-01-01 00:00:00')


n = 0 # number of GRB before data


ssh = paramiko.SSHClient()
ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
ssh.connect('mastercr1.igfae.usc.es', username="yago.lema", password="r6P#!2!s", banner_timeout=200) #'mastercr1.igfae.usc.es' 'nodo014.inv.usc.es'
sftp = ssh.open_sftp()

#%%

print(' ')
print('=============================================================================')
print(' ')

print('Start DGH')
print(' ')
        

# DGH_File       = 'eff_area_data/eff_area_energy_zenith_DGH_periods_1Jan04-31Aug18.dat'
DGH_File       = 'eff_area_data/eff_area_energy_zenith_DGH_periods_1Jan04-31Dec21.dat'
#lines_DGH_file = 1328784 
lines_DGH_file = 1629360
#change_Energy  = 7144
change_Energy = 8760
#limit = 1786
limit = 2190

DGH_limits_k_tot = []
DGH_limits_k_NC  = []
DGH_limits_k_CC_e  = []
DGH_limits_k_CC_mu = []
DGH_limits_k_CC_tau= []

DGH_T90 = []
DGH_Alt = []

DGH_stack_k_tot    = 0
DGH_stack_k_CC_mu  = 0
DGH_stack_k_CC_e   = 0
DGH_stack_k_CC_tau = 0
DGH_stack_k_NC     = 0


n_no_complete = 0

n_DGH = 0

for key in GRB_DGH:
    
    Date = GRB_DGH[key].UTCtime()
    Alt  = 90 - GRB_DGH[key].Alt(observatory) # This is what is called theta
    T90  = GRB_DGH[key].T90()
    
    actual_leap_seconds = 13 + sum(leap_seconds <= Date.mjd)
    
    
    #print(Date)
    if Date > astropy.time.Time('2017-06-27') and Date < astropy.time.Time('2017-06-29'): continue
    '''
    if Date == astropy.time.Time('2008-11-18 21:00:53.536'): continue
    if Date > astropy.time.Time('2017-06-27') and Date < astropy.time.Time('2017-06-29'): continue
    
    if Date == astropy.time.Time('2009-04-19 23:55:05.051'): continue
    if Date == astropy.time.Time('2009-05-18 01:54:44.517'): continue
    if Date > astropy.time.Time('2009-05-03') and Date < astropy.time.Time('2009-06-17'): continue
    if Date > astropy.time.Time('2009-08-12') and Date < astropy.time.Time('2009-08-19'): continue
    if Date == astropy.time.Time('2011-09-20 08:07:16.410') : continue
    '''
    #ntank    = ntanks(Date) * T90 
    #complete = True
    av_ntank = av_ntanks(Date)
    
    Date = astropy.time.Time(Date)
    ntank, complete = list_ntanks(Date, T90, actual_leap_seconds)  # Contains the Time factor of the exposure, and the number of actual tanks
    
    print(Date, complete)
    
    if complete == False:
        n_no_complete += 1
    '''
    if ntank != tank:
        print(' ')
        print(Date)
        print(' next')
    
    print(ntank, tank, T90)
    '''
    
    
    if ntank == 0 or av_ntank == 0: 
        continue # we skip all events during dead-times
    
    DGH_T90.append(T90)
    DGH_Alt.append(Alt)
    
    if Date > First_Date :
        
        # TimeDelta astropy object
        dt = Date - First_Date
        # absolute difference in days between dates
        dt = dt.value 
        # 3 days periods
        period = int(dt/3) 
        if period + 1 > limit : break #there´s no more exposure since this limit period
        # first 4 lines to be read
        lines_index = period*4 + np.array([1,2,3,4]) 
        aux = lines_index.copy()
        lines_index = list(lines_index)
        
        n_DGH += 1
        
        n = lines_index[0]
        
        while n < lines_DGH_file :
            aux = aux + change_Energy
            n = n + change_Energy
            lines_index += list(aux)
            
        #print(lines_index)
        
        
        data_NC     = {}
        data_CC_e   = {}
        data_CC_mu  = {}
        data_CC_tau = {}
        
        aux_NC     = [[],[]]
        aux_CC_e   = [[],[]]
        aux_CC_mu  = [[],[]]
        aux_CC_tau = [[],[]]
        
        theta = '75.00' # lowest altitude value
        
        altitudes = [75.]
        
        # Now when reading the file, we access just to the lines that we want
        # using linecache.getline, then we create 4 dics of lists, once per 
        # kind of interaction, every key in the dictionaires is the angle, we 
        # have energy and exposure data for each key
        
        for line in lines_index:
            Line = linecache.getline(DGH_File, line, module_globals=None)
            
            if Line != '':
                
                Line = Line.split()
                
                energy   = float(Line[3])
                exposure = float(Line[4])
                
                if Line[2] == theta :

                    if Line[1] == '1':
                        aux_NC[0].append(energy)
                        aux_NC[1].append(exposure*ntank/av_ntank)
                    elif Line[1] == '2':
                        aux_CC_mu[0].append(energy)
                        aux_CC_mu[1].append(exposure*ntank/av_ntank)
                    elif Line[1] == '3':
                        aux_CC_e[0].append(energy)
                        aux_CC_e[1].append(exposure*ntank/av_ntank)
                    elif Line[1] == '4':
                        aux_CC_tau[0].append(energy)
                        aux_CC_tau[1].append(exposure*ntank/av_ntank)
                
                else:

                    data_NC[theta]     = aux_NC
                    data_CC_e[theta]   = aux_CC_e
                    data_CC_mu[theta]  = aux_CC_mu
                    data_CC_tau[theta] = aux_CC_tau
                    
                    theta = Line[2]
                    altitudes.append(float(Line[2]))
                    
                    aux_NC     = [[],[]]
                    aux_CC_e   = [[],[]]
                    aux_CC_mu  = [[],[]]
                    aux_CC_tau = [[],[]]
        
        # We have to save the last data too
                    
        data_NC[theta]     = aux_NC
        data_CC_e[theta]   = aux_CC_e
        data_CC_mu[theta]  = aux_CC_mu
        data_CC_tau[theta] = aux_CC_tau
                    
        # Now that the data is readed and separated, we have to do the computation
        # of the limit for each GRB, we create lists
        
        k_NC     = []
        k_CC_e   = []
        k_CC_mu  = []
        k_CC_tau = []
        k_tot    = []
        
        # let's suppose conservatively that in theta = 90º, the exposure decreases
        # a 50% of the exposure at 89º, we do this in order to consider 
        # GRB between 89 and 90 too, so
        
        data_NC['90.00']     = [np.array(data_NC['89.00'][0]), np.array(data_NC['89.00'][1])*10**(-3)]
        data_CC_e['90.00']   = [np.array(data_CC_e['89.00'][0]), np.array(data_CC_e['89.00'][1])*10**(-3)]
        data_CC_mu['90.00']  = [np.array(data_CC_mu['89.00'][0]), np.array(data_CC_mu['89.00'][1])*10**(-3)]
        data_CC_tau['90.00'] = [np.array(data_CC_tau['89.00'][0]), np.array(data_CC_tau['89.00'][1])*10**(-3)]
                    
        
        for key_NC,key_CC_e,key_CC_mu,key_CC_tau in zip(data_NC,data_CC_e, data_CC_mu, data_CC_tau):
            
            data = data_NC[key_NC]
            energy   = np.array(data[0])
            exposure = np.array(data[1])
            
            Integral_NC = Num_Integral(exposure, energy)
            k_NC.append(2.44/Integral_NC)
            k1 = 2.44/Integral_NC
            
            data = data_CC_e[key_CC_e]
            energy   = np.array(data[0])
            exposure = np.array(data[1])
            
            Integral_CC_e = Num_Integral(exposure, energy)
            k_CC_e.append(2.44/Integral_CC_e)
            k2 = 2.44/Integral_CC_e

            
            data = data_CC_mu[key_CC_mu]
            energy   = np.array(data[0])
            exposure = np.array(data[1])
            
            Integral_CC_mu = Num_Integral(exposure, energy)
            k_CC_mu.append(2.44/Integral_CC_mu)
            k3 = 2.44/Integral_CC_mu

            
            data = data_CC_tau[key_CC_tau]
            energy   = np.array(data[0])
            exposure = np.array(data[1])
        
            Integral_CC_tau = Num_Integral(exposure, energy)
            k_CC_tau.append(2.44/Integral_CC_tau)
            k4 = 2.44/Integral_CC_tau
            
            
            k_tot.append(1/(1/k1+1/k2+1/k3+1/k4)) 
            
        
            
        altitudes.append(90)
            
        f_k_tot    = interp1d(altitudes,k_tot)
        f_k_NC     = interp1d(altitudes,k_NC)
        f_k_CC_e   = interp1d(altitudes,k_CC_e)
        f_k_CC_mu  = interp1d(altitudes,k_CC_mu)
        f_k_CC_tau = interp1d(altitudes,k_CC_tau)
        
        # We are storing the k_grb * T90 factor of the fluence
        
        DGH_limits_k_tot.append(float(f_k_tot(Alt))*T90)
        DGH_limits_k_NC.append(float(f_k_NC(Alt))*T90)
        DGH_limits_k_CC_e.append(float(f_k_CC_e(Alt))*T90)
        DGH_limits_k_CC_mu.append(float(f_k_CC_mu(Alt))*T90)
        DGH_limits_k_CC_tau.append(float(f_k_CC_tau(Alt))*T90)
        
        print(key, GRB_DGH[key].Alt(observatory), float(f_k_tot(Alt))*T90)
        
        
        DGH_stack_k_tot    += 1/(float(f_k_tot(Alt))*T90)
        DGH_stack_k_NC     += 1/(float(f_k_NC(Alt))*T90)
        DGH_stack_k_CC_tau += 1/(float(f_k_CC_tau(Alt))*T90)
        DGH_stack_k_CC_mu  += 1/(float(f_k_CC_mu(Alt))*T90)
        DGH_stack_k_CC_e   += 1/(float(f_k_CC_e(Alt))*T90)
        
        
            
    else:
        n = n + 1
        
        
DGH_stack_k_NC = 1/DGH_stack_k_NC
DGH_stack_k_CC_e = 1/DGH_stack_k_CC_e
DGH_stack_k_CC_mu = 1/DGH_stack_k_CC_mu
DGH_stack_k_CC_tau = 1/DGH_stack_k_CC_tau
DGH_stack_k_tot = 1/DGH_stack_k_tot
        
#%%        

print(' ')
print('=============================================================================')
print(' ')

print('Start ES')
print(' ')
        
# ES_File       = 'eff_area_data/eff_area_energy_zenith_ES_periods_1Jan04-31Aug18.dat'
ES_File       = 'eff_area_data/eff_area_energy_zenith_ES_periods_1Jan04-31Dec21.dat'

#lines_ES_file = 3734526
lines_ES_file = 4579290
#change_Energy  = 1786
change_Energy = 2190

#GRB_ES = {'GRB150721732' : GRB_ES['GRB150721732']}

ES_limits_k = []
ES_T90     = []
ES_Alt     = []
ES_Stack_k = 0

n_ES = 0
        
for key in GRB_ES:
    
    Date = GRB_ES[key].UTCtime()
    Alt  = 90 - GRB_ES[key].Alt(observatory) # This is what is called theta
    T90  = GRB_ES[key].T90()
    
    actual_leap_seconds = 13 + sum(leap_seconds <= Date.mjd)
    
    
    
    #ntank = ntanks(Date) * T90
    ntank, complete = list_ntanks(Date, T90, actual_leap_seconds)  # Contains the Time factor of the exposure, and the number of actual tanks
    
    if complete == False:
        n_no_complete += 1
    
    av_ntank = av_ntanks(Date)
    
    print(Date)
    
    if Date > astropy.time.Time('2017-06-27') and Date < astropy.time.Time('2017-06-29'): continue
    if ntank == 0 or av_ntank == 0: continue # we skip all events during dead-times
    
    ES_T90.append(T90)
    ES_Alt.append(Alt)
    
    if Date > First_Date :
        
        # TimeDelta astropy object
        dt = Date - First_Date
        # absolute difference in days between dates
        dt = dt.value 
        # 3 days periods
        period = int(dt/3) 
        if period + 1 > limit : break #there´s no more exposure since this limit period
        # first line to be read
        lines_index = [period - 1]
        aux = lines_index[0]
        n = lines_index[0]
        
        n_ES += 1
        
        while n < lines_ES_file :
            aux = aux + change_Energy
            n = n + change_Energy
            lines_index.append(aux)
        
        data = {}
        
        aux  = [[],[]]
        
        theta = '90.05' # lowest altitude value
        
        altitudes = [90.05]
        
        # Now when reading the file, we access just to the lines that we want
        # using linecache.getline, then we create 4 dics of lists, once per 
        # kind of interaction, every key in the dictionaires is the angle, we 
        # have energy and exposure data for each key
        
        print(lines_index)
       
        for line in lines_index:
            Line = linecache.getline(ES_File, line, module_globals=None)
            
            
            if Line != '':
                
                Line = Line.split()
                
                energy   = float(Line[2])
                exposure = float(Line[3])
                
                if Line[1] == theta :

                    aux[0].append(energy)
                    aux[1].append(exposure*ntank/av_ntank)
                    
                else:
                    data[theta] = aux
                    
                    theta = Line[1]
                    altitudes.append(float(Line[1]))
                    
                    aux = [[],[]]
        
        # We have to save the last data too
                    
        data[theta] = aux
                    
        # Now that the data is readed and separated, we have to do the computation
        # of the limit for each GRB, we create lists
        
        k = []
        
        # let's suppose conservatively that in theta = 90º, the exposure is 0
        # in order to make the actual computation we will take 10**-6
        data['90.00']  = [np.array(data['90.05'][0]), np.array(data['90.05'][1])*10**(-3)]
        
        for key_inf in data :
            
            inf = data[key_inf]
            energy   = np.array(inf[0])
            exposure = np.array(inf[1])
            
            Integral = Num_Integral(exposure, energy)
            k.append(2.44/Integral)
            
        altitudes.append(90)
        
        f_k = interp1d(altitudes,k)

        ES_limits_k.append(float(f_k(Alt))*T90)
        ES_Stack_k += 1/(float(f_k(Alt))*T90)

    else:
        n = n + 1
        
ES_Stack_k = 1/ES_Stack_k

#%%

print(' ')
print('=============================================================================')
print(' ')

print('Start DGL')
print(' ')
        

n_DGL = 0
        
for key in GRB_DGL:
    
    Date = GRB_DGL[key].UTCtime()
    
    actual_leap_seconds = 13 + sum(leap_seconds <= Date.mjd)
    
    if complete == False:
        n_no_complete += 1
    
    print(Date,complete)
    
    ntank, complete = list_ntanks(Date, T90, actual_leap_seconds)  # Contains the Time factor of the exposure, and the number of actual tanks
    
    if complete == False:
        n_no_complete += 1
        
    av_ntank = av_ntanks(Date)
    
    if Date > astropy.time.Time('2017-06-27') and Date < astropy.time.Time('2017-06-29'): continue
    if ntank == 0 or av_ntank == 0: continue # we skip all events during dead-times
    
    if Date > First_Date :

        n_DGL += 1        
        
    else:
        n = n + 1


#%%     

sftp.close()
ssh.close()




Stack_k = 1/(1/ES_Stack_k + 1/DGH_stack_k_tot)


#%% Plot example to see the Fluence limits result and stacking limit

mini = 10**(-7)
maxi = 10**6

E = np.logspace(17, 20, num=1000, endpoint=True, base=10.0)

def Fluence(E,k):
    return E**(-2)*k

plt.figure(1, figsize = [10,10])

plt.subplot(2,3,1)  
plt.loglog(E, Fluence(E/10**9, DGH_limits_k_CC_e[0]), 'r-', label = 'DGH_CC_e')
plt.loglog(E, Fluence(E/10**9, DGH_stack_k_CC_e), 'r--', label = 'DGH_CC_e_Stack')
for i in range(1,len(DGH_limits_k_CC_e)):

    plt.loglog(E, Fluence(E/10**9, DGH_limits_k_CC_e[i]), 'r-')

plt.yticks(10.**np.arange(np.log10(mini),np.log10(maxi)))
plt.ylim(mini,maxi)
plt.ylabel('Fluence $\\left[GeV\\cdot cm^{-2}\\right]$')

plt.legend()
    
plt.subplot(2,3,2)   
plt.loglog(E, Fluence(E/10**9, DGH_limits_k_CC_mu[0]), 'g-', label = 'DGH_CC_mu')
plt.loglog(E, Fluence(E/10**9, DGH_stack_k_CC_mu), 'g--', label = 'DGH_CC_e_Stack')
for i in range(1,len(DGH_limits_k_CC_mu)):

    plt.loglog(E, Fluence(E/10**9, DGH_limits_k_CC_mu[i]), 'g-')
    
plt.yticks(10.**np.arange(np.log10(mini),np.log10(maxi)))
plt.ylim(mini,maxi)

plt.legend()
    
plt.subplot(2,3,3) 
plt.loglog(E,Fluence(E/10**9,DGH_limits_k_CC_tau[0]), 'c-', label = 'DGH_CC_tau')
plt.loglog(E, Fluence(E/10**9, DGH_stack_k_CC_tau), 'c--', label = 'DGH_CC_tau_Stack')
for i in range(1,len(DGH_limits_k_CC_tau)):
    plt.loglog(E,Fluence(E/10**9,DGH_limits_k_CC_tau[i]),'c-')
    
plt.yticks(10.**np.arange(np.log10(mini),np.log10(maxi)))
plt.ylim(mini,maxi)

plt.legend()
    
plt.subplot(2,3,4) 
plt.loglog(E,Fluence(E/10**9,DGH_limits_k_NC[0]), 'b-', label = 'DGH_NC')
plt.loglog(E, Fluence(E/10**9, DGH_stack_k_NC), 'b--', label = 'DGH_NC_Stack')
for i in range(1,len(DGH_limits_k_NC)):
    plt.loglog(E,Fluence(E/10**9,DGH_limits_k_NC[i]),'b-')
    
plt.yticks(10.**np.arange(np.log10(mini),np.log10(maxi)))
plt.ylim(mini,maxi)
plt.ylabel('Fluence $\\left[GeV\\cdot cm^{-2}\\right]$')

plt.legend()

plt.subplot(2,3,5) 
plt.loglog(E,Fluence(E/10**9,DGH_limits_k_tot[0]), 'k-', label = 'DGH_tot')
plt.loglog(E, Fluence(E/10**9, DGH_stack_k_tot), 'k--', label = 'DGH_tot_Stack')
for i in range(1,len(DGH_limits_k_NC)):
    plt.loglog(E,Fluence(E/10**9,DGH_limits_k_tot[i]),'k-')
    
plt.yticks(10.**np.arange(np.log10(mini),np.log10(maxi)))
plt.ylim(mini,maxi)
plt.xlabel('$\\nu$ energy $\\left[eV\\right]$')

plt.legend()

plt.subplot(2,3,6) 
plt.loglog(E, Fluence(E/10**9,ES_limits_k[0]), 'y-', label = 'ES')
plt.loglog(E, Fluence(E/10**9, ES_Stack_k), 'y--', label = 'ES_Stack')
for i in range(1,len(ES_limits_k)):
    plt.loglog(E,Fluence(E/10**9,ES_limits_k[i]),'y-')
    
plt.yticks(10.**np.arange(np.log10(mini),np.log10(maxi)))
plt.ylim(mini,maxi)

plt.legend()
plt.tight_layout()

plt.savefig('Fluence_limits.png')


# =============================================================================
#%% Stacking limits plot 
# =============================================================================

with open('Waxman-Bahcall.txt','r') as file:
    
    lines = file.readlines()
    
    energy_WB = []
    fluence_WB  = []
    
    for line in lines:
        linea = line.split()
        if line[0] == '#': continue
        
        energy_WB.append(float(linea[0]))
        fluence_WB.append(float(linea[1]))
        

energy_WB = np.array(energy_WB)*10**9 # in eV
fluence_WB  = (n_DGH+n_ES)*np.array(fluence_WB)


    
plt.figure(2,figsize = [7,7])

#plt.loglog(E, Fluence(E/10**9, DGH_stack_k_CC_e), 'r-', label = 'DGH_CC_e', linewidth  = 3)
#plt.loglog(E, Fluence(E/10**9, DGH_stack_k_CC_mu), 'g-', label = 'DGH_CC_mu', linewidth  = 3)
#plt.loglog(E, Fluence(E/10**9, DGH_stack_k_CC_tau), 'c-', label = 'DGH_CC_tau', linewidth  = 3)
#plt.loglog(E, Fluence(E/10**9, DGH_stack_k_NC), 'b-', label = 'DGH_NC', linewidth  = 3)
plt.loglog(E, Fluence(E/10**9, DGH_stack_k_tot), 'b-', label = 'DGH_tot', linewidth  = 3) # m
plt.loglog(E, Fluence(E/10**9, ES_Stack_k), 'r-', label = 'ES', linewidth  = 3) # y
plt.loglog(E, Fluence(E/10**9, Stack_k), 'k--', label = 'Tot', linewidth  = 3)

plt.loglog(energy_WB, fluence_WB, color = 'gray', label = 'Waxman-Bahcall '+str(n_DGH+n_ES)+' GRBs', linewidth  = 3)

box = {'facecolor': 'none',
       'boxstyle': 'round'
      }

plt.text(10**18.5,10**0, str(n_DGH) +  ' DGH GRBs\n ' + str(n_ES) + ' ES GRBs', bbox = box )
    
#plt.yticks(10.**np.arange(np.log10(mini),np.log10(maxi)))



plt.legend(loc = 'lower left')

plt.xlabel('$\\nu$ energy $\\left[eV\\right]$')
plt.ylabel('Fluence $\\left[GeV\\cdot cm^{-2}\\right]$')

plt.yticks(10.**np.arange(-8,2))
plt.ylim(10**(-6.5),10**2)
plt.xlim(10**13,10**19.99)

plt.grid()

plt.tight_layout()
plt.savefig('Stacked_limits.png')

# =============================================================================
#%% Scatter plot k vs T90
# =============================================================================


plt.figure(3, figsize = [7,7])

plt.loglog(DGH_T90, DGH_limits_k_tot, '.b')
plt.loglog(ES_T90, ES_limits_k, 'r.')

plt.ylabel('$k_{GRB}$ $\\left[GeV^{-3}cm^{-2}s^{-1}\\right]$')


# =============================================================================
#%% Scatter plot k vs Thetta 
# =============================================================================
    
plt.figure(4, figsize = [7,7])

plt.semilogy(DGH_Alt, DGH_limits_k_tot, '.b', label = 'DGH')
plt.semilogy(ES_Alt, ES_limits_k, 'r.', label = 'ES')


plt.ylabel('$k_{GRB}$ $\\left[GeV^{-3}cm^{-2}s^{-1}\\right]$')
plt.xlabel('$\\theta_{GRB}(t_{0})$ $\\left[degrees\\right]$')

plt.grid()
plt.legend()
plt.savefig('Scatter_limits.png')

# =============================================================================
#%% Histogram kappa 
# =============================================================================
    
plt.figure(5, figsize = [7,7])

#plt.semilogx([],[])
plt.hist(np.array(np.log10((DGH_limits_k_tot)*np.array(DGH_T90))), bins = 50, label = 'DGH', color = 'blue')
plt.hist(np.array(np.log10((ES_limits_k)*np.array(ES_T90))), bins = 50, label = 'ES', color = 'red')

plt.xlabel('$log10(\\kappa_{GRB})$ $\\left[log10(GeV^{-3}cm^{-2})\\right]$')
plt.ylabel('number of GRBs')

plt.grid()
plt.legend()
plt.savefig('Histogram_kappa.png')





    
