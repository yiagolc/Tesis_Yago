# -*- coding: utf-8 -*-
"""
Created on Fri Dec 16 11:25:47 2022

@author: yagol
"""
import numpy as np
import subprocess
from gcn_kafka import Consumer
import pycurl
import os 
import signal
import SendTheMail as SDM


currentdir = os.getcwd()
multi_SUFFIX = ".multiorder.fit"
Flat_resolution_SUFFIX = ".fits.gz"

# Subscribe to topics and receive alerts, here I am suscribed to all alerts
# by Fermi, IceCube and Ligo-Virgo-Kagra

consumer = Consumer(client_id='40n0r4p4tv7ka9n9ns6dgnmtq',
                    client_secret='1utkg01f1m2gsq44jsgjfd5p56ve1qpu7j71jggfnjijq1j8rfdl')


'''
We can choose if we want to receive text, VOEevent or binary changing the third word

gcn.classic.text     # Just binary text, which is the approach chosen
gcn.classic.voevent  # The lxml.etree
gcn.classic.binary
'''

consumer.subscribe(['gcn.classic.text.FERMI_GBM_ALERT',
                    'gcn.classic.text.FERMI_GBM_FIN_POS',
                    'gcn.classic.text.FERMI_GBM_FLT_POS',
                    'gcn.classic.text.FERMI_GBM_GND_POS',
                    'gcn.classic.text.FERMI_GBM_LC',
                    'gcn.classic.text.FERMI_GBM_POS_TEST',
                    'gcn.classic.text.FERMI_GBM_SUBTHRESH',
                    'gcn.classic.text.FERMI_GBM_TRANS',
                    'gcn.classic.text.FERMI_LAT_GND',
                    'gcn.classic.text.FERMI_LAT_MONITOR',
                    'gcn.classic.text.FERMI_LAT_OFFLINE',
                    'gcn.classic.text.FERMI_LAT_POS_DIAG',
                    'gcn.classic.text.FERMI_LAT_POS_INI',
                    'gcn.classic.text.FERMI_LAT_POS_TEST',
                    'gcn.classic.text.FERMI_LAT_POS_UPD',
                    'gcn.classic.text.FERMI_LAT_TRANS',
                    'gcn.classic.text.FERMI_POINTDIR',
                    'gcn.classic.text.FERMI_SC_SLEW',
                    'gcn.classic.text.ICECUBE_ASTROTRACK_BRONZE',
                    'gcn.classic.text.ICECUBE_ASTROTRACK_GOLD',
                    'gcn.classic.text.ICECUBE_CASCADE',
                    'gcn.classic.text.LVC_COUNTERPART',
                    'gcn.classic.text.LVC_EARLY_WARNING',
                    'gcn.classic.text.LVC_INITIAL',
                    'gcn.classic.text.LVC_PRELIMINARY',
                    'gcn.classic.text.LVC_RETRACTION',
                    'gcn.classic.text.LVC_TEST',
                    'gcn.classic.text.LVC_UPDATE'])

processDict = {}


while True:
    for message in consumer.consume():
        value = message.value()
        value = value.decode("utf-8")
        
        print(value)
    
        lines = value.split('\n')
        
        dic = {}
        
        
        spline = [lines[0][:17].strip()[:-1],lines[0][17:].strip()]
        
        
        aux = spline[0]        
        aux_list = []
        
        for line in lines:
            if line.strip():    #para no leer filas sin nada
                
                # strip method deletes border blanks
                spline = [line[:17].strip()[:-1],line[17:].strip()]
                
                #print(spline[1])
            
                if aux == spline[0]:
                    aux_list.append(spline[1]) 
                    
                else:
                    if len(aux_list) == 1:
                        dic[aux] = str(aux_list)[2:-2] # in order to avoid the "[", "]"
                    else:
                        dic[aux] = aux_list
                    
                    aux = spline[0]
                    aux_list = [spline[1]]
            
            dic[aux] = aux_list # to add the last one
        
        value += '\n\n'
        
        #print(dic)
        
        #np.save('my_file.npy', dic) 
        
# =============================================================================
# Once we have the diccionary with the notice data we start to run the scripts
# =============================================================================

        if len(dic) > 1: # To avoid first messages
        
            if 'LVC' in dic['TITLE']:
                
                if 'Retraction' not in dic['NOTICE_TYPE']:
            
                    FILE_DEST = currentdir + "/Data_FITS/"+dic['TRIGGER_NUM']
                    
                    IDThere = False
                    n = 0        
                    for fi in os.listdir('Data_FITS'):
                        if dic['TRIGGER_NUM'] in fi:
                                n = n + 1
                                IDThere = True
                                
                    value += 'IDThere: ' + str(IDThere) + '\n'    
                    if not IDThere:
                        
                        processDict.update( {dic['TRIGGER_NUM']:[]} )
                        
                        value += 'No superevent with this ID has plots yet.\n'
                        
                        FILE_DEST = FILE_DEST + '_0'
                    
                    else:
                        
                        value += 'The plots of a superevent with this ID has already been calculated, this is the ' + str(n) + ' actualization\n'
                        
                        FILE_DEST = FILE_DEST + '_' + str(n)
                        
                    FILE_DEST += ".fits.gz"
                        
                    try:
                        url = dic['SKYMAP_FITS_URL']
                             
                        url = url.replace(multi_SUFFIX, Flat_resolution_SUFFIX)
                        
                        print(url)
                        
                        with open(FILE_DEST, 'wb') as f:
                                
                             c = pycurl.Curl()
                        
                             # cojemos los datos del url
                             c.setopt(c.URL, url)
                             # los escribimos en f
                             c.setopt(c.WRITEDATA, f)
                             # que lo realice
                             c.perform()
                             # cerramos la sesion de pycurl
                             c.close()
                             
                             value += 'Saved GW info in '+FILE_DEST+'\n'
                             
                    except: print("No hay skymaps") # si no tenemos LVK pasamos
                    
                    
                    p = subprocess.Popen([ 'python3', 'Contour.py', "yes", FILE_DEST , dic['TRIGGER_NUM']],preexec_fn=os.setpgrp)
            
                    processDict[dic['TRIGGER_NUM']].append(p)
                    
                    
                else:
                
                    # If appears a retraction we kill all subprocesses (and son subprocesses) and FITS_files
                    if dic['TRIGGER_NUM'] in processDict:
                        for pi in processDict[dic['TRIGGER_NUM']]:
                            
                            os.killpg(os.getpgid(pi.pid), signal.SIGTERM)
                            pi.kill()
                            pi.communicate()
            
                        
                    for fi in os.listdir('Data_FITS'):
                        if dic['TRIGGER_NUM'] in fi:
                        
                            os.remove(fi) 
                            value += 'Killed fits files and processes for this GW event as it was retracted.'       
                        
            
        
            else: pass
        
            attachments = []
                    
            SDM.SendTheMail(SDM.sender, SDM.sendermail, SDM.receivers, SDM.receivermails, value, attachments)         
                
            
                
    print(" ")
    print("=============================================================================")
    print(" ")

            
            
            
            


