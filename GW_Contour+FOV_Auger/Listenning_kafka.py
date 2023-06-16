# -*- coding: utf-8 -*-
"""
Created on Fri Dec 16 11:25:47 2022

This script once launched, it waits for CGN alerts from LVC collaboration (and IceCube), 
it runs subprocesses that perform plots and computations for each event, and finally
it sends a Circular automatically
"""
import numpy as np
import subprocess
from gcn_kafka import Consumer
import pycurl
import os 
import signal
import SendTheMail as SDM

# If we want to include LVC tests in the annalysis TEST = 'yes', otherwise
# TEST = 'no'

TEST = 'yes'

# If we want to send the email after some dt with the last FITS, then Delay = True
# if we want to send an email for each alert, then Delay= False

Delay = True

dt = str(60*10) # 4 hours for instance


currentdir = os.getcwd()
multi_SUFFIX = ".multiorder.fits"
Flat_resolution_SUFFIX = ".fits.gz"


# Diretory where we will store the FITS files
FITS_destination = currentdir + "/Data_FITS/"

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

# In order to see avaiable alerts sing in in https://gcn.nasa.gov/ , and then
# in client credentials you can add or modify.

consumer.subscribe(['gcn.classic.text.ICECUBE_ASTROTRACK_BRONZE',
                    'gcn.classic.text.ICECUBE_ASTROTRACK_GOLD',
                    'gcn.classic.text.ICECUBE_CASCADE',
                    'gcn.classic.text.LVC_COUNTERPART',
                    'gcn.classic.text.LVC_INITIAL',
                    'gcn.classic.text.LVC_PRELIMINARY',
                    'gcn.classic.text.LVC_RETRACTION',
                    'gcn.classic.text.LVC_UPDATE'])

# 'gcn.classic.text.LVC_EARLY_WARNING',     'gcn.classic.text.LVC_TEST',

# In this Dict all subprocesses are stored with the trigger as a key
processDict = {}


while True:
    for message in consumer.consume():
        
        # Once a alert is received this starts, in value we store the GCN alert
        # as a str. We need to decode to change from binary
        value = message.value()
        value = value.decode("utf-8")
        
        print(value)
    
        lines = value.split('\n')
        
        # In this dic we store all the information from the Alert
        
        dic = {}
        
        # It is customary in Alerts, that the 18th first characters are describing
        # a feature, for instance TRIGGER_NUM:, NOTICE_TYPE:, SKYMAP_FITS_URL:.
        # After that 18 characters comes the information
        
        # strip method deletes border blanks
        
        spline = [lines[0][:17].strip()[:-1],lines[0][17:].strip()]
        
        # The first feature is basically the kind of GCN Alert, GCN/LVC NOTICE,
        # GCN/AMON NOTICE (icecube)... 
        
        # Let's create two auxiliary variables, the first containing one feature and
        # the following one, all the entries related with that feature. For instance
        # the COMMENTS: feature in LVC contains several lines. We store all the lines
        # as different elements in aux_list
        
        aux = spline[0]        
        aux_list = []
        
        for line in lines:
            # In order to skip reading blank lines
            if line.strip():    
                
                spline = [line[:17].strip()[:-1],line[17:].strip()]
                
                #print(spline[1])
                
                # If we continue having the same feature, we include it in the
                # auxiliary list
            
                if aux == spline[0]:
                    aux_list.append(spline[1]) 
                    
                else:
                    
                    # If we have just one element for the feature, we keep that
                    # element just as a str
                    
                    if len(aux_list) == 1:
                        # in order to avoid the "[", "]"
                        dic[aux] = str(aux_list)[2:-2] 
                        
                    # Otherwise we keep the whole list
                    else:
                        dic[aux] = aux_list
                        
                    # Change the auxiliary variables
                    aux = spline[0]
                    aux_list = [spline[1]]
            
            dic[aux] = aux_list # to add the last one
            
        # Value will be a main text for a local email to Wuppertal and 
        # Santiago, but we will include more information
        value += '\n\n'
        
        #print(dic)
        
        #np.save('my_file.npy', dic) 
        
# =============================================================================
# Once we have the dictionary with the alert data we start to run the scripts
# =============================================================================

        if len(dic) > 1: # To avoid first messages (subscribed topic no avaiable:)
        
        
            if 'LVC' in dic['TITLE']:
                
                ######################## LVK EVENTS ##########################
                
                # First we treat the case of no retraction
                
                if 'TEST' in dic['NOTICE_TYPE']:
                    print('This is a test notice')
                    if TEST == 'no':
                        continue
                
                if 'Retraction' not in dic['NOTICE_TYPE']:
                    
                    # we add the number trigger to the path of the file (it is 
                    # just part of the name of the file)
                    
                    FILE_DEST = FITS_destination + dic['TRIGGER_NUM']
                    
                    # This is just in order to check if we have already any 
                    # subproccess running for this concrete trigger, for instance
                    # a preeliminary when we receive an initial alert
                    IDThere = False
                    n = 0        
                    for fi in os.listdir(currentdir + '/Data_FITS'):
                        if dic['TRIGGER_NUM'] in fi:
                                n = n + 1  # number of active subprocesses related to this event
                                IDThere = True
            
                    # We include this information in value
                    value += 'IDThere: ' + str(IDThere) + '\n'    
                    
                    # If it's a new trigger, we just create in processDict a 
                    # new entry with the number of the trigger as a key like in
                    # the dic dictionary
                    
                    Delay_file = 'Data_FITS/file_'+dic['TRIGGER_NUM']+ '.txt'
                    
                    
                    if not IDThere:
                        
                        processDict.update( {dic['TRIGGER_NUM']:[]} )
                        # We include this information in value
                        value += 'No superevent with this ID has plots yet.\n'
                        
                        FILE_DEST = FILE_DEST + '_0' # if it's the first, add a 0 to indicate it
                        
                        # The suffix was missing
                        FILE_DEST += ".fits.gz"
                        
                        if Delay == True :
                            
                            p = subprocess.Popen([ 'python3', 'Delay.py', dt, Delay_file],preexec_fn=os.setpgrp)
                            processDict[dic['TRIGGER_NUM']].append(p)
                            
                            
                    
                    # If it's not the first subprocess related to this trigger,
                    # then just change the FILE_DEST
                    
                    else:
                        
                        value += 'The plots of a superevent with this ID has already been calculated, this is the ' + str(n) + ' actualization\n'
                        
                        FILE_DEST = FILE_DEST + '_' + str(n)
                        
                        # The suffix was missing    
                        FILE_DEST += ".fits.gz"
                        
                    
                    
                    
                    # from now on, we extract the FITS file from the 'SKYMAPS_FITS_URL'
                    # we use the pycurl package for that
                        
                    try:
                        
                        # Now it's customary that alerts contain the .mutiorder.fits
                        # format as the preferred file, so the link includes 
                        # information to this file. We are more used to  the old format  
                        # .fits.gz. Since both formats areincluded, we just change the suffix

                        url = dic['SKYMAP_FITS_URL'] 
                        url = url.replace(multi_SUFFIX, Flat_resolution_SUFFIX)
                        
                        print(url)
                        
                        with open(FILE_DEST, 'wb') as f:
                                
                             c = pycurl.Curl()
                        
                             # we take date from url
                             c.setopt(c.URL, url)
                             # we write it in f
                             c.setopt(c.WRITEDATA, f)
                             c.perform()
                             # close session
                             c.close()
                             
                             value += 'Saved GW info in '+FILE_DEST+'\n'
                             
                    except: print("No hay skymaps") # if we don't have FITS file, then we skip
                    
                    
                    # Here we start all the GCN automation for the event by running 
                    # Contour.py as a subprocess. We also include this information
                    # in the processDic
                    
                    
                    if Delay == False:
                    
                        p = subprocess.Popen([ 'python3', 'Contour.py', "yes", FILE_DEST , dic['TRIGGER_NUM']],preexec_fn=os.setpgrp)
                        processDict[dic['TRIGGER_NUM']].append(p)
                    
                    else:
                        p = [ 'python3', 'Contour.py', "yes", FILE_DEST , dic['TRIGGER_NUM']]
                        
                        with open(Delay_file, 'w') as f:
                            for i in p:
                                f.write(i + ',  ')

                  
                
                # This lines start when we have a RETRACTION event, then we kill 
                # all subprocesses and files related to the trigger
                    
                else:
                    
                    Test = 1
                    # If appears a retraction we kill all subprocesses (and son subprocesses)
                    # In the case of not having any subprocess opened or any file
                    # then it was a Test notice, we don't send any mail, or open any calculation
                    if dic['TRIGGER_NUM'] in processDict:
                        for pi in processDict[dic['TRIGGER_NUM']]:
                            
                            if processDict[dic['TRIGGER_NUM']] == []:
                                Test = 1 
                            else:
                                Test = 0
                            
                            os.killpg(os.getpgid(pi.pid), signal.SIGTERM)
                            pi.kill()
                            pi.communicate()
            
                
                    # with this for, we eliminate all the files and plots
                    # associated with this retracted trigger
                    for dirpath, dirnames, filenames in os.walk(currentdir): 
                         for filename in [f for f in filenames if dic['TRIGGER_NUM'] in f]:
                             
                                 if [f for f in filenames if dic['TRIGGER_NUM'] in f] == []:
                                     Test = 1
                                 else:
                                     Test = 0
                                 
                                 path_with_filename_e = os.path.join(dirpath, filename)
                                 os.remove(path_with_filename_e)
                    
                    if Test == 1:
                        print('This is a test notice')
                        
                        if TEST == 'no':
                            continue
                            
                            
                            
                    value += 'Killed all files and processes for this GW event as it was retracted.'       
                    
            
            
            ##################### ICECUBE EVENTS #############################
        
            elif 'AMON' in dic['TITLE']:
                
                # We basically consider GOLD and Bronze events, cascades have very 
                # low precission measuring the sky arrival direction
                
                if ('ICECUBE Astrotrack Bronze' in dic['NOTICE_TYPE']) or ('ICECUBE Astrotrack Gold' in dic['NOTICE_TYPE']):
                    
                    # Comments[0] basically contains the following text:
                    #     IceCube (Gold/Bronze) event
                    # Then we replace blanks by '_' in order to create a name
                    # for the needed data file that we have to pass to the 
                    # fov_Auger_checker_transient.py
                    
                    dic['COMMENTS'][0] = dic['COMMENTS'][0].replace(' ', '_')
                    
                    # This could be a very long name, but this is something like
                    # an ID for the event.
                    
                    # Now we store the information
                    
                    name   = dic['COMMENTS'][0] + '_'+ dic['STREAM']+ '_' + dic['RUN_NUM'] + '_' + dic['EVENT_NUM']
                    
                    for i in range(len(dic['SRC_RA'])):
                        if (dic['SRC_RA'][i] == '{'):
                            RA_hr  = dic['SRC_RA'][i + 2:i + 4]
                            RA_min = dic['SRC_RA'][i + 6:i + 8]
                            RA_sec = dic['SRC_RA'][i + 10:i + 12]
                            
                    for i in range(len(dic['SRC_DEC'])):
                        if (dic['SRC_DEC'][i] == '{'):
                            DEC_deg  = dic['SRC_DEC'][i + 2:i + 4]
                            DEC_min = dic['SRC_DEC'][i + 6:i + 8]
                            DEC_sec = dic['SRC_DEC'][i + 11:i + 13] 
                            
                    fecha  = dic['DISCOVERY_DATE'][25:33].replace("/","-")
                    fecha = '20'+ fecha[-2:] + '-' + fecha[:-3]
                    time   = fecha + ' ' + dic['DISCOVERY_TIME'][12:22]
                    
                    # We create a data file with all this information 
                    
                    dirfile = currentdir + '/Data_transient/' + name + '.dat'
                    
                    with open(dirfile, 'w') as file:
    
                        file.write('NAME  ' + name + '\n')
                        file.write('RA_hr  ' + RA_hr + '\n')
                        file.write('RA_min ' + RA_min + '\n')
                        file.write('RA_sec ' + RA_sec+ '\n')
                        file.write('DEC_deg   ' + DEC_deg + '\n')
                        file.write('DEC_min   ' + DEC_min+ '\n')
                        file.write('DEC_sec   '+ DEC_deg + '\n')
                        file.write('UTC   ' + time)
                    
                    # We launch the checker transient script.
                        
                    p = subprocess.Popen([ 'python3', 'fov_Auger_checker_transient.py', dirfile])
            
            if Delay == True: continue
                        
            # We send the email with the preeliminary iformation to Santiago and
            # Wuppertal
                
            attachments = []
                    
            SDM.SendTheMail(SDM.sender, SDM.sendermail, SDM.receivers, SDM.receivermails, value, attachments)         
                
            
                
    print(" ")
    print("=============================================================================")
    print(" ")

            
            
            
            


