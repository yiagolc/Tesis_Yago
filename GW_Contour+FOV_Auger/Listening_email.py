# -*- coding: utf-8 -*-
"""
This program has the aim of listening GCN events and in the case of a GW event
it reads the xml file, extracts the FITS_file and starts Contour.py and 
CL_coverage.py. when finished the Script send two emails, the first containing
information of the event extracted from the xml file and the second email 
(the circular) contain some results extracted from the calculation
including the two plots.

In order to work well the Scripts CL_coverage.py, Contour.py  and fov_Auger_GW_builder.py
should have the 'directory' part commented and the Input part uncommented, also
Contour.py should have run_fov_auger = 'yes' at the begining of the script
"""


# =============================================================================
# modules
# =============================================================================

import gcn
import gcn.handlers
import gcn.notice_types

import pycurl

from string import Template

import subprocess

import healpy as hp

import astropy.time

import smtplib

import email.utils
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
from email.mime.image import MIMEImage
import os
import sys



sendermail = 'yago.lema@rai.usc.es'
sender     = 'Yago Lema Capeans'
    
receivermails = ['yagolemacapeans@outlook.es', 'yago.lema@rai.usc.es']#, 'jaime.alvarez@usc.es']
receivers     = ['Yago Lema Capeans', 'Yago Lema Capeans']#, 'Jaime Alvarez Muñiz']
    


fitname = 'FITS_files/GW170814_skymap.fits.gz'
GWname  = 'GW170814'


#%%

# GCN Servers                 PUBLIC 1.1/2.0    PRIVATE 1.1/2.0
  # 209.208.78.170 (Atlantic_2)  8099             8092  8096
  # 45.58.43.186   (Atlantic_3)       8099
  # 50.116.49.68   (Linode)      8099                   8096
  # 68.169.57.253  (eApps)            8099        8092  8096


# =============================================================================
# Send email function
# =============================================================================

def SendTheMail(sender, sendermail, receivers, receivermails, MT, attachments):
    
    

    msg = MIMEMultipart()
    
    msg['Subject'] = 'Test mail con resultados'
    
    msg.attach(MIMEText(MT))
    
    
    msg["From"] = email.utils.formataddr((sender, sendermail))   
    
    for att in attachments:
        
        if ".png" in att:
        
            img = att
            img_data = open(img, "rb").read()
            image = MIMEImage(img_data, name = os.path.basename(img))
            msg.attach(image)
            
        if '.txt' in att:
            
            txtname = att
            txt = MIMEText(open(txtname, "r", encoding = 'utf-8').read(), txtname)
            # con esta linea hacemos que el fichero enviado se llame igual que el que leemos
            txt.add_header('Content-Disposition', 'attachment', filename=txtname) 
            msg.attach(txt)


    
    msg["To"] = ', '.join([ email.utils.formataddr((i,j)) for i,j  in zip(receivers, receivermails) ])
    
    
    password = 'r6P#!2!s'
    
    with smtplib.SMTP(host='smtp-mail.outlook.com', port=587) as server:
        
            server.ehlo()
            server.starttls()
            server.ehlo()
            server.login(sendermail, password)
            server.send_message(msg)
            print("Successfully sent email")  
            
    

# =============================================================================
# This is the handler for the GCN Listening
# =============================================================================


@gcn.handlers.include_notice_types(
  gcn.notice_types.LVC_EARLY_WARNING,
  gcn.notice_types.LVC_PRELIMINARY,
  gcn.notice_types.LVC_INITIAL,
  gcn.notice_types.LVC_UPDATE,
  gcn.notice_types.LVC_RETRACTION)
def process_gcn(payload, root):

    # VERY IMPORTANT! Use the following code
    # to respond to only real 'observation' events:
    if root.attrib['role'] != 'observation':
        return

    params = { elem.attrib['name']: elem.attrib['value'] for elem in root.iterfind('.//Param') }

    MT = 'We\'re having a '+str(params['AlertType'])+' '+str(root.attrib['role'])+' event notice!\n'
    
    for key, value in params.items(): # bucle a keys y values del diccionario como un zip
        MT += str(key)+' = ' + str(value)+'\n'
        MT += '\n'
    
    if params['AlertType'] != 'Retraction' and 'skymap_fits' in params:
        
        FILE_DEST = 'FITS_files/'+params['GraceID']+'-'+str(params['Pkt_Ser_Num'])+'-'+str(params['AlertType'])+'.fits.gz'

        #diff
        IDThere = False
        
        n = 0
        for fi in os.listdir('FITS_files'):
            if params['GraceID'] in fi:
                n = n + 1
                IDThere = True
        MT += 'IDThere: ' + str(IDThere) + '\n'
        if not IDThere:
            
            MT += 'No superevent with this ID has plots yet.\n'
            
            FILE_DEST = FILE_DEST + '_0'
        
        else:
            
            MT += 'The plots of a superevent with this ID has already been calculated, this is the ' , n, ' actualization\n'
            
            FILE_DEST = FILE_DEST + '_' + str(n)
            
        with open(FILE_DEST, 'wb') as f:
             
             c = pycurl.Curl()
        
             # cojemos los datos del url
             c.setopt(c.URL, params['skymap_fits'])
             # los escribimos en f
             c.setopt(c.WRITEDATA, f)
             # que lo realice
             c.perform()
             # cerramos la sesion de pycurl
             c.close()
             
             MT += 'Saved GW info in '+FILE_DEST+'\n'
                        
        # Now we run the main scripts for this new fits file
        
        coverage = subprocess.run([ 'python3', 'CL_coverage.py', f ], capture_output=True, text=True)
        
        a = coverage.stdout.split('\n')
        
        
        GCN_ID = params['GraceID']
        
        GWtime = a[0]
        GWname = a[1]
        
        fov_prob_first = round(float(a[2]),2) 
        t_in  = round(float(a[3])) 
        t_fin = round(float(a[4]))
        
        subprocess.run([ 'python3', 'Contour.py', f ])
        
        
        with open('MT_Circular.txt', 'r', encoding='utf-8') as template_file:
            template_file_content = template_file.read()
            
            template_file = Template(template_file_content)
            
            MT_Circular = template_file.substitute(GWtime = GWtime, GWname = GWname, fov_prob_first = fov_prob_first, name_trigger = GWname, GCN_ID = GCN_ID, t_in = t_in, t_fin = t_fin )
            
        
        attachments = ["CL_Coverage/GW_confidence_region_coverage" + GWname + ".png", 
                   "Fov_Contours/fov_Auger_{0}_mollweide.png".format(GWname) ]
        
        
        SendTheMail(sender, sendermail, receivers, receivermails, MT_Circular, attachments)
            
             
    elif params['AlertType'] == 'Retraction':
        for fi in os.listdir('FITS_files'):
            if params['GraceID'] in fi:
            
                os.remove(fi) 
                MT += 'Killed all follower process for this GW event as it was retracted.'       
                
                attachments = []
             
    else:
        MT += 'skymap_fits NOT in params! Look at the GCN notice in particular!'
        
        attachments = [] 
      
    
    SendTheMail(sender, sendermail, receivers, receivermails, MT, attachments)
        
# =============================================================================
# finally we send the email
# =============================================================================
    
                   
           
            
PORT = '8099'
HOST = '68.169.57.253'

if len(sys.argv) > 1 and '.' in sys.argv[1]:
    HOST = sys.argv[1]
    if len(sys.argv) > 2:
      try:
        int(sys.argv[2])
        PORT = sys.argv[2]
      except Exception:
        pass

print('listening to host',HOST,'port',PORT)
gcn.listen(host=HOST, port=int(PORT), handler=process_gcn)