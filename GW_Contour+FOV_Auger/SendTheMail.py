# -*- coding: utf-8 -*-
"""
Created on Sat Nov 26 12:34:14 2022

@author: yagol
"""

import email.utils
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
from email.mime.image import MIMEImage 
from email.mime.base  import MIMEBase
from email import encoders
import os

import smtplib

fitname = 'FITS_files/GW170814_skymap.fits.gz'
GWname  = 'GW170814'

sendermail = 'yago.lema@rai.usc.es'
sender     = 'Yago Lema Capeans'
    
receivermails = ['yagolemacapeans@outlook.es', 'yago.lema@rai.usc.es']#, 'jaime.alvarez@usc.es']
receivers     = ['Yago Lema Capeans', 'Yago Lema Capeans']#, 'Jaime Alvarez Muñiz']



def SendTheMail(sender, sendermail, receivers, receivermails, MT, attachments):
    
    

    msg = MIMEMultipart()
    
    msg['Subject'] = 'Circular follow-up GW Auger, NO CANDIDATES INFO'
    
    msg.attach(MIMEText(MT))
    
    
    msg["From"] = email.utils.formataddr((sender, sendermail))   
    
    for att in attachments:
        
        if ".png" in att:
        
            img = att
            img_data = open(img, "rb").read()
            image = MIMEImage(img_data, name = os.path.basename(img))
            msg.attach(image)
            
        if ".pdf" in att:
        
            pdfname = att
            binary_pdf = open(pdfname, 'rb')
            payload = MIMEBase('application', 'octate-stream', Name=pdfname)
            payload.set_payload((binary_pdf).read())
            encoders.encode_base64(payload)
            payload.add_header('Content-Decomposition', 'attachment', filename=pdfname)
            msg.attach(payload)
            
        if '.txt' in att:
            
            txtname = att
            txt = MIMEText(open(txtname, "r", encoding = 'utf-8').read(), txtname)
            # con esta linea hacemos que el fichero enviado se llame igual que el que leemos
            txt.add_header('Content-Disposition', 'attachment', filename=txtname) 
            msg.attach(txt)
            
    msg["To"] = ', '.join([ email.utils.formataddr((i,j)) for i,j  in zip(receivers, receivermails) ])
    
    
    password = 'r7P#!2!s'
    
    with smtplib.SMTP('smtp-mail.outlook.com',587, '[127.0.0.1]') as server: 
        
            server.connect('smtp-mail.outlook.com',587)
        
            server.ehlo()
            server.starttls()
            server.ehlo()
            server.login(sendermail, password)
            server.send_message(msg)
            print("Successfully sent email") 