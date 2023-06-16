# -*- coding: utf-8 -*-
"""
Created on Thu Jun 15 10:24:09 2023

@author: yagol
"""

import subprocess
import sys
import os
import time

dt = int(sys.argv[1])
Delay_file = sys.argv[2]

# after dt, the file is runned

time.sleep(dt)

with open(Delay_file, 'r') as f:
    
    linea = f.readlines()[0].split(',  ')

subprocess.run(linea, preexec_fn=os.setpgrp)

os.remove(Delay_file)
