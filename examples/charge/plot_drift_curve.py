#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 15 09:37:37 2015

@author: rob
"""
import os
import glob
import csv
from matplotlib import pyplot as plt

try:
    suffix = os.environ['G4CMP_HIT_SUFFIX']
except KeyError:
        print("Need to set G4CMP_HIT_SUFFIX env variable to match "
              "the drift_curve macro")
        exit(1)

# iZip is 2.54 cm thick and each track begins in the middle of the zip:
dx = 1.27  # cm

files = glob.glob(''.join(('epos_*', suffix, '.txt')))

# File names should all look like: "epos_###v-" + suffix + ".txt",
# where ### can represent any number of digits, with or without a decimal
field = [float(e)/(2.0*dx)  # convert to volt/cm
         for e in [txt[5:txt.find(''.join(('v-', suffix, '.txt')))]
         for txt in files]]

dt_e = []
dt_h = []
for file in files:
    temp_h = 0.0
    temp_e = 0.0
    count_e = 0
    count_h = 0
    with open(file) as text:
        reader = csv.DictReader(text)
        for line in reader:
            if line["Particle Name"] == "G4CMPDriftHole":
                count_h += 1
                temp_h += (float(line["Final Time [ns]"]) - 
                          float(line["Start Time [ns]"])) * 1e-9
            elif line["Particle Name"] == "G4CMPDriftElectron":
                count_e += 1
                temp_e += (float(line["Final Time [ns]"]) - 
                          float(line["Start Time [ns]"])) * 1e-9
    dt_h.append(temp_h/float(count_h))
    dt_e.append(temp_e/float(count_e))

speed_h = [dx/t/1e5 for t in dt_h]  # convert to km/s
speed_e = [dx/t/1e5 for t in dt_e]  # convert to km/s

plt.plot(field, speed_e, 'rs', label='G4CMP Electrons')
plt.plot(field, speed_h, 'ro', label='G4CMP Holes')
plt.legend(loc='upper left', fancybox=True, shadow=True)
plt.xlabel("Electric Field [V/cm]")
plt.ylabel("Drift Speed [km/s]")
#plt.ylim([5, 45])
plt.title("Charge Drift Speeds in Germanium")
plt.savefig(''.join(('Charge_Drift_Speed_', suffix, '.png')))

exit(0)
