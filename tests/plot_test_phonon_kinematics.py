#!/usr/bin/env python
"""
Run this script after the phononKinematics test binary to plot the resulting
slowness surfaces and group velocities.
"""

# Get command line arguments first, with a usage line

import os, sys
from optparse import OptionParser

parser = OptionParser("Usage: %prog [options] Material (Ge or Si)")
parser.add_option("-v", "--velocity",
                  choices=("group_vel","phase_vel","slowness"),
                  default="all")
parser.add_option("-p", "--phonon",
                  choices=("long","trans_fast","trans_slow"),
                  default="all")

(opt, mat) = parser.parse_args();
if len(args)!=1 or not args in ['Ge','Si']:
    parser.error("Material must be specified as Ge or Si.")

# Configure consistent format for all plots

import numpy as np
import matplotlib
from matplotlib import pyplot as plt
from matplotlib import rc
from mpl_toolkits.mplot3d import Axes3D

matplotlib.rcParams.update({'font.size': 14})
rc('font',**{'family':'serif','serif':['Computer Modern Roman']})
rc('text', usetex=True)

def do_plot(data, label, units, size, color, filename):
    fig = plt.figure()
    ax = Axes3D(fig)
    ax.view_init(30, -20)
    ax.scatter(data[:,0], data[:,1], data[:,2], s=size, c=color)
    ax.set_xlabel(f'${label}_x$ [{units}]')
    ax.set_ylabel(f'${label}_y$ [{units}]')
    ax.set_zlabel(f'${label}_z$ [{units}]')
    #ax.set_xlim([0, 5])
    #ax.set_ylim([0, 5])
    #ax.set_zlim([0, 5])
    #[t.set_rotation(10) for t in ax.get_yticklabels()]
    #[t.set_rotation(-30) for t in ax.get_xticklabels()]
    plt.savefig(filename, format='png')
    plt.clf()

def make_plot(material, velname, modename):
    label = { "group_vel": "vg", "phase_vel": "vp", "slowness" : "s" }
    scale = { "group_vel": 1000, "phase_vel": 1000, "slowness" : 1/1000 }
    units = { "group_vel": "km/s", "phase_vel": "km/s", "slowness" : "s/km" }
    color = { "long": "b", "trans_fast": "g", "trans_slow", "r" }
    size = 1

    fname = f'{material}_{velname}_{modename}'
    vdata = np.genfromtxt(fname, delimiter=',')

    do_plot(data=vdata*scale[velname], label=label[velname],
            units=units[velname], size=size, color=color[modename],
            filename=f'{fname}.png')

# Process all the data

for vel in opt.velocity.choices:
    if opt.velocity == vel or opt.velocity == None:
        for mode in opt.phonon.choices:
            if opt.phonon == mode or opt.phonon == None:
                make_plot(mat, vel, mode)
