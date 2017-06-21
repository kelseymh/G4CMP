#!/usr/bin/python
"""
Run this script after the phononKinematics test binary to plot the resulting
slowness surfaces and group velocities.
"""
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
from matplotlib import rc
from mpl_toolkits.mplot3d import Axes3D

matplotlib.rcParams.update({'font.size': 14})
rc('font',**{'family':'serif','serif':['Computer Modern Roman']})
rc('text', usetex=True)

def do_plot(data, s, c, filename):
    fig = plt.figure()
    ax = Axes3D(fig)
    ax.view_init(30, -20)
    ax.scatter(data[:,0], data[:,1], data[:,2], s=s, c=c)
    ax.set_xlabel('$v_x$ [km/s]')
    ax.set_ylabel('$v_y$ [km/s]')
    ax.set_zlabel('$v_z$ [km/s]')
    #ax.set_xlim([0, 5])
    #ax.set_ylim([0, 5])
    #ax.set_zlim([0, 5])
    #[t.set_rotation(10) for t in ax.get_yticklabels()]
    #[t.set_rotation(-30) for t in ax.get_xticklabels()]
    plt.savefig(filename, format='png')
    plt.clf()

# Load up the data 
# commented lines remove any negative values, so we can just
# plot octant 1
long = np.genfromtxt('phonon_group_vel_long', delimiter=',')
ft = np.genfromtxt('phonon_group_vel_trans_fast', delimiter=',')
st = np.genfromtxt('phonon_group_vel_trans_slow', delimiter=',')

s = 1 # marker size
do_plot(long/1000, s, 'b', 'phonon_group_vel_long.png')
do_plot(ft/1000, s, 'g', 'phonon_group_vel_fast_trans.png')
do_plot(st/1000, s, 'r', 'phonon_group_vel_slow_trans.png')

# Do the same for slowness surfaces
long = np.genfromtxt('phonon_slowness_longi', delimiter=',')
ft = np.genfromtxt('phonon_slowness_trans_fast', delimiter=',')
st = np.genfromtxt('phonon_slowness_trans_slow', delimiter=',')

do_plot(long/1000, s, 'b', 'phonon_slowness_long.png')
do_plot(ft/1000, s, 'g', 'phonon_slowness_fast_trans.png')
do_plot(st/1000, s, 'r', 'phonon_slowness_slow_trans.png')
