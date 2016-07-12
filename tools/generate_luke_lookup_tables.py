# -*- coding: utf-8 -*-
"""
Created on Mon Apr 18 17:10:15 2016

@author: Rob Agnese

Generate look-up table for G4CMP's LukeScatter process
"""

import luke_emission_functions as lef
import numpy as np
import sys
import getopt
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
import warnings

warnings.simplefilter('always', RuntimeWarning)


def main(argv):
    try:
        opts, args = getopt.getopt(argv, "p", ["plots"])
    except getopt.GetoptError:
        sys.exit(2)

    do_plots = False
    for opt, arg in opts:
        if opt in ("-p", "--plots"):
            do_plots = True

    if len(args) != 10:
        print('usage: python [-p | --plots] directory min_field max_field '
              'field_bins min_mach max_mach mach_bins min_theta max_theta '
              'theta_bins')
        sys.exit(2)

    crystal_dir = args[0]

    try:
        conf_file = open(crystal_dir+'/config.txt')
        conf_data = conf_file.readlines()
        for line in conf_data:
            # Remove comment and tokenize:
            split_line = line.partition('#')[0].split()
            if len(split_line) == 0:
                continue
            if split_line[0] == 'vsound':
                vl = float(split_line[1]) / 1000  # from m/s to um/ns
            elif split_line[0] == 'l0_e':
                l0_e = float(split_line[1]) * 1e6  # from m to um
            elif split_line[0] == 'l0_h':
                l0_h = float(split_line[1]) * 1e6  # from m to um
            elif split_line[0] == 'emass':
                m_components = np.array([float(m) for m in split_line[1:]]) *\
                                9.11e-31  # kg
                m_e = len(m_components)/sum(1/m_components)
            elif split_line[0] == 'hmass':
                m_components = np.array([float(m) for m in split_line[1:]]) *\
                                9.11e-31  # kg
                m_h = len(m_components)/sum(1/m_components)
    except:
        print("Couldn't read config.txt from crystal directory")
        sys.exit(1)

    field_list = np.linspace(float(args[1]), float(args[2]), num=int(args[3]))
    mach_list  = np.linspace(float(args[4]), float(args[5]), num=int(args[6]))
    theta_list = np.linspace(float(args[7]), float(args[8]), num=int(args[9]))

    fill_look_up_table(crystal_dir+'/luke_scatter_times_elec.dat', vl, l0_e,
                       m_e, -1, field_list, mach_list, theta_list, do_plots)
    fill_look_up_table(crystal_dir+'/luke_scatter_times_hole.dat', vl, l0_h,
                       m_h, +1, field_list, mach_list, theta_list, do_plots)


def fill_look_up_table(file_path, vl, l0, m, qsign,
                       field_list, mach_list, theta_list, do_plots):
    # Physical Constants:
    alpha = 1e-12  # convert Efield from V/m to terms of um and ns
    hbar = 1.0545718e-31  # in um and ns units
    q = 1.602e-19  # coulombs - NOTE: sign doesn't matter

    # Derived physical parameters:
    kl = vl * m / hbar
    k0_list = mach_list * kl
    dkdt_list = q * alpha * field_list / hbar
    cos_theta_list = qsign * np.cos(theta_list)

    outfile = create_file(file_path, field_list, theta_list, mach_list)

    # p0 is the initial guess for the fitting coefficients (A, mu and sigma)
    p0 = [1., 1., 1.]

    # initial guess for time window of pdf
    t_test_list = np.linspace(0, 50, 500)

    for i, dkdt in enumerate(dkdt_list):
        for j, cos_theta in enumerate(cos_theta_list):
            for k, k0 in enumerate(k0_list):
                params = [l0, vl, kl, k0, dkdt, cos_theta]
                t_list = adjust_time_bounds(t_test_list, *params)
                pdf_values = [lef.pdf(t, *params) for t in t_list]
                if dkdt == 0 and k0 < kl:
                    coeffs, _ = curve_fit(flat, t_list, pdf_values, p0=p0)
                    fit = flat(t_list, *coeffs)
                    fit_type = 'flat'
                else:
                    try:
                        coeffs, _ = curve_fit(gauss, t_list, pdf_values, p0=p0)
                        fit = gauss(t_list, *coeffs)
                        fit_type = 'gaus'
                    except RuntimeError:
                        coeffs, _ = curve_fit(expon, t_list, pdf_values, p0=p0)
                        fit = expon(t_list, *coeffs)
                        fit_type = 'exp'

                outfile.write('%s %f %f\n' % (fit_type, coeffs[1], coeffs[2]))

                if do_plots:
                    plt.plot(t_list, pdf_values, '.b', label='Data')
                    plt.plot(t_list, fit, '.r', label='Fit: ' + fit_type)
                    plt.legend()
                    plt.xlabel('Time [ns]')
                    plt.title('PDF of Luke Phonon Scatter Times')
                    if qsign < 0:
                        carr = 'elec'
                    else:
                        carr = 'hole'
                    plt.savefig("pdf_" + carr +
                                "_E-" + str(round(field_list[i])) +
                                "_theta-" + str(round(theta_list[j])) +
                                "_mach-" + str(round(mach_list[k])) +
                                ".png")
                    plt.clf()


def create_file(file_path, fields, thetas, machs):
    outfile = open(file_path, 'w')

    outfile.write("number of fields, first field [V/m], last field [V/m], " +
                  "number of thetas, first theta [rad], last theta [rad], " +
                  "number of machs, first mach, last mach\n")
    outfile.write(' '.join([str(len(fields)), str(fields[0]), str(fields[-1]),
                            str(len(thetas)), str(thetas[0]), str(thetas[-1]),
                            str(len(machs)), str(machs[0]), str(machs[-1])]))
    outfile.write("\nfit type [0 = Gaus., 1 = Expon., 2 = Zero], mean [ns], std [ns]\n")

    return outfile


# Define model functions to be used to fit to the data we generate:
def gauss(x, *p):
    A, mu, sigma = p
    return A*np.exp(-(x-mu)**2/(2.*sigma**2))


def expon(x, *p):
    A, mu, _ = p
    return A*np.exp(-x/mu)


def flat(x, *p):
    A, _, _ = p
    return A*np.ones_like(x)


def adjust_time_bounds(tvals, *params):
    """Given a time range, expand or shrink it to fit the PDF data"""
    thresh = 1e-6
    yvals = [lef.pdf(t, *params) for t in tvals]
    nonzero_idx = [i for i, v in enumerate(yvals) if v > thresh]
    try:
        max_idx = nonzero_idx[-1]
    except IndexError:
        max_idx = len(tvals) - 1
    tmax = tvals[max_idx] + 1
    if max_idx == len(tvals) - 1:
        while lef.pdf(tmax, *params) > thresh:
            tmax = tmax + 1

    return np.linspace(0, tmax, 500)


if __name__ == "__main__":
    main(sys.argv[1:])
