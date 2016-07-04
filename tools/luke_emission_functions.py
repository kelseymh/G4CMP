# -*- coding: utf-8 -*-
"""
Created on Mon Mar 21 12:21:48 2016

@author: Rob Agnese

Calculate Luke phonon emission rate, pdf of emission times, and mean free path
"""

from scipy.integrate import quad
import numpy as np
import warnings


def kmag(t, k0, dkdt, cos_theta):
    arg = k0**2 + 2.0 * k0 * dkdt * cos_theta * t + (dkdt * t)**2
    if arg < 0:
        warnings.warn("Argument to numpy.sqrt(arg) is negative. \n"
                      "Returning zero instead of NaN. \n"
                      "    t = " + str(t) + "\n"
                      "    k0 = " + str(k0) + "\n"
                      "    dk/dt = " + str(dkdt) + "\n"
                      "    cos(theta) = " + str(cos_theta) + "\n"
                      "    arg = " + str(arg),
                      category=RuntimeWarning)
        return 0
    else:
        return np.sqrt(arg)


def xmag(t, v0, dvdt, cos_theta):
    arg = v0**2 + v0 * dvdt * cos_theta * t + 0.25 * (dvdt * t)**2
    if arg < 0:
        warnings.warn("Argument to numpy.sqrt(arg) is negative. \n"
                      "Returning zero instead of NaN. \n"
                      "    t = " + str(t) + "\n"
                      "    v0 = " + str(v0) + "\n"
                      "    dv/dt = " + str(dvdt) + "\n"
                      "    cos(theta) = " + str(cos_theta) + "\n"
                      "    arg = " + str(arg),
                      category=RuntimeWarning)
        return 0
    else:
        return t * np.sqrt(arg)


def rate(t, l0, vl, kl, k0, dkdt, cos_theta):
    k1 = kmag(t, k0, dkdt, cos_theta)
    if k1 <= kl:
        dndt = 0.
    else:
        dndt = vl/3.0/l0 * (k1/kl)**2.0 * (1.0 - kl/k1)**3.0
    return dndt


def pdf(t, l0, vl, kl, k0, dkdt, cos_theta):
    return rate(t, l0, vl, kl, k0, dkdt, cos_theta) * \
        np.exp(-quad(rate, 0., t, args=(l0, vl, kl, k0, dkdt, cos_theta))[0])


def mfp_integrand(t, l0, vl, v0, dvdt, kl, k0, dkdt, cos_theta):
    return xmag(t, v0, dvdt, cos_theta) * pdf(t, l0, vl, kl, k0, dkdt, cos_theta)


def mfp(mach, Efield, theta, alpha=1e-12, hbar=1.0545718e-31, q=1.602e-19,
        l0=108, vl=5.4, m=.35*9.11e-31, full_integration_info=0):
    kl = vl * m / hbar
    v0 = mach * vl
    dvdt = q * alpha * Efield / m
    k0 = mach * kl
    dkdt = q * alpha * Efield / hbar
    cos_theta = np.cos(theta)
    denom = quad(pdf, 0., np.inf,
                 args=(l0, vl, kl, k0, dkdt, cos_theta))[0]
    return quad(mfp_integrand, 0., np.inf,
                args=(l0, vl, v0, dvdt, kl, k0, dkdt, cos_theta))[0]/denom
