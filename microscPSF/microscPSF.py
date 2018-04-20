#!/usr/bin/env python
"""
Generate a PSF using the Gibson and Lanni model.

Note: All distance units are microns.

This is slightly reworked version of the Python code provided by Kyle
Douglass, "Implementing a fast Gibson-Lanni PSF solver in Python".

http://kmdouglass.github.io/posts/implementing-a-fast-gibson-lanni-psf-solver-in-python.html


References:

1. Li et al, "Fast and accurate three-dimensional point spread function computation 
   for fluorescence microscopy", JOSA, 2017.

2. Gibson, S. & Lanni, F. "Experimental test of an analytical model of
   aberration in an oil-immersion objective lens used in three-dimensional
   light microscopy", J. Opt. Soc. Am. A 9, 154–166 (1992), [Originally 
   published in J. Opt. Soc. Am. A 8, 1601-1613 (1991)].

3. Kirshner et al, "3-D PSF fitting for fluorescence microscopy: implementation 
   and localization application", Journal of Microscopy, 2012.

Hazen 04/18
"""
import math
import numpy
import scipy
import scipy.integrate
import scipy.interpolate
import scipy.special


# Internal constants.
num_basis = 100     # Number of rescaled Bessels that approximate the phase function.
rho_samples = 1000  # Number of pupil sample along the radial direction.
               
# Microscope parameters.
m_params = {"M" : 100.0,             # magnification
            "NA" : 1.4,              # numerical aperture
            "ng0" : 1.515,           # coverslip RI design value
            "ng" : 1.515,            # coverslip RI experimental value
            "ni0" : 1.515,           # immersion medium RI design value
            "ni" : 1.515,            # immersion medium RI experimental value
            "ns" : 1.33,             # specimen refractive index (RI)
            "ti0" : 150,             # microns, working distance (immersion medium thickness) design value
            "tg" : 170,              # microns, coverslip thickness experimental value
            "tg0" : 170,             # microns, coverslip thickness design value
            "zd0" : 200.0 * 1.0e+3}  # microscope tube length (in microns).


def configure(mp, wvl):
    # Scaling factors for the Fourier-Bessel series expansion
    min_wavelength = 0.436 # microns
    scaling_factor = mp["NA"] * (3 * numpy.arange(1, num_basis + 1) - 2) * min_wavelength / wvl

    # Not sure this is completely correct for the case where the axial
    # location of the flourophore is 0.0.
    #
    max_rho = min([mp["NA"], mp["ng0"], mp["ng"], mp["ni0"], mp["ni"], mp["ns"]]) / mp["NA"]

    return [scaling_factor, max_rho]


def deltaFocus(mp, zd):
    """
    Return focal offset needed to compensate for the camera being at zd.

    mp - The microscope parameters dictionary.
    zd - Actual camera position in microns.
    """
    a = mp["NA"] * mp["zd0"] / mp["M"]  # Aperture radius at the back focal plane.
    return a*a*(mp["zd0"] - zd)/(2.0*mp["zd0"]*zd)


def gLZRFocalScan(mp, rv, zv, pz = 0.0, wvl = 0.6, zd = None):
    """
    Calculate radial G-L at specified radius and ti values. This is models the PSF
    you would measure by scanning the microscopes focus.

    mp - The microscope parameters dictionary.
    rv - A numpy array containing the radius values.
    zv - A numpy array containing the (relative) z offset values of the coverslip (negative is closer to the objective).
    px - Particle z position above the coverslip (positive values only).
    wvl - Light wavelength in microns.
    zd - Actual camera position in microns. If not specified the microscope tube length is used.
    """
    if zd is None:
        zd = mp["zd0"]

    [scaling_factor, max_rho] = configure(mp, wvl)
    rho = numpy.linspace(0.0, max_rho, rho_samples)

    ti = zv.reshape(-1,1) + mp["ti0"]
    opdt = OPD(mp, rho, ti, pz, wvl, zd)

    # Sample the phase
    #phase = numpy.cos(opdt) + 1j * numpy.sin(opdt)
    phase = numpy.exp(1j * opdt)

    # Define the basis of Bessel functions
    # Shape is (number of basis functions by number of rho samples)
    J = scipy.special.jv(0, scaling_factor.reshape(-1, 1) * rho)

    # Compute the approximation to the sampled pupil phase by finding the least squares
    # solution to the complex coefficients of the Fourier-Bessel expansion.
    # Shape of C is (number of basis functions by number of z samples).
    # Note the matrix transposes to get the dimensions correct.    
    C, residuals, _, _ = numpy.linalg.lstsq(J.T, phase.T)

    b = 2 * numpy.pi * rv.reshape(-1, 1) * mp["NA"] / wvl

    # Convenience functions for J0 and J1 Bessel functions
    J0 = lambda x: scipy.special.jv(0, x)
    J1 = lambda x: scipy.special.jv(1, x)
    
    # See equation 5 in Li, Xue, and Blu
    denom = scaling_factor * scaling_factor - b * b
    R = (scaling_factor * J1(scaling_factor * max_rho) * J0(b * max_rho) * max_rho - b * J0(scaling_factor * max_rho) * J1(b * max_rho) * max_rho)
    R /= denom

    # The transpose places the axial direction along the first dimension of the array, i.e. rows
    # This is only for convenience.
    PSF_rz = (numpy.abs(R.dot(C))**2).T

    # Normalize to the maximum value
    PSF_rz /= numpy.max(PSF_rz)    

    return PSF_rz


def OPD(mp, rho, ti, pz, wvl, zd):
    """
    Calculate phase aberration term.

    mp - The microscope parameters dictionary.
    rho - Rho term.
    ti - Coverslip z offset in microns.
    pz - Particle z position above the coverslip in microns.
    wvl - Light wavelength in microns.
    zd - Actual camera position in microns.
    """
    NA = mp["NA"]
    ns = mp["ns"]
    ng0 = mp["ng0"]
    ng = mp["ng"]
    ni0 = mp["ni0"]
    ni = mp["ni"]
    ti0 = mp["ti0"]
    tg = mp["tg"]
    tg0 = mp["tg0"]
    zd0 = mp["zd0"]
        
    a = NA * zd0 / mp["M"]  # Aperture radius at the back focal plane.
    k = 2.0 * numpy.pi/wvl  # Wave number of emitted light.
    
    OPDs = pz * numpy.sqrt(ns * ns - NA * NA * rho * rho) # OPD in the sample.
    OPDi = ti * numpy.sqrt(ni * ni - NA * NA * rho * rho) - ti0 * numpy.sqrt(ni0 * ni0 - NA * NA * rho * rho) # OPD in the immersion medium.
    OPDg = tg * numpy.sqrt(ng * ng - NA * NA * rho * rho) - tg0 * numpy.sqrt(ng0 * ng0 - NA * NA * rho * rho) # OPD in the coverslip.
    OPDt = a * a * (zd0 - zd) * rho * rho / (2.0 * zd0 * zd) # OPD in camera position.
    
    return k * (OPDs + OPDi + OPDg + OPDt)


