#!/usr/bin/env python
"""
Test focal scan G-L PSF.
"""
import numpy

import microscPSF
import microscPSF.microscPSF as msPSF

def test_01():
    """
    Particle on surface.
    """
    mp = msPSF.m_params
    rv = numpy.arange(0.0, 1.01, 0.1)
    zv = numpy.arange(-1.0, 1.01, 0.2)

    fast_rz = msPSF.gLZRFocalScan(mp, rv, zv)
    slow_rz = msPSF.gLZRFocalScanSlow(mp, rv, zv)

    assert (numpy.allclose(fast_rz, slow_rz))


def test_02():
    """
    Particle above surface.
    """
    mp = msPSF.m_params
    rv = numpy.arange(0.0, 1.01, 0.1)
    zv = numpy.arange(-1.0, 1.01, 0.2)

    fast_rz = msPSF.gLZRFocalScan(mp, rv, zv, pz = 0.5)
    slow_rz = msPSF.gLZRFocalScanSlow(mp, rv, zv, pz = 0.5)

    assert (numpy.allclose(fast_rz, slow_rz, atol = 1.0e-4, rtol = 1.0e-4))


def test_03():
    """
    Detector offset.
    """
    mp = msPSF.m_params
    rv = numpy.arange(0.0, 1.01, 0.1)
    zv = numpy.arange(-1.0, 1.01, 0.2)

    zd = mp["zd0"] + 1000
    fast_rz = msPSF.gLZRFocalScan(mp, rv, zv, zd = zd)
    slow_rz = msPSF.gLZRFocalScanSlow(mp, rv, zv, zd = zd)

    assert (numpy.allclose(fast_rz, slow_rz))


def test_04():
    """
    Particle scan.
    """
    mp = msPSF.m_params
    rv = numpy.arange(0.0, 1.01, 0.1)
    pv = numpy.arange(0.0, 2.01, 0.1)

    fast_rz = msPSF.gLZRParticleScan(mp, rv, pv)
    slow_rz = msPSF.gLZRParticleScanSlow(mp, rv, pv)

    assert (numpy.allclose(fast_rz, slow_rz, rtol = 1.0e-4, atol = 1.0e-4))


def test_05():
    """
    Particle scan, focus offset.
    """
    mp = msPSF.m_params
    rv = numpy.arange(0.0, 1.01, 0.1)
    pv = numpy.arange(1.0, 3.01, 0.2)

    fast_rz = msPSF.gLZRParticleScan(mp, rv, pv, zv = -2.0)
    slow_rz = msPSF.gLZRParticleScanSlow(mp, rv, pv, zv = -2.0)

    assert (numpy.allclose(fast_rz, slow_rz, rtol = 1.0e-3, atol = 1.0e-3))


if (__name__ == "__main__"):
    test_01()
    test_02()
    test_03()
    test_04()
    test_05()
