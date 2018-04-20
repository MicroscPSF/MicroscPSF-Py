#!/usr/bin/env python
"""
Test focal scan G-L PSF.
"""
import numpy

import microscPSF
import microscPSF.microscPSF as msPSF

def test_01():
    mp = msPSF.m_params
    rv = numpy.arange(0.0, 1.01, 0.1)
    zv = numpy.arange(-1.0, 1.01, 0.2)

    fast_rz = msPSF.gLZRFocalScan(mp, rv, zv)
    slow_rz = msPSF.gLZRFocalScanSlow(mp, rv, zv)

    print(numpy.allclose(fast_rz, slow_rz))
    
    print(numpy.max(numpy.abs(fast_rz - slow_rz)))


if (__name__ == "__main__"):
    test_01()
    
