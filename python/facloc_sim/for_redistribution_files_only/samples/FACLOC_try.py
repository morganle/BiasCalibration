#!/usr/bin/env python
"""
Sample script that uses the facloc_sim module created using
MATLAB Compiler SDK.

Refer to the MATLAB Compiler SDK documentation for more information.
"""

from __future__ import print_function
import facloc_sim
import matlab

my_facloc_sim = facloc_sim.initialize()

thetaIn = matlab.double([0.3333333333333333, 1.26, 0.83, 1.96, 0.83], size=(1, 5))
xIn = matlab.double([0.84, 0.84, 0.14, 0.68, 0.35, 0.5, 0.54, 0.2, 0.68, 0.36, 0.93, 0.05, 0.23, 0.93], size=(1, 14))
numTrucksIn = matlab.double([5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0], size=(1, 7))
runlengthIn = matlab.double([1000.0], size=(1, 1))
seedIn = matlab.double([1.0], size=(1, 1))
ontimeOut, DaysOut, TraceOut = my_facloc_sim.FACLOC(thetaIn, xIn, numTrucksIn, runlengthIn, seedIn, nargout=3)
print(ontimeOut, DaysOut, TraceOut, sep='\n')

my_facloc_sim.terminate()
