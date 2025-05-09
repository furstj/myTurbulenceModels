#!/usr/bin/env python

from numpy import *
U = 20
nu = 1.5e-5

ws = loadtxt("postProcessing/surfaces/10000/wallShearStress_wall.raw")
ws = ws[ws[:,0].argsort()]
cf = -ws[:,3]*2/U**2
Rex = ws[:,0]*U/nu
for k in range(1,len(cf)):
    if cf[k] > cf[k-1]:
        print("Re_xt = ", Rex[k])
        break
