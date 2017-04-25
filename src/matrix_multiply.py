#!/usr/bin/env python
# encoding: utf-8
"""
matrix_multiply.py

Created by Kyle Mandli on 2009-03-25.
Copyright (c) 2009 University of Washington. All rights reserved.

Simple matrix multiply script with timings
"""

import time
import numpy as np

N_array = [100,1000,1500,2000,2500,3000]
time_dot = []
time_loop = []

# Python version
for N in N_array:
    A = np.random.rand(N,N)
    B = np.random.rand(N,N)
    
    # Dot method
    print "Dot method ",N
    C = np.zeros((N,N))
    timing = time.time()
    C = np.dot(A,B)
    time_dot.append(time.time() - timing)
    
    # Loop method
    if N <= 1000:
        print "Loop method ",N
        C = np.zeros((N,N))
        timing = time.time()
        for i in xrange(0,N):
            for j in xrange(0,N):
                for k in xrange(0,N):
                    C[i,j] = C[i,j] + A[i,k] * B[k,j]
    
        time_loop.append(time.time() - timing)

print "Dot times = ", time_dot
print "Loop times =", time_loop
