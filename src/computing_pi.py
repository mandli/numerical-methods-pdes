#!/usr/bin/env python
# encoding: utf-8

"""Note passing from Python"""

from __future__ import absolute_import
from __future__ import print_function

from mpi4py import MPI
import numpy

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

# This does not seem to work well in Python
# if rank == 0:
#     print("Using %s processes" % size)
#     N = input("  Input N...")

N = 100
N = comm.bcast(N, root=0)

dx = 1.0 / N

points_per_proc = (N + size - 1) / size

if rank == 0:
    print("Points/process = %s" % points_per_proc)

start = int(rank * points_per_proc + 1)
end = int(min((rank + 1) * points_per_proc, N))

print("Process %s will take i = %s through i = %s" % (rank, start, end))

pi_sum_local = 0.0
for i in range(start, end + 1):
    x = (i - 0.5) * dx
    pi_sum_local += 1.0 / (1.0 + x**2)
pi_sum_local = numpy.array(pi_sum_local, dtype='d')

pi_sum = numpy.zeros(1)
comm.Reduce([pi_sum_local, MPI.DOUBLE], [pi_sum, MPI.DOUBLE],
            op=MPI.SUM, root=0)

if rank == 0:
    pi = 4.0 * dx * pi_sum
    print("The approximation to pi is %s" % (pi))
    print("Difference = %s" % (numpy.abs(numpy.pi - pi)))
