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

tag = 42
N = 10

if rank == 0:
    data = numpy.arange(N, dtype=numpy.float64)
    print("Process %s note = %s" % (rank, data))
    # Note here that MPI datatype discovery is automatic
    comm.Send(data, dest=rank + 1, tag=tag)

elif rank < size - 1:
    data = numpy.empty(N, dtype=numpy.float64)
    comm.Recv(data, source=rank - 1, tag=tag)

    print("Process %s note = %s" % (rank, data))

    comm.Send(data, dest=rank + 1, tag=tag)

elif rank == size - 1:
    data = numpy.empty(N, dtype=numpy.float64)
    comm.Recv(data, source=rank - 1, tag=tag)

    print("Process %s note = %s" % (rank, data))

else:
    raise Exception("Invalid rank.")
