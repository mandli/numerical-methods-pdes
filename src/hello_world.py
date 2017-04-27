#!/usr/bin/env python
# encoding: utf-8

"""MPI Hello World from Python"""

from __future__ import absolute_import
from __future__ import print_function

from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

print('Hello from Process number %s of %s processes.' % (rank + 1, size))
