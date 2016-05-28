import os
import shutil
import datetime
import random
import string
import time
import numpy as np
import commands

import subprocess

from mpi4py import MPI


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
host = commands.getoutput("hostname")

nProcs = comm.Get_size()

print "We have %d processes available"%(nProcs)
