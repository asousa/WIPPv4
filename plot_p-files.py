# plot p-files.py

import matplotlib       
matplotlib.use("Agg")       # For headless plotting on subnodes

import numpy as np
import pandas as pd
import pickle
#from build_database import flux_obj
from scipy import interpolate
import matplotlib.pyplot as plt
import os
import itertools
import random
import os
# %matplotlib inline
import re
from load_sim_constants import load_sim_constants


import time
from mpi4py import MPI  # parallel toolbox
import commands
from partition import partition



from load_pfiles import load_pfiles_latlon as load_pfiles
from plots import plot_pN_pS

# Initialize MPI:
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
host = commands.getoutput("hostname")


inlats = np.arange(15, 51, step=5)
outlats = np.arange(20,60)
outlons = [2, 4, 6, 8] #$np.arange(0, 10, step=2)
pwr = -10000
root_dir = '/shared/users/asousa/WIPP/WIPPv4/outputs/agu2016_kp0_v3/'

if rank == 0:
    if not os.path.exists(os.path.join(root_dir, 'figures')):
        os.mkdir(os.path.join(root_dir,'figures'))


    tasklist = [(w,x,y) for w,x,y in itertools.product(inlats, outlats, outlons)]
    # Adjacent frequencies take similar time to complete... shuffle to distribute nicely
    np.random.shuffle(tasklist)
    print "total tasks to do: ", len(tasklist)

    sc = load_sim_constants(os.path.join(root_dir, 'pwr_%d/consts.h'%pwr))

    # Prep output filetree... dammit
    if not os.path.exists(os.path.join(root_dir,'figures')):
        os.mkdir(os.path.join(root_dir, 'figures'))

    for inlat in inlats:
        for outlon in outlons:
            figdir = os.path.join(root_dir, 'figures','in_%d/lon_%d'%(inlat, outlon))
            if not os.path.exists(os.path.join(root_dir, 'figures','in_%d'%inlat)):
                os.mkdir(os.path.join(root_dir,'figures','in_%d'%inlat))
            if not os.path.exists(figdir):
                os.mkdir(figdir)




else:
    tasklist = None
    sc = None

tasklist = comm.bcast(tasklist, root=0)
sc = comm.bcast(sc, root=0)

nTasks = 1.0*len(tasklist)
nProcs = 1.0*comm.Get_size()
nSteps = np.ceil(nTasks/nProcs).astype(int)

chunks = partition(tasklist, nProcs)

if (rank < len(chunks)):
    print "Process %d on host %s, doing %g jobs"%(rank, host, len(chunks[rank]))
    for job in chunks[rank]:
        inlat = job[0]
        outlat = job[1]
        outlon = job[2]

        figdir = os.path.join(root_dir, 'figures','in_%d/lon_%d'%(inlat, outlon))

        # if not os.path.exists(os.path.join(root_dir, 'figures','in_%d'%inlat)):
        #     os.mkdir(os.path.join(root_dir,'figures','in_%d'%inlat))
        # if not os.path.exists(figdir):
        #     os.mkdir(figdir)



        # sc = load_sim_constants(os.path.join(root_dir, 'pwr_%d/consts.h'%pwr))
        pN, pS = load_pfiles(os.path.join(root_dir,'pwr_%d/in_%d/lat_%d'%(pwr, inlat,outlat)), outlat, outlon, sc, fmin=0, fmax=40000, zipped=True)
        print np.shape(pN)
        if len(pN) > 0:
            print "Max N (deg): ", np.max(180./np.pi*np.sqrt(pN))
            print "Max S (deg): ", np.max(180./np.pi*np.sqrt(pS))
            plot_pN_pS(np.sqrt(np.sum(pN, axis=0)), np.sqrt(np.sum(pS, axis=0)), sc)
            plt.suptitle('in %g deg, lat = %g, lon = %g'%(inlat, outlat, outlon))
            plt.savefig(os.path.join(figdir,'p_in%g_lat%g_lon%g.png'%(inlat, outlat, outlon )) )
            plt.close('all')
        else:
            "no files found"


