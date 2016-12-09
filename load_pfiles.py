import re
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
from load_sim_constants import load_sim_constants

# ------------------------------------------------
# Loads all pN, pS files in a directory, sorted by frequency
# returns: pN, pS (nFreqs x nTimes x nEnergies numpy arrays)
#       This version for files with L-shell in name
# ------------------------------------------------
def load_pfiles_L(directory, L, fmin=0, fmax=10e10):
    # Get files, frequencies:
    d = os.listdir(directory)


    freqs_pN = []
    freqs_pS = []
    l_pN = []
    l_pS = []

    p = re.compile("\d+")
    for s in d:
        if s.endswith(".dat"):
            if s.startswith('pN'):
                tmp = p.findall(s)
                freqs_pN.append(int(tmp[0]))
                l_pN.append(float(tmp[1]))

            elif s.startswith('pS'):
                tmp = p.findall(s)
                freqs_pS.append(int(tmp[0]))
                l_pS.append(float(tmp[1]))

    freqs_pN.sort()
    freqs_pS.sort()

    # print l_pN
    assert freqs_pN == freqs_pS, "frequency mismatch!"

    freqs = freqs_pN

    #Pre-allocate
    pN = np.zeros((len(freqs), sc.NUM_E, sc.NUM_STEPS))
    pS = np.zeros((len(freqs), sc.NUM_E, sc.NUM_STEPS))


    # Binary files -- little endian, four-byte floats
    dt = np.dtype('<f4')

    # for f, l in zip(freqs, l_pN):
    for f_ind, f in enumerate([ff for ff in freqs if ff > fmin and ff < fmax]):

        # for binary-formatted files
        tmp_N = np.fromfile(os.path.join(directory,'pN%d_%g.dat'%(f,L)),dtype=np.dtype('<f4'))

        if (np.shape(tmp_N)[0]==sc.NUM_E*sc.NUM_STEPS):
            pN[f_ind, :,:] = tmp_N.reshape(sc.NUM_E, sc.NUM_STEPS, order='c')
        else:
            print "no N data at %d"%f 
            
        tmp_S = np.fromfile(os.path.join(directory,'pS%d_%g.dat'%(f,L)),dtype=np.dtype('<f4'))
        if (np.shape(tmp_S)[0]==sc.NUM_E*sc.NUM_STEPS):
            pS[f_ind, :,:] = tmp_S.reshape(sc.NUM_E, sc.NUM_STEPS, order='c')
        else:
            print "no S data at %d"%f

        # For ASCII-formatted files
#         pN.append(np.loadtxt(os.path.join(directory,"pN%d_%d.dat"%(f,L))))
#         pS.append(np.loadtxt(os.path.join(directory,"pS%d_%d.dat"%(f,L))))

    return pN, pS




import re
from load_sim_constants import load_sim_constants

# ------------------------------------------------
# Loads all pN, pS files in a directory, sorted by frequency
# returns: pN, pS (nFreqs x nTimes x nEnergies numpy arrays)
# ------------------------------------------------
def load_pfiles_latlon(directory, out_lat, out_lon, sc, fmin=0, fmax=10e10):
    # Get files, frequencies:
    d = os.listdir(directory)


    freqs_pN = []
    freqs_pS = []
    lat_pN = []
    lat_pS = []
    lon_pN = []
    lon_pS = []

    p = re.compile("\d+")
    for s in d:
        if s.endswith(".dat"):
            if s.startswith('pN'):
                tmp = p.findall(s)
                freqs_pN.append(int(tmp[2]))
                lat_pN.append(float(tmp[0]))
                lon_pN.append(float(tmp[1]))

            elif s.startswith('pS'):
                tmp = p.findall(s)
                freqs_pS.append(int(tmp[2]))
                lat_pS.append(float(tmp[0]))
                lon_pS.append(float(tmp[1]))

    freqs_pN.sort()
    freqs_pS.sort()

    # print l_pN
    assert freqs_pN == freqs_pS, "frequency mismatch!"

    freqs = np.unique(freqs_pN)

    #Pre-allocate
    pN = np.zeros((len(freqs), sc.NUM_E, sc.NUM_STEPS))
    pS = np.zeros((len(freqs), sc.NUM_E, sc.NUM_STEPS))


    # Binary files -- little endian, four-byte floats
    dt = np.dtype('<f4')

    # for f, l in zip(freqs, l_pN):
    for f_ind, f in enumerate([ff for ff in freqs if ff > fmin and ff < fmax]):
        # for binary-formatted files
        tmp_N = np.fromfile(os.path.join(directory,'pN_%g_%g_%d.dat'%(out_lat, out_lon, f)),dtype=np.dtype('<f4'))

        if (np.shape(tmp_N)[0]==sc.NUM_E*sc.NUM_STEPS):
            pN[f_ind, :,:] = tmp_N.reshape(sc.NUM_E, sc.NUM_STEPS, order='c')
        else:
            print "no N data at %d"%f 
            
        tmp_S = np.fromfile(os.path.join(directory,'pS_%g_%g_%d.dat'%(out_lat, out_lon, f)),dtype=np.dtype('<f4'))
        if (np.shape(tmp_S)[0]==sc.NUM_E*sc.NUM_STEPS):
            pS[f_ind, :,:] = tmp_S.reshape(sc.NUM_E, sc.NUM_STEPS, order='c')
        else:
            print "no S data at %d"%f

        # For ASCII-formatted files
#         pN.append(np.loadtxt(os.path.join(directory,"pN%d_%d.dat"%(f,L))))
#         pS.append(np.loadtxt(os.path.join(directory,"pS%d_%d.dat"%(f,L))))

    return pN, pS
