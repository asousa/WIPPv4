# run_wipp_MPI.py
#
# Builds and runs the full WIPP code. 
# I can't believe I'm writing this shit yet again
#
# Mmm... spaghetti.
#

import os
import shutil
import datetime
import random
import string
import time
import numpy as np
import itertools
import commands
import subprocess

from mpi4py import MPI

# Frequencies to do:
# freqs = [200,240,289,347,418,502,603,725,872,1048,1259,1514,1819,2187,2629,3160,3798,4565,5487,6596,7928,9530,11455,13769,16550,19893,23912,28742,34549,41528,49916,60000]
freqs_log = np.linspace(np.log10(200), np.log10(60000), 130)
freqs = np.round(pow(10, freqs_log))

# freqs = freqs[0:3]
# freqs = np.arange(200,10000,step=1000)

freq_pairs = zip(freqs[0:], freqs[1:])
# print "freqs:", freqs
# print "Freq_pairs:", freq_pairs

# print "Frequency pairs:",np.shape(freq_pairs)

# L-shells to calculate at:
# 20 to 70 degrees (magnetic), uniform
#L_targ_list = [1.33, 1.37, 1.41, 1.45, 1.49, 1.54, 1.59, 1.64, 1.70, 1.77, 1.84, 1.92, 2.00, 2.09, 2.19, 2.30, 2.42, 2.55, 2.70, 2.86, 3.04, 3.24, 3.46, 3.72, 4.00, 4.32, 4.69, 5.11, 5.60, 6.17, 6.83, 7.61, 8.55]
L_targ_list = [3]
# print "L targets:", np.shape(L_targ)

# Input latitudes to launch from:
# in_lat = 35.0
center_lats = np.arange(40,60,step=10)

# print "Center lats:", np.shape(center_lats)

# Input power to do:
I0_list    = [-100000.0]

# Relevant paths:
root_dir = '/shared/users/asousa/WIPP/WIPPv4/'
code_dir = os.path.join(root_dir,'codesrc')
ray_dir  = '/shared/users/asousa/WIPP/WIPPv4/rays/130f_60s' # Location of ray files
out_dir  = os.path.join(root_dir,'out_full')                # Where to assemble final results
log_dir  = os.path.join(out_dir, 'logs')                    # Where to write log files
local_dir= os.path.join(os.path.expanduser("~"),'scatter_tmp') # Local directory on each node


# if not os.path.exists(local_dir):
#     os.mkdir(local_dir)


# Initialize MPI:
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
host = commands.getoutput("hostname")

# Generate task list:
tasklist = [(w,x,y,z) for w,x,y,z in itertools.product(I0_list, center_lats, L_targ_list, freq_pairs)]

# Segment task list for each node:
nTasks = 1.0*len(tasklist)
nProcs = 1.0*comm.Get_size()
nSteps = np.ceil(nTasks/nProcs).astype(int)

chunks = [tasklist[i:i+nSteps] for i in range(0, len(tasklist), nSteps)]

# See what each node is calling dibs on:
# print chunks[rank]


# Prep output directory tree (root node only):
if rank==0:


  if not os.path.exists(out_dir):
    os.mkdir(out_dir)
  if not os.path.exists(log_dir):
    os.mkdir(log_dir)


  for in_pwr, in_lat, out_L in itertools.product(I0_list, center_lats, L_targ_list):
    pwr_path = os.path.join(out_dir,'pwr_%g'%(in_pwr))
    if not os.path.exists(pwr_path):
      os.mkdir(pwr_path) 

    lat_path = os.path.join(pwr_path,'in_%g'%(in_lat))
    if not os.path.exists(lat_path):
      os.mkdir(lat_path)

    out_L_path = os.path.join(lat_path,'out_%g'%(out_L))
    if not os.path.exists(out_L_path):
      os.mkdir(out_L_path)



# # Compile C code for local machine:
# if not os.path.exists(os.path.join(local_dir,'calc_scattering')):
#   print "---- COMPILING: Node %s ----"%(host)
#   os.system("gcc -o %s/calc_scattering %s/calc_scattering.c -lm"%(local_dir,code_dir))



if (rank < len(chunks)):
  print "Process %d on host %s, doing %g jobs"%(rank, host, len(chunks[rank]))

  for job in chunks[rank]:
    I0 = job[0]
    center_lat = job[1]
    L_targ = job[2]
    f_low  = job[3][0]
    f_high = job[3][1]


    if not os.path.exists(local_dir):
      print "making local dir on %s"%(host)
      os.mkdir(local_dir)

    working_dir = os.path.join(local_dir,'%g_%g_%g_%g'%(np.abs(job[0]),job[1],job[2],job[3][0]))

    if os.path.exists(working_dir):
      os.system("rm -r %s"%(working_dir))

    os.mkdir(working_dir)

    # Move to working directory
    os.chdir(working_dir)

    # Compile C code:
    # print "---- COMPILING: Node %s ----"%(host)
    os.system("gcc -o calc_scattering %s/calc_scattering.c -lm"%(code_dir))

    # Run it:
    job_cmd = "./calc_scattering %s %s %s %s %s %s" %(ray_dir, I0, center_lat, f_low, f_high, L_targ)
    print "Job command: ",job_cmd

    runlog = subprocess.check_output(job_cmd,shell=True)
    # print runlog
    # os.system(job_cmd)
    # write output log:
    file = open(os.path.join(log_dir,'ray_%g_%g_%g_%g.log'%(I0, center_lat, L_targ, f_low)),'w+')
    file.write(runlog)
    file.close()

    # Move compiled files to out dir
    job_out_path = os.path.join(out_dir,"pwr_%g/in_%g/out_%g"%(I0, center_lat, L_targ))
    os.system("mv pN* %s"%(job_out_path))
    os.system("mv pS* %s"%(job_out_path))

    os.chdir(root_dir)

    # Keep it neat, yo
    os.system("rm -r %s"%(working_dir))


  # Remove local directory:
  # os.system('rm -r %s'%(local_dir))
else:
  print "Process %d on host %s has nothing to do..."%(rank, host)
