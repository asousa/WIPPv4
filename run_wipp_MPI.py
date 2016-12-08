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
import subprocess # For calling shell commands with returned output
from coordinate_structure import coordinate_structure # My coordinate transformations
from partition import partition
from mpi4py import MPI  # parallel toolbox

# Initialize MPI:
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
host = commands.getoutput("hostname")


# Enable or disable each section
do_scattering = True
do_flux       = False


# Frequencies to do:
# freqs = [200,240,289,347,418,502,603,725,872,1048,1259,1514,1819,2187,2629,3160,3798,4565,5487,6596,7928,9530,11455,13769,16550,19893,23912,28742,34549,41528,49916,60000]
# freqs = [1000, 1100]
#freqs_log = np.linspace(np.log10(200), np.log10(60000), 130)
#freqs = np.round(pow(10, freqs_log))

f1 = 200; f2 = 30000;
num_freqs = 33
flogs = np.linspace(np.log10(f1), np.log10(f2), num_freqs)
freqs = np.round(pow(10, flogs)/10.)*10


freq_pairs = zip(freqs[0:], freqs[1:])

# L-shells to calculate at:
# out_lats = np.linspace(30,70,32)
out_lats = np.arange(20, 60)
cs = coordinate_structure(out_lats, [0], [100],'geomagnetic')
cs.transform_to('L_dipole')
L_targ_list = np.round(100.*cs.L())/100.


# L_targ_list = [2.422]
# L_targ_list = [ 1.33, 1.36, 1.4,  1.45, 1.49, 1.54, 1.59, 1.65,
#                 1.72, 1.78, 1.86, 1.94, 2.03, 2.13, 2.23, 2.35,
#                 2.48, 2.63, 2.79, 2.96, 3.16, 3.38, 3.63, 3.92,
#                 4.24, 4.61, 5.03, 5.53, 6.1,  6.78, 7.58, 8.54 ]
# L_targ_list = [1.33, 1.37, 1.41, 1.45, 1.49, 1.54, 1.59, 1.64, 1.70, 1.77, 1.84, 1.92, 2.00, 2.09, 2.19, 2.30, 2.42, 2.55, 2.70, 2.86, 3.04, 3.24, 3.46, 3.72, 4.00, 4.32, 4.69, 5.11, 5.60, 6.17, 6.83, 7.61, 8.55]

# Input latitudes to launch from:
# in_lat = 35.0
# center_lats = np.arange(5,65,step=5)
#center_lats = [5, 10, 15, 20, 30]
center_lats = [15, 20, 25,      35, 40, 45, 50, 55]
#center_lats = [20, 30, 40, 50]
# Input power to do:
I0_list    = [-10e3]

# Relevant paths:
root_dir = '/shared/users/asousa/WIPP/WIPPv4/'
code_dir = os.path.join(root_dir,'codesrc')
ray_dir  = '/shared/users/asousa/WIPP/WIPPv4/rays/33f_kp0' # Location of ray files
out_dir  = os.path.join(root_dir,'outputs','agu2016_kp0_v2')                 # Where to assemble final results
log_dir  = os.path.join(out_dir, 'logs')                    # Where to write log files
local_dir= os.path.join(os.path.expanduser("~"),'scatter_tmp') # Local directory on each node
# local_dir = '/tmp/scatter_tmp';


if not os.path.exists(local_dir):
  try:
    os.mkdir(local_dir)
  except:
    pass
comm.Barrier()

# ----------------------------- calc_scattering.c ---------------------------------------

if do_scattering:
  # Generate task list:
  if rank == 0:
    tasklist = [(w,x,y,z) for w,x,y,z in itertools.product(I0_list, center_lats, L_targ_list, freq_pairs)]
    # Adjacent frequencies take similar time to complete... shuffle to distribute nicely
    np.random.shuffle(tasklist)
  else:
    tasklist = None

  tasklist = comm.bcast(tasklist, root=0)

  # tasklist = [(w,x,y,z) for w,x,y,z in itertools.product(I0_list, center_lats, L_targ_list, freq_pairs)]

  # Segment task list for each node:
  nTasks = 1.0*len(tasklist)
  nProcs = 1.0*comm.Get_size()
  nSteps = np.ceil(nTasks/nProcs).astype(int)

  # chunks = [tasklist[i:i+nSteps] for i in range(0, len(tasklist), nSteps)]
  chunks = partition(tasklist, nProcs)

  # See what each node is calling dibs on:
  # print chunks[rank]


  # # Prep output directory tree (root node only):
  if rank==0:
    print "----------- SCATTERING ------------"

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

      # Copy the constants file into the first directory, for reference
      os.system("cp %s/consts.h %s"%(code_dir, pwr_path))
      

  comm.Barrier()

  time.sleep(5)
  # Run each set of jobs on current node:
  if (rank < len(chunks)):
    print "Process %d on host %s, doing %g jobs"%(rank, host, len(chunks[rank]))

    # if not os.path.exists(local_dir):
    #   print "making local dir on %s"%(host)
    #   os.mkdir(local_dir)

    # os.chdir(local_dir)

    # # Compile C code:
    # print "---- COMPILING: Node %s ----"%(host)
    # os.system("gcc -o calc_scattering %s/calc_scattering.c -lm"%(code_dir))


    for job in chunks[rank]:
      I0 = job[0]
      center_lat = job[1]
      L_targ = job[2]
      f_low  = job[3][0]
      f_high = job[3][1]

      working_dir = os.path.join(local_dir,'%s_%d_%g_%g_%g'%(host, rank, L_targ, center_lat, f_low))

      # if os.path.exists(working_dir):
      #   os.system("rm %s/*"%(working_dir))
      # else:
      if not os.path.exists(working_dir):
          os.mkdir(working_dir)

      # Move to working directory
      os.chdir(working_dir)

      # Compile C code:
      # print "---- COMPILING: Node %s ----"%(host)
      os.system("gcc -o calc_scattering %s/calc_scattering.c -lm"%(code_dir))

      # Run it:
      job_cmd = "./calc_scattering %s %s %s %s %s %g" %(ray_dir, I0, center_lat, f_low, f_high, L_targ)
      print "%s/%d: Job command: %s"%(host, rank, job_cmd)

      runlog = subprocess.check_output(job_cmd,shell=True)
      # print runlog
      # os.system(job_cmd)
      # write output log:
      file = open(os.path.join(working_dir,'scatter_%g_%g_%g_%g.log'%(I0, center_lat, L_targ, f_low)),'w')
      file.write(runlog)
      file.close()

      # file = open(os.path.join(log_dir,'scatter_%g_%g_%g_%g.log'%(I0, center_lat, L_targ, f_low)),'w+')
      # subprocess.call(job_cmd, shell=True, stdout=file)
      # file.close()

      # Move compiled files to out dir

      job_out_path = os.path.join(out_dir,"pwr_%g/in_%g/out_%g"%(I0, center_lat, L_targ))

      for attempt in range(3):
        try:
          n_filename = "pN%d_%g.dat"%(f_low, L_targ)
          shutil.copy(n_filename, "%s/%s"%(job_out_path, n_filename))
          # os.system("cp pN%d_%g.dat %s"%(f_low, L_targ, job_out_path))
          s_filename = "pS%d_%g.dat"%(f_low, L_targ)
          shutil.copy(s_filename, "%s/%s"%(job_out_path, s_filename))
          # shutil.copy("pS%d_%g.dat %s"%(f_low, L_targ), job_out_path)
          # os.system("cp pS%d_%g.dat %s"%(f_low, L_targ, job_out_path))
          log_filename ="scatter_%g_%g_%g_%g.log"%(I0, center_lat, L_targ, f_low)
          shutil.copy(log_filename, "%s/%s"%(log_dir, log_filename))
          # os.system("cp scatter_%g_%g_%g_%g.log %s"%(I0, center_lat, L_targ, f_low, log_dir))
        except:
          print "[%s/%d:Failed to copy files for (%d, %g) (attempt %d of 3)"%(host, rank, f_low, L_targ, attempt)
          time.sleep(5)
          continue
        else:
          print "[%s/%d:Successfully copied files for (%d, %g)"%(host, rank, f_low, L_targ)
          break
      
      # os.chdir(root_dir)

      time.sleep(5)

      # Keep it neat, yo
      # os.system("rm -r %s"%(working_dir))

      # # pause for a smidge
      # time.sleep(np.random.randint(1,5))


    # Remove local directory:
    # os.system('rm -r %s'%(local_dir))
  else:
    print "Process %d on host %s has nothing to do..."%(rank, host)


  # Wait until all scatter jobs are complete
  comm.Barrier()

  os.system("rm -r %s/*"%local_dir)

  if rank==0:
    print "------------Finished with scattering-------------"



# ----------------------------- calc_flux.c ---------------------------------------
# if do_flux:
#   print "----------- CALCULATING FLUX ------------"
#   # Generate task list:

#   tasklist = [(w,x,y) for w,x,y in itertools.product(I0_list, center_lats, L_targ_list)]

#   # Segment task list for each node:
#   nTasks = 1.0*len(tasklist)
#   nProcs = 1.0*comm.Get_size()
#   nSteps = np.ceil(nTasks/nProcs).astype(int)

#   chunks =  comm.bcast([tasklist[i:i+nSteps] for i in range(0, len(tasklist), nSteps)])

#   if (rank < len(chunks)):

#     print "Process %d on host %s, doing %g jobs"%(rank, host, len(chunks[rank]))


#     if not os.path.exists(local_dir):
#       print "making local dir on %s"%(host)
#       os.mkdir(local_dir)

#     for job in chunks[rank]:
#       I0 = job[0]
#       center_lat = job[1]
#       L_targ = job[2]


#       working_dir = os.path.join(local_dir,'flux_%g_%g_%g/'%(np.abs(I0),center_lat,L_targ))
#       if os.path.exists(working_dir):
#         print "clearing old working directory"
#         os.system("rm -r %s"%(working_dir))

#       os.mkdir(working_dir)
#       os.chdir(working_dir)

#       os.system("gcc -o calc_flux %s/calc_flux.c -lm"%(code_dir))

#       scatter_dir = os.path.join(out_dir,'pwr_%g/in_%g/out_%g'%(I0, center_lat, L_targ))
#       flux_file   = os.path.join(code_dir,'EQFLUXMA.dat')

#       job_cmd = "./calc_flux %s %g %s" %(scatter_dir, job[2], flux_file)
#       print "Job command: ", job_cmd

#       # os.system(job_cmd)
#       # runlog = subprocess.check_output(job_cmd, shell=True)

#       # file = open(os.path.join(log_dir,'flux_%g_%g_%g.log'%(I0, center_lat, L_targ)),'w+')
#       # file.write(runlog)
#       # file.close()

#       file = open(os.path.join(log_dir,'flux_%g_%g_%g.log'%(I0, center_lat, L_targ)),'w+')
#       subprocess.call(job_cmd, shell=True, stdout=file)
#       file.close()

#       print "Completed: %s, %s"%(center_lat, L_targ)

#       # Move completed files over to shared directory:
#       # os.system("mv phi* %s"%(scatter_dir))

#       os.chdir(root_dir)

#       # Keep it neat, yo
#       os.system("rm -r %s"%(working_dir)) 

#   else:
#     print "Process %d on host %s has nothing to do..."%(rank, host)



