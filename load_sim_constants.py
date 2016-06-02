# Load simulation constants from "consts.h" into Python.
# use:
# sc = load_sim_constants(path/to/const.h)


import numpy as np

class empty_struct(object):
    pass

def load_sim_constants(directory):

    sc = empty_struct
    with open(directory,'r+') as file:
        for line in iter(file):
            tmp = line.split()
            if len(tmp) >= 2:
                if tmp[0] == '#define':
                    try:
                        tmpstr = 'sc.%s = %s'%(tmp[1], tmp[2])
                        exec tmpstr
                        # print "succeeded: ", tmpstr
                    except:
                        print "failed:", tmpstr

    # Do the terms we know it won't pick up:

    try:
        sc.T_STEP = (1.0*((1.0*sc.T_MAX)/sc.NUM_STEPS))
        sc.E_EXP_BOT = np.log10(sc.E_MIN)
        sc.E_EXP_TOP = np.log10(sc.E_MAX)
        sc.DE_EXP = ( (sc.E_EXP_TOP - sc.E_EXP_BOT)/sc.NUM_E)

        # Generate energy and velocity arrays
        sc.E_tot_arr = pow(10,sc.E_EXP_BOT + sc.DE_EXP*np.arange(0,sc.NUM_E))
        sc.v_tot_arr = sc.C*np.sqrt(1 - pow(sc.E_EL/(sc.E_EL + sc.E_tot_arr),2))
    except:
        print "...shrug."

    return sc
