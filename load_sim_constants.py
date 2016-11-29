# Load simulation constants from "consts.h" into Python.
# use:
# sc = load_sim_constants(path/to/const.h)


import numpy as np

class consts_struct(object):
    # pass
    def __init__(self):
        params = ['Q_EL','M_EL','E_EL','EPS0','C','Z0','R_E',
                  'H_MAGNETO', 'H_IONO', 'A','B','H_E','P_DIST','Q_DIST',
                  'AN_CM_DIST','V0_DIST','M_RES','E_MIN','E_MAX','NUM_E',
                  'SQUARE','EALimS',
                  'EALimN','EAIncr','dL0','DF','DT','WAVE_PWR_THRESH',
                  'T_MAX','T_STEP','NUM_STEPS','RES_DT','RES_FINT','LK','I0']
        for p in params:
            exec 'self.%s = None'%p


def load_sim_constants(directory, old_format = False):
    params = ['Q_EL','M_EL','E_EL','EPS0','C','Z0','R_E',
              'H_MAGNETO', 'H_IONO', 'A','B','H_E','P_DIST','Q_DIST',
              'AN_CM_DIST','V0_DIST','M_RES','E_MIN','E_MAX','NUM_E',
              'SQUARE','EALimS',
              'EALimN','EAIncr','dL0','DF','DT','WAVE_PWR_THRESH',
              'T_MAX','T_STEP','NUM_STEPS','RES_DT','RES_FINT','LK','I0']

    sc = consts_struct()
    with open(directory,'r+') as file:
        for line in iter(file):
            tmp = line.split()
            if len(tmp) >= 2:
                if tmp[0] == '#define':
                    if tmp[1] in params:
                        try:
                            tmpstr = 'sc.%s = %s'%(tmp[1], tmp[2])
                            exec tmpstr
                            # print "succeeded: ", tmpstr
                        except:
                            print "failed:", tmpstr

    # Do the terms we know it won't pick up:

    if old_format:
        try:
            sc.MU0 = np.pi*4e-7
            sc.T_MAX = sc.RES_FINT
            sc.T_STEP= sc.RES_DT
            sc.NUM_STEPS = np.floor(sc.T_MAX/sc.T_STEP)

            sc.E_EXP_BOT = np.log10(sc.E_MIN)
            sc.E_EXP_TOP = np.log10(sc.E_MAX)
            sc.DE_EXP = ( (sc.E_EXP_TOP - sc.E_EXP_BOT)/sc.NUM_E)

            # Generate energy and velocity arrays
            sc.E_tot_arr = pow(10,sc.E_EXP_BOT + sc.DE_EXP*np.arange(0,sc.NUM_E))
            sc.v_tot_arr = sc.C*np.sqrt(1 - pow(sc.E_EL/(sc.E_EL + sc.E_tot_arr),2))

        except:
            print "...failed (old format)"
    else:        
        try:
            sc.MU0 = np.pi*4e-7

            sc.T_STEP = (1.0*((1.0*sc.T_MAX)/sc.NUM_STEPS))
            sc.E_EXP_BOT = np.log10(sc.E_MIN)
            sc.E_EXP_TOP = np.log10(sc.E_MAX)
            sc.DE_EXP = ( (sc.E_EXP_TOP - sc.E_EXP_BOT)/sc.NUM_E)

            # Generate energy and velocity arrays
            sc.E_tot_arr = pow(10,sc.E_EXP_BOT + sc.DE_EXP*np.arange(0,sc.NUM_E))
            sc.v_tot_arr = sc.C*np.sqrt(1 - pow(sc.E_EL/(sc.E_EL + sc.E_tot_arr),2))
        except:
            print "...failed (new format)"

    return sc
