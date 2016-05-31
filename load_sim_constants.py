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
                        print "succeeded: ", tmpstr
                    except:
                        print "failed:", tmpstr
    return sc
