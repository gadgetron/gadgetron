# -*- coding: utf-8 -*-
#%% imports
import os
import sys
import ismrmrd
import ismrmrd.xsd
import numpy as np
from ismrmrdtools import show

sys.path.append(os.environ['GADGETRON_HOME'] + '/share/gadgetron/python')

def gadget_wait_function(first_gadget):
    g = first_gadget;
    while (g):
        g.wait()
        g = g.next_gadget

def gadget_config(first_gadget, conf):
    g = first_gadget;
    while (g):
        g.process_config(conf)
        g = g.next_gadget

def get_last_gadget(first_gadget):
    g = first_gadget
    while (g.next_gadget):
        g = g.next_gadget
    return g

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print "Usage: " + sys.argv[0] + " <gadgetron_python_chain.py> <ismrmrd_out.h5> <ismrmrd_out.h5>"
        raise "Invalid number of arguments."

    python_function_file = sys.argv[1]
    filename = sys.argv[2]
    filename_out = sys.argv[3]

    if not os.path.isfile(python_function_file):
        print("%s is not a valid file" % python_function_file)
        raise SystemExit
    
    namespace = {}
    execfile(python_function_file, namespace)
    globals().update(namespace)

    g0 = define_gadget_chain() #Call function from imported file

    #%% Load file
    if not os.path.isfile(filename):
        print("%s is not a valid file" % filename)
        raise SystemExit
    dset = ismrmrd.Dataset(filename, 'dataset', create_if_needed=False)

    #%% Send in data
    #First ISMRMRD XML header
    gadget_config(g0,dset.read_xml_header())
            
    # Loop through the rest of the acquisitions and stuff
    for acqnum in range(0,dset.number_of_acquisitions()):
        acq = dset.read_acquisition(acqnum)
        g0.process(acq.getHead(),acq.data.astype('complex64'))

    #%%
    gadget_wait_function(g0)

    res = get_last_gadget(g0).get_results()
    print "Received " + str(len(res)) + " result items"

    out_dset = ismrmrd.Dataset(filename_out, "out")
    for o in res:
        print "Appending image to out file"
        img = ismrmrd.Image(head=o[0])
        img.data.ravel()[:] = o[1].ravel()[:] 
        out_dset.append_image("image_%d" % img.image_series_index, img)

    
