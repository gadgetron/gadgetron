import numpy as np
import GadgetronPythonMRI as g
import kspaceandimage as ki

myGadgetReference = g.GadgetReference()

def set_gadget_reference(gadref):
    global myGadgetReference
    myGadgetReference = gadref

def config_function(conf):
    print "Configuration received"
    print str(conf)

def recon_function(acq, data):
    global myVariable
    global myGadgetReference

    orig_size = list(data.shape);
    data2 = data.reshape([(data.size/data.shape[data.ndim-1]), data.shape[data.ndim-1]])
    new_length = data2.shape[1]>>1
    data2 = ki.itok(ki.ktoi(data2,[1])[:,(0+(new_length>>1)):(new_length+(new_length>>1))],[1])
    orig_size[data.ndim-1] = new_length
    data2.reshape(tuple(orig_size))
    acq.samples = new_length
    myGadgetReference.return_acquisition(acq,data2.astype('complex64'))

    

