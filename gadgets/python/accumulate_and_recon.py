import numpy as np
import GadgetronPythonMRI as g


myLocalGadgetReference = g.GadgetReference()

def set_gadget_reference(gadref):
    global myGadgetReference
    myLocalGadgetReference = gadref

def config_function(conf):
    print "Configuration received"
    print str(conf)

def recon_function(acq, data):
    global myLocalGadgetReference

    print "Received line: " + str(acq.idx.line)


