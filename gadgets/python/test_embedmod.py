import numpy as np
import GadgetronPythonMRI

myVariable = 1
myGadgetReference = 0

def set_GadgetReference(gadref):
    global myGadgetReference
    myGadgetReference = gadref


def simple_function(acq, data):
    global myVariable
    global myGadgetReference

    #print "Received line: " + str(acq.idx.line)
    #print "Received types: " + str(type(data))
    #print "Hello"

    myVariable = myVariable + 1

    if (type(myGadgetReference) == GadgetronPythonMRI.GadgetReference):
        t2 = GadgetronPythonMRI.GadgetMessageAcquisition()
        t2.idx.line = 100
        arr = np.ndarray(shape=(2,2), dtype=float, order='F')
        myGadgetReference.return_acquisition(acq, data)

    #if (myVariable > 3):
    #    return_val = np.ndarray(shape=(2,2), dtype=float, order='F')
    #    gadgetron.result_ready(return_val)

    #print "Current sum: " + str(myVariable)


