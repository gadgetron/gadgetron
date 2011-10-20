import GadgetronPythonMRI as g
import numpy as np

gadget_ref = g.GadgetReference()

def set_gadget_reference(ref):
    global gadget_ref
    gadget_ref = ref

def config_function(cfg):
    print "RMS Coil Combine, Config ignored"

def recon_function(h,im):
    global gadget_ref
    combined_image = np.sqrt(np.sum(np.square(np.abs(im)),axis=0))
    h.channels = 1
    gadget_ref.return_image(h,combined_image.astype('complex64'))


