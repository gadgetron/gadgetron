import GadgetronPythonMRI as g
import kspaceandimage as ki
import numpy as np
import matplotlib.pyplot as plt
import pylab

gadget_ref = g.GadgetReference()
myCounter = 0;
def set_gadget_reference(ref):
    global gadget_ref
    gadget_ref = ref

def config_function(cfg):
    print "Configuration Ignored"

def itok_function(h,im):
    global gadget_ref
    global myCounter
    myCounter = myCounter + 1
    #print str(im.nbytes/im.size)
    #pl = plt.imshow(abs(np.squeeze(im)))
    #res = np.fft.fftn(np.fft.helper.ifftshift(im),axes=[2]);#ki.itok(im)
    #res = ki.itok(im);
    #print "This is call number: " + str(myCounter) + ", " + str(im.shape)
    #print im[0,0,0]
    #res = np.fft.fftn(np.fft.helper.ifftshift(im,axes=None),axes=None)
    #res = im
    #res[:,97:,:] = 0
    #print res[0,0,0]
    #gadget_ref.return_image(h,ki.itok(im))
    gadget_ref.return_image(h,ki.itok(im).astype('complex64'))


