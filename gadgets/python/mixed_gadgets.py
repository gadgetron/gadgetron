# -*- coding: utf-8 -*-
"""
Created on Tue Mar 10 11:11:39 2015

@author: Michael S.Hansens
"""


#%% imports
import os
import sys

sys.path.append(os.environ['GADGETRON_HOME'] + '/share/gadgetron/python')

import ismrmrd
import ismrmrd.xsd
import numpy as np
from ismrmrdtools import show
from remove_2x_oversampling import Remove2xOversampling
from accumulate_and_recon import AccumulateAndRecon
from rms_coil_combine import RMSCoilCombine
from gadgetron import WrapperGadget
import GadgetronPythonMRI as g


#%% Setup gadgets
#g5 = WrapperGadget("gadgetron_mricore","ImageFinishGadgetFLOAT")
#g4 = WrapperGadget("gadgetron_mricore","ExtractGadget", next_gadget=g5)
g3 = RMSCoilCombine()
g2 = AccumulateAndRecon(g3)
g1 = Remove2xOversampling(g2)
g0 = WrapperGadget("gadgetron_mricore","NoiseAdjustGadget",next_gadget=g1)



#%% Load file
filename = 'testdata.h5'
if not os.path.isfile(filename):
    print("%s is not a valid file" % filename)
    raise SystemExit
dset = ismrmrd.Dataset(filename, 'dataset', create_if_needed=False)

#%% Send in data
#First ISMRMRD XML header
g1.process_config(dset.read_xml_header())
g2.process_config(dset.read_xml_header())
g3.process_config(dset.read_xml_header())

# Loop through the rest of the acquisitions and stuff
for acqnum in range(0,dset.number_of_acquisitions()):
    acq = dset.read_acquisition(acqnum)
    g0.process(acq.getHead(),acq.data.astype('complex64'))

g = g0;
while (g):
    g.wait()
    g = g.next_gadget

print g3.get_results()
#%%
#Get result and display    
#res = g3.get_results()
#show.imshow(np.squeeze(abs(res[0][1])))
