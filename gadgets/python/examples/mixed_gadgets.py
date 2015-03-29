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
import threading
import time

#%% Setup gadgets
g4 = WrapperGadget("gadgetron_mricore","ExtractGadget", next_gadget=None)
g3 = RMSCoilCombine(g4)
g2 = AccumulateAndRecon(g3)
g1 = Remove2xOversampling(g2)
g0 = WrapperGadget("gadgetron_mricore","NoiseAdjustGadget",next_gadget=g1)


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
    

#%% Load file
filename = 'testdata.h5'
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

print g4.get_results()
