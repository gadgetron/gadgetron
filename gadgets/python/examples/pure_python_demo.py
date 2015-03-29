# -*- coding: utf-8 -*-
"""
Created on Sun Mar  1 19:59:54 2015

This is a simple example of using Python Gadgets in a standalone Python 
environment. 

In order to run this example, you need the ismrmrd-python API and also the 
ismrmrd-python-tools toolbox. 

These tools depend on h5py (version 2.3 or higher)

To install the tools:

ISMRMRD Python API:

git clone https://github.com/ismrmrd/ismrmrd-python.git
cd ismrmrd-python
sudo python setup.py install

ISMRMRD Python tools:

git clone https://github.com/ismrmrd/ismrmrd-python-tools.git
cd ismrmrd-python-tools
sudo python setup.py install


@author: Michael S. Hansen
"""

#%% imports
import os
import sys
import ismrmrd
import ismrmrd.xsd
import numpy as np
from ismrmrdtools import show
from remove_2x_oversampling import Remove2xOversampling
from accumulate_and_recon import AccumulateAndRecon
from rms_coil_combine import RMSCoilCombine

#%% Setup gadgets
g3 = RMSCoilCombine()
g2 = AccumulateAndRecon(g3)
g1 = Remove2xOversampling(g2)


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
    g1.process(acq.getHead(),acq.data)
    
#%%
#Get result and display    
res = g3.get_results()
show.imshow(np.squeeze(abs(res[0][1])))