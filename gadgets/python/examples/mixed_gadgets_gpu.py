#%% imports
import os
import sys

sys.path.append(os.environ['GADGETRON_HOME'] + '/share/gadgetron/python')

import ismrmrd
import ismrmrd.xsd
import numpy as np
from ismrmrdtools import show
from gadgetron import WrapperGadget
import GadgetronPythonMRI as g


  # <gadget>
  #   <name>NoiseAdjust</name>
  #   <dll>gadgetron_mricore</dll>
  #   <classname>NoiseAdjustGadget</classname>
  # </gadget>

g1 = WrapperGadget("gadgetron_mricore","NoiseAdjustGadget")

  
  # <gadget>
  #   <name>PCA</name>
  #   <dll>gadgetron_mricore</dll>
  #   <classname>PCACoilGadget</classname>
  # </gadget>
  
g2 = WrapperGadget("gadgetron_mricore","PCACoilGadget", next_gadget=None)
g1.next_gadget = g2

  # <gadget>
  #   <name>CoilReduction</name>
  #   <dll>gadgetron_mricore</dll>
  #   <classname>CoilReductionGadget</classname>
  #   <property><name>coils_out</name><value>16</value></property>
  # </gadget>

g3 = WrapperGadget("gadgetron_mricore","CoilReductionGadget", next_gadget=None)
g3.set_parameter("CoilReductionGadget","coils_out","16");
g2.next_gadget = g3

  # <gadget>
  #   <name>gpuSpiralSensePrepGadget</name>
  #   <dll>gadgetron_spiral</dll>
  #   <classname>gpuSpiralSensePrepGadget</classname>
  #   <property><name>deviceno</name><value>0</value></property>
  #   <property><name>use_multiframe_grouping</name><value>true</value></property>
  #   <property><name>propagate_csm_from_set</name><value>0</value></property>
  #   <property><name>buffer_convolution_kernel_width</name><value>5.5</value></property>
  #   <property><name>buffer_convolution_oversampling_factor</name><value>1.25</value></property>
  #   <property><name>reconstruction_os_factor_x</name><value>1.5</value></property>
  #   <property><name>reconstruction_os_factor_y</name><value>1.5</value></property>
  # </gadget>
  
  # <gadget>
  #   <name>gpuCgSenseGadget</name>
  #   <dll>gadgetron_gpuparallelmri</dll>
  #   <classname>gpuCgSenseGadget</classname>
  #   <property><name>pass_on_undesired_data</name>  <value>true</value></property>
  #   <property><name>deviceno</name>                <value>0</value></property>
  #   <property><name>setno</name>                   <value>0</value></property>
  #   <property><name>number_of_iterations</name>    <value>10</value></property>
  #   <property><name>cg_limit</name>                <value>1e-6</value></property>
  #   <property><name>oversampling_factor</name>     <value>1.25</value></property>
  #   <property><name>kernel_width</name>            <value>5.5</value></property>
  #   <property><name>kappa</name>                   <value>0.3</value></property>
  # </gadget>

  # <gadget>
  #   <name>gpuCgSenseGadget</name>
  #   <dll>gadgetron_gpuparallelmri</dll>
  #   <classname>gpuCgSenseGadget</classname>
  #   <property><name>pass_on_undesired_data</name>  <value>true</value></property>
  #   <property><name>deviceno</name>                <value>0</value></property>
  #   <property><name>setno</name>                   <value>1</value></property>
  #   <property><name>number_of_iterations</name>    <value>10</value></property>
  #   <property><name>cg_limit</name>                <value>1e-6</value></property>
  #   <property><name>oversampling_factor</name>     <value>1.25</value></property>
  #   <property><name>kernel_width</name>            <value>5.5</value></property>
  #   <property><name>kappa</name>                   <value>0.3</value></property>
  # </gadget>

g4 = WrapperGadget("gadgetron_gpuparallelmri","gpuCgSenseGadget",gadgetname="gpuCgSenseGadget1", next_gadget=None)
g4.prepend_gadget("gadgetron_gpuparallelmri","gpuCgSenseGadget", gadgetname="gpuCgSenseGadget2")
g4.prepend_gadget("gadgetron_spiral","gpuSpiralSensePrepGadget",gadgetname="gpuSpiralSensePrepGadget")

g4.set_parameter("gpuSpiralSensePrepGadget","deviceno","0")
g4.set_parameter("gpuSpiralSensePrepGadget","use_multiframe_grouping","true")
g4.set_parameter("gpuSpiralSensePrepGadget","propagate_csm_from_set","0")
g4.set_parameter("gpuSpiralSensePrepGadget","buffer_convolution_kernel_width","5.5")
g4.set_parameter("gpuSpiralSensePrepGadget","buffer_convolution_oversampling_factor","1.25")
g4.set_parameter("gpuSpiralSensePrepGadget","reconstruction_os_factor_x","1.5")
g4.set_parameter("gpuSpiralSensePrepGadget","reconstruction_os_factor_y","1.5")

g4.set_parameter("gpuCgSenseGadget1","pass_on_undesired_data","true")
g4.set_parameter("gpuCgSenseGadget1","deviceno","0")
g4.set_parameter("gpuCgSenseGadget1","setno","0")
g4.set_parameter("gpuCgSenseGadget1","number_of_iterations","10")
g4.set_parameter("gpuCgSenseGadget1","cg_limit","1e-6")
g4.set_parameter("gpuCgSenseGadget1","oversampling_factor","1.25")
g4.set_parameter("gpuCgSenseGadget1","kernel_width","5.5")
g4.set_parameter("gpuCgSenseGadget1","kappa","0.3")

g4.set_parameter("gpuCgSenseGadget2","pass_on_undesired_data","true")
g4.set_parameter("gpuCgSenseGadget2","deviceno","0") #Think this should be "1"
g4.set_parameter("gpuCgSenseGadget2","setno","1")
g4.set_parameter("gpuCgSenseGadget2","number_of_iterations","10")
g4.set_parameter("gpuCgSenseGadget2","cg_limit","1e-6")
g4.set_parameter("gpuCgSenseGadget2","oversampling_factor","1.25")
g4.set_parameter("gpuCgSenseGadget2","kernel_width","5.5")
g4.set_parameter("gpuCgSenseGadget2","kappa","0.3")

g3.next_gadget = g4

  # <gadget>
  #   <name>PhaseSubtraction</name>
  #   <dll>gadgetron_mricore</dll>
  #   <classname>FlowPhaseSubtractionGadget</classname>
  # </gadget>
  
g5 = WrapperGadget("gadgetron_mricore","FlowPhaseSubtractionGadget", next_gadget=None)

g4.next_gadget = g5
  # <gadget>
  #   <name>MaxwellCorrection</name>
  #   <dll>gadgetron_mricore</dll>
  #   <classname>MaxwellCorrectionGadget</classname>
  # </gadget>

g6 = WrapperGadget("gadgetron_mricore","MaxwellCorrectionGadget", next_gadget=None)

g5.next_gadget = g6;
  
  # <gadget>
  #   <name>Extract</name>
  #   <dll>gadgetron_mricore</dll>
  #   <classname>ExtractGadget</classname>
  #   <property><name>extract_mask</name><value>9</value></property>
  # </gadget>

g7 = WrapperGadget("gadgetron_mricore","ExtractGadget", next_gadget=None)
g7.set_parameter("ExtractGadget","extract_mask","9")

g6.next_gadget = g7

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
filename = '/home/hansenms/temp/simple_spiral.h5'
if not os.path.isfile(filename):
    print("%s is not a valid file" % filename)
    raise SystemExit
dset = ismrmrd.Dataset(filename, 'dataset', create_if_needed=False)

#%% Send in data
#First ISMRMRD XML header
gadget_config(g1,dset.read_xml_header())

# Loop through the rest of the acquisitions and stuff
for acqnum in range(0,dset.number_of_acquisitions()):
    print "Sending in acquisition " + str(acqnum) + " of " + str(dset.number_of_acquisitions())
    acq = dset.read_acquisition(acqnum)
    g1.process(acq.getHead(),acq.data.astype('complex64'))

# #%%
gadget_wait_function(g1)

print g7.get_results()
