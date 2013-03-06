import numpy as np
import GadgetronPythonMRI as g
import kspaceandimage as ki
import libxml2

myLocalGadgetReference = g.GadgetReference()
myBuffer = 0
myParameters = 0
myCounter = 1;
mySeries = 1;

def set_gadget_reference(gadref):
    global myLocalGadgetReference
    myLocalGadgetReference = gadref

def config_function(conf):
    global myBuffer
    global myParameters

    myParameters = dict()

    doc = libxml2.parseDoc(str(conf))
    context = doc.xpathNewContext()
    context.xpathRegisterNs("ismrm", "http://www.ismrm.org/ISMRMRD")
    myParameters["matrix_x"] = int((context.xpathEval("/ismrm:ismrmrdHeader/ismrm:encoding/ismrm:encodedSpace/ismrm:matrixSize/ismrm:x")[0]).content)
    myParameters["matrix_y"] = int((context.xpathEval("/ismrm:ismrmrdHeader/ismrm:encoding/ismrm:encodedSpace/ismrm:matrixSize/ismrm:y")[0]).content)
    myParameters["matrix_z"] = int((context.xpathEval("/ismrm:ismrmrdHeader/ismrm:encoding/ismrm:encodedSpace/ismrm:matrixSize/ismrm:z")[0]).content)
    myParameters["channels"] = int((context.xpathEval("/ismrm:ismrmrdHeader/ismrm:acquisitionSystemInformation/ismrm:receiverChannels")[0]).content)
    myParameters["slices"] = int((context.xpathEval("/ismrm:ismrmrdHeader/ismrm:encoding/ismrm:encodingLimits/ismrm:slice/ismrm:maximum")[0]).content)+1
    myParameters["center_line"] = int((context.xpathEval("/ismrm:ismrmrdHeader/ismrm:encoding/ismrm:encodingLimits/ismrm:kspace_encoding_step_1/ismrm:center")[0]).content)

    myBuffer = (np.zeros((myParameters["channels"],myParameters["slices"],myParameters["matrix_z"],myParameters["matrix_y"],(myParameters["matrix_x"]>>1)))).astype('complex64')

def recon_function(acq, data):
    global myLocalGadgetReference
    global myBuffer
    global myParameters
    global myCounter
    global mySeries

    line_offset = (myParameters["matrix_y"]>>1)-myParameters["center_line"];
    myBuffer[:,acq.idx.slice,acq.idx.kspace_encode_step_2,acq.idx.kspace_encode_step_1+line_offset,:] = data
    
    if (acq.flags & (1<<7)): #Is this the last scan in slice
        image = ki.ktoi(myBuffer,(2,3,4))
        image = image * np.product(image.shape)*100 #Scaling for the scanner
        #Create a new image header and transfer value
        img_head = g.ImageHeader()
        img_head.channels = acq.active_channels
        img_head.slice = acq.idx.slice
        g.img_set_matrix_size(img_head, 0, myBuffer.shape[4])
        g.img_set_matrix_size(img_head, 1, myBuffer.shape[3])
        g.img_set_matrix_size(img_head, 2, myBuffer.shape[2])
        g.img_set_position(img_head, 0,g.acq_get_position(acq,0))
        g.img_set_position(img_head, 1,g.acq_get_position(acq,1))
        g.img_set_position(img_head, 2,g.acq_get_position(acq,2))
        g.img_set_read_dir(img_head, 0, g.acq_get_read_dir(acq, 0))
        g.img_set_read_dir(img_head, 1, g.acq_get_read_dir(acq, 1))
        g.img_set_read_dir(img_head, 2, g.acq_get_read_dir(acq, 2))
        g.img_set_phase_dir(img_head, 0, g.acq_get_phase_dir(acq, 0))
        g.img_set_phase_dir(img_head, 1, g.acq_get_phase_dir(acq, 1))
        g.img_set_phase_dir(img_head, 2, g.acq_get_phase_dir(acq, 2))
        g.img_set_slice_dir(img_head, 0, g.acq_get_slice_dir(acq, 0))
        g.img_set_slice_dir(img_head, 1, g.acq_get_slice_dir(acq, 1))
        g.img_set_slice_dir(img_head, 2, g.acq_get_slice_dir(acq, 2))
	g.img_set_patient_table_position(img_head, 0, g.acq_get_patient_table_position(acq,0))
	g.img_set_patient_table_position(img_head, 1, g.acq_get_patient_table_position(acq,1))
	g.img_set_patient_table_position(img_head, 2, g.acq_get_patient_table_position(acq,2))
        img_head.acquisition_time_stamp = acq.acquisition_time_stamp
	img_head.image_index = myCounter;
	img_head.image_series_index = mySeries;

	myCounter = myCounter + 1
	if (myCounter > 5):
		mySeries = mySeries + 1
		myCounter = 1

        #Return image to Gadgetron
	return myLocalGadgetReference.return_image(img_head,image.astype('complex64'))
	
    #print "Returning to Gadgetron"
    return 0 #Everything OK

