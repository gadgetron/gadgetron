import numpy as np
import ismrmrd
import ismrmrd.xsd

from gadgetron import Gadget
from gadgetron.util.cfft import cifftn


class AccumulateAndRecon(Gadget):

    def __init__(self, next_gadget=None):
        Gadget.__init__(self, next_gadget)
        self.myBuffer = None
        self.myCounter = 1
        self.mySeries = 1
        self.header = None
        self.enc = None

    def process_config(self, conf):
        self.header = ismrmrd.xsd.CreateFromDocument(conf)
        self.enc = self.header.encoding[0]

    def process(self, acq, data, *args):
        if self.myBuffer is None:
            channels = acq.active_channels
            if self.enc.encodingLimits.slice != None:
                nslices = self.enc.encodingLimits.slice.maximum + 1
            else:
                nslices = 1
            eNz = self.enc.encodedSpace.matrixSize.z
            eNy = self.enc.encodedSpace.matrixSize.y
            eNx = self.enc.encodedSpace.matrixSize.x
        
            self.myBuffer = np.zeros(( int(eNx/2),eNy,eNz,nslices,channels),dtype=np.complex64)

        line_offset = self.enc.encodedSpace.matrixSize.y/2 - self.enc.encodingLimits.kspace_encoding_step_1.center             
        self.myBuffer[:,int(acq.idx.kspace_encode_step_1+line_offset), int(acq.idx.kspace_encode_step_2), int(acq.idx.slice),:] = data

        if (acq.flags & (1<<7)): #Is this the last scan in slice
            image = cifftn(self.myBuffer, axes=(0, 1, 2))
            image = image * np.product(image.shape)*100 #Scaling for the scanner
            #Create a new image header and transfer value
            img_head = ismrmrd.ImageHeader()
            img_head.version = 1
            img_head.measurement_uid = acq.measurement_uid
            img_head.field_of_view = (
                self.enc.reconSpace.fieldOfView_mm.x,
                self.enc.reconSpace.fieldOfView_mm.y,
                self.enc.reconSpace.fieldOfView_mm.z
            )
            img_head.channels = acq.active_channels
            img_head.matrix_size[0] = self.myBuffer.shape[0]
            img_head.matrix_size[1] = self.myBuffer.shape[1]
            img_head.matrix_size[2] = self.myBuffer.shape[2]
            img_head.position = acq.position
            img_head.read_dir = acq.read_dir
            img_head.phase_dir = acq.phase_dir
            img_head.slice_dir = acq.slice_dir
            img_head.patient_table_position = acq.patient_table_position
            img_head.average = acq.idx.average
            img_head.slice = acq.idx.slice
            img_head.contrast = acq.idx.contrast
            img_head.phase = acq.idx.phase
            img_head.repetition = acq.idx.repetition
            img_head.set = acq.idx.set
            img_head.acquisition_time_stamp = acq.acquisition_time_stamp
            img_head.physiology_time_stamp = acq.physiology_time_stamp
            img_head.image_index = self.myCounter
            img_head.image_series_index = self.mySeries
            img_head.data_type = ismrmrd.DATATYPE_CXFLOAT
            self.myCounter += 1
            if self.myCounter > 5:
                    self.mySeries += 1
                    self.myCounter = 1

            #Return image to Gadgetron
            self.put_next(img_head, image.astype('complex64'), *args)
            
        return 0 #Everything OK
