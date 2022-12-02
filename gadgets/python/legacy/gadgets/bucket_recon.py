import itertools
import numpy as np
import sys 

from gadgetron import Gadget
from gadgetron.util.cfft import cifftn

import ismrmrd
import ismrmrd.xsd



class BucketRecon(Gadget):
    def __init__(self, next_gadget=None):
        Gadget.__init__(self,next_gadget)
        self.header = None
        self.enc = None

        self.image_indices = itertools.count(1)

    def process_config(self, conf):
        self.header = ismrmrd.xsd.CreateFromDocument(conf)
        self.enc = self.header.encoding[0]

    def process(self, recondata):
        print(np.shape(recondata[0].data.data), file=sys.stderr)

        image = cifftn(recondata[0].data.data, axes=(0, 1, 2))
        image = np.reshape(image,(image.shape[0],image.shape[1],image.shape[2],image.shape[3]))

        #Create a new image header and transfer value
        acq = np.ravel(recondata[0].data.headers)[0]

        img_head = ismrmrd.ImageHeader()
        img_head.version = 1
        img_head.measurement_uid = acq.measurement_uid
        img_head.channels = acq.active_channels
        img_head.slice = acq.idx.slice
        img_head.matrix_size = (image.shape[0],image.shape[1],image.shape[2])
        img_head.field_of_view = (
            self.enc.reconSpace.fieldOfView_mm.x,
            self.enc.reconSpace.fieldOfView_mm.y,
            self.enc.reconSpace.fieldOfView_mm.z
        )
        img_head.position = acq.position
        img_head.read_dir = acq.read_dir
        img_head.phase_dir = acq.phase_dir
        img_head.slice_dir = acq.slice_dir
        img_head.patient_table_position = acq.patient_table_position
        img_head.acquisition_time_stamp = acq.acquisition_time_stamp
        img_head.average = acq.idx.average
        img_head.slice = acq.idx.slice
        img_head.contrast = acq.idx.contrast
        img_head.phase = acq.idx.phase
        img_head.repetition = acq.idx.repetition
        img_head.set = acq.idx.set
        img_head.image_index = next(self.image_indices)
        img_head.image_series_index = 0
        img_head.data_type = ismrmrd.DATATYPE_CXFLOAT

        #Return image to Gadgetron
        self.put_next(img_head,image)
        print(f"Slice {img_head.slice}", file=sys.stderr)
        print("----------------------------------------------", file=sys.stderr)
        return 0    
        #print "Returning to Gadgetron"
 
