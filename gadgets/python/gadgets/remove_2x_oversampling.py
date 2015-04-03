import numpy as np
from ismrmrdtools import transform
from gadgetron import Gadget

class Remove2xOversampling(Gadget):

    def process_config(self, conf):
        #print "remove 2x oversampling: Configuration received"
        #print str(conf)
        return

    def process(self, acq, data):
        orig_size = list(data.shape);
        data2 = data.reshape([(data.size/data.shape[data.ndim-1]), data.shape[data.ndim-1]])
        new_length = data2.shape[1]>>1
        data2 = transform.transform_image_to_kspace(transform.transform_kspace_to_image(data2,dim=(1,))[:,(0+(new_length>>1)):(new_length+(new_length>>1))],dim=(1,))
        orig_size[data.ndim-1] = new_length
        data2.reshape(tuple(orig_size))
        acq.samples = new_length

        self.put_next(acq,data2)
        return 0
