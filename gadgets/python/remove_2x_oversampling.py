import numpy as np
import kspaceandimage as ki
import libxml2
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
        data2 = ki.itok(ki.ktoi(data2,[1])[:,(0+(new_length>>1)):(new_length+(new_length>>1))],[1])
        orig_size[data.ndim-1] = new_length
        data2.reshape(tuple(orig_size))
        acq.samples = new_length

        self.put_next(acq,data2)
        return 0
