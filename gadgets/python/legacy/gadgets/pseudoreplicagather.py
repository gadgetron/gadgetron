import sys

import ismrmrd
import ismrmrd.xsd
from numpy import *

from gadgetron import Gadget



def calcPseudoreplica(original, pseudoreplicas):
    stdImg = std(pseudoreplicas,axis=pseudoreplicas.ndim-1)
    return stdImg

class PseudoreplicaGather(Gadget):
    def __init__(self, next_gadget=None):
        Gadget.__init__(self,next_gadget)
        self.imageBuffer = None
        self.counter = 0
        self.header = None
        self.enc = None
        self.original = None

    def process_config(self, conf):
        #self.header = ismrmrd.xsd.CreateFromDocument(conf)
        #self.enc = self.header.encoding[0]
        print(self.params, file=sys.stderr)
        self.repetitions = int(self.params["repetitions"])

    def process(self, header, img,*args):
        if self.imageBuffer is None:
            s = shape(img)
            s2 = s + (self.repetitions,) 
            self.imageBuffer = zeros(s2,dtype=complex64)

        if self.counter == 0: #First image is without added noise.
            self.original = img
        else:
            self.imageBuffer[...,self.counter - 1] = img

        print(f"Counter: {self.counter} Repetitions: {self.repetitions}")
        if (self.counter == self.repetitions):
            img_head = header
            img_head.data_type = ismrmrd.DATATYPE_FLOAT
            pseudoreplica = calcPseudoreplica(self.original,self.imageBuffer)
            #Return image to Gadgetron
            self.counter = 0
            self.imageBuffer = None
            self.original = None
            print(f"Putting image on stream", file=sys.stderr)
            self.put_next(img_head,pseudoreplica.astype('float32'),*args)
        else:
            self.counter += 1
   
        #print "Returning to Gadgetron"
        return 0 #Everything OK
