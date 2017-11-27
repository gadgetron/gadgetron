import ismrmrd
import ismrmrd.xsd
import numpy as np
# from ismrmrdtools import transform, coils, grappa
from gadgetron import Gadget
import copy 
import math

class CineContouring(Gadget):

    def process_config(self, cfg):
        print("Process config of cine contouring ... ")
        self.images = []
        self.headers = []
        self.metas = []

        self.header = ismrmrd.xsd.CreateFromDocument(cfg)

        # first encoding space
        self.enc = self.header.encoding[0]

        #Parallel imaging factor
        self.acc_factor = self.enc.parallelImaging.accelerationFactor.kspace_encoding_step_1

        self.slc = self.enc.encodingLimits.slice.maximum+1
        self.phs = self.enc.encodingLimits.phase.maximum+1
        self.phs_retro = int(self.params["phase"])

        print("CineContouring, maximal number of slice ", self.slc)

    def process(self, header, image, metadata=None):

        print("Receiving image, phase ", header.phase, ", slice ", header.slice)

        # buffer incoming images
        self.images.append(image)
        self.headers.append(header)
        if metadata is not None:
            # deserialize metadata
            curr_meta = ismrmrd.Meta.deserialize(metadata)
            self.metas.append(curr_meta)
            print(metadata)

        # if all images are received
        if header.slice<self.slc-1 or header.phase<self.phs_retro-1:
            self.put_next(header,image,curr_meta)
            return 0

        # send out the last image
        self.put_next(header,image,curr_meta)

        # enough images received
        print("Sufficient images are received ... ")
        print(len(self.headers))

        # set the Meta

        return 0
