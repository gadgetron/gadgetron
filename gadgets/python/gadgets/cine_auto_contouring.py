import ismrmrd
import ismrmrd.xsd
import numpy as np
# from ismrmrdtools import transform, coils, grappa
from gadgetron import Gadget
import copy 
import math

class CineContouring(Gadget):

    def __init__(self):
        self.images = []
        self.headers = []
        self.metas = []

    def process_config(self, cfg):
        print("Process config of cine contouring ... ")
        self.header = ismrmrd.xsd.CreateFromDocument(cfg)

        # first encoding space
        self.enc = self.header.encoding[0]

        #Parallel imaging factor
        self.acc_factor = self.enc.parallelImaging.accelerationFactor.kspace_encoding_step_1

        slc = self.enc.encodingLimits.slice.maximum+1
        phs = self.enc.encodingLimits.phase.maximum+1
        phs_retro = int(self.params["phase"])

        print("CineContouring, maximal number of slice ", slc)

    def process(self, header, image, metadata=None):

        # buffer incoming images
        self.images.append(image)
        self.headers.append(header)
        if metadata is not None:
            self.metas.append(metadata)
            # deserialize metadata
            curr_meta = ismrmrd.Meta.deserialize(metadata)
            print(metadata)

        # if all images are received
        if header.slice<slc-1:
            return 0

        if header.phase<phase_retro-1:
            return 0

        # enough images received
        print("Sufficient images are received ... ")
        print(header)

        # set the Meta
        for x in self.metas:
            


        return 0
