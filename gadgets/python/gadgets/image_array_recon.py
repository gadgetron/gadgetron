import ismrmrd
import ismrmrd.xsd
import numpy as np
# from ismrmrdtools import transform, coils, grappa
from gadgetron import Gadget,IsmrmrdImageArray
import copy 
import math

class ImageArrayRecon(Gadget):

    def process_config(self, cfg):
        print("ImageArrayRecon, process config ... ")
        self.header = ismrmrd.xsd.CreateFromDocument(cfg)

        # first encoding space
        self.enc = self.header.encoding[0]

        #Parallel imaging factor
        self.acc_factor = self.enc.parallelImaging.accelerationFactor.kspace_encoding_step_1

        self.slc = self.enc.encodingLimits.slice.maximum+1
        self.phs = self.enc.encodingLimits.phase.maximum+1
        self.phs_retro = int(self.params["phase"])

        self.array_data = IsmrmrdImageArray()

        print("ImageArrayRecon, maximal number of slice ", self.slc)

    def process(self, array_data):

        ndim = array_data.data.shape

        RO = ndim[0]
        E1 = ndim[1]
        E2 = ndim[2]
        CHA = ndim[3]
        PHS = ndim[4]
        S = ndim[5]
        SLC = ndim[6]

        print("ImageArrayRecon, receiving image array, ", ndim)

        has_suffcient_image=False

        if self.slc == 1 and self.phs_retro==PHS:
            print("ImageArrayRecon, maximal 1 slice is expected ... ")
            self.array_data=array_data
            self.put_next(array_data)
            has_sufficient_image=True
        elif SLC<self.slc:
            print("ImageArrayRecon, pass this array down the chain ... ")
            self.put_next(array_data)
            has_sufficient_image=False
        else:
            print("ImageArrayRecon, sufficient images are received ... ")
            self.array_data=array_data
            self.put_next(array_data)
            has_sufficient_image=True

        if has_sufficient_image is False:
            return 0

        # enough images received
        # convert Meta
        mt = list()
        for x in self.array_data.meta:
            curr_meta = ismrmrd.Meta.deserialize(x)
            mt.append(curr_meta)

        print("ImageArrayRecon, convert %d meta containers ... ", len(mt))
        print("ImageArrayRecon, receive %d waveforms ... ", len(self.array_data.waveform))

        if len(self.array_data.waveform)>0:
            print(self.array_data.waveform[0].version)

        for slc in range(0,SLC-1):
            for phs in range(0,PHS-1):
                h = np.ravel(self.array_data.headers)[phs, 0, slc]
                print(h)
                h.image_series_index += 100

        print(mt[0])
        print(self.array_data.headers[0, 0, 0])

        return 0
