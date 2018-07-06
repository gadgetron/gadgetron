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

        print("\nImageArrayRecon, receiving image array, ", ndim)
        
        # self.put_next(array_data)

        print("ImageArrayRecon, parse meta ... ")
        mt = list()
        for x in array_data.meta:
            curr_meta = ismrmrd.Meta.deserialize(x)
            mt.append(curr_meta)

        print("ImageArrayRecon, convert %d meta containers ... ", len(mt))

        if array_data.waveform is not None:
            print("ImageArrayRecon, receive %d waveforms ... ", len(array_data.waveform))

            if len(array_data.waveform)>0:
                print(array_data.waveform[0].version)
            
        for slc in range(0,SLC):
            for s in range(0,S):
                for phs in range(0,PHS):
                    print("send out image %d-%d-%d" % (phs, s, slc))
                    a = array_data.data[:,:,:,:,phs,s,slc]
                    print(a.shape)
                    self.put_next(array_data.headers[phs,s,slc], a)
        
        """
        print("ImageArrayRecon, receiving image array, ", ndim)

        print("ImageArrayRecon, buffer incoming array ... ")
        self.array_data=array_data

        print("ImageArrayRecon, pass down incoming array ... ") 
        self.put_next(array_data)

        print("ImageArrayRecon, parse meta ... ")
        mt = list()
        for x in self.array_data.meta:
            curr_meta = ismrmrd.Meta.deserialize(x)
            mt.append(curr_meta)

        print("ImageArrayRecon, convert %d meta containers ... ", len(mt))
        print("ImageArrayRecon, receive %d waveforms ... ", len(self.array_data.waveform))

        if len(self.array_data.waveform)>0:
            print(self.array_data.waveform[0].version)

        for slc in range(0,SLC-1):
            h = np.ravel(self.array_data.headers)[0, 0, slc]
            print(h)
            h.image_series_index += 100

        if len(mt)>0:
            print(mt[0])
            
        print(self.array_data.headers[0, 0, 0])
        """
        return 0
