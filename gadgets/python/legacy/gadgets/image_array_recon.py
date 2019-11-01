import os
import sys
import pickle
import ismrmrd
import ismrmrd.xsd
import numpy as np
import platform

from gadgetron import Gadget, IsmrmrdImageArray
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

        try:
            self.debug_folder = self.params["debug_folder"]

            if (len(self.debug_folder)==0):
                if platform.system() == "Windows":
                    self.debug_folder = "C:/temp/gadgetron"
                else:
                    self.debug_folder = "/tmp/gadgetron"
        except:
            self.debug_folder = None

        self.num_processed_ = 0

        print("ImageArrayRecon, maximal number of slice ", self.slc)
        print("ImageArrayRecon, find debug folder ", self.debug_folder)

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

        if (self.debug_folder is not None):
            save_file = os.path.join(self.debug_folder, "image_array"+str(self.num_processed_)+".dat")
            with open(save_file, "wb") as f:
                pickle.dump(array_data, f, pickle.HIGHEST_PROTOCOL)
                print("Save incoming array data to %s" % save_file)

        # self.put_next(array_data)

        print("ImageArrayRecon, parse meta ... ")
        mt = list()
        for x in array_data.meta:
            curr_meta = ismrmrd.Meta.deserialize(x)
            mt.append(curr_meta)

        print("ImageArrayRecon, convert %d meta containers ... " % len(mt))

        if array_data.waveform is not None:
            print("ImageArrayRecon, receive %d waveforms ... " % len(array_data.waveform))
            print("ImageArrayRecon, receive waveforms ", array_data.waveform)

        if array_data.acq_headers is not None:
            print("ImageArrayRecon, receive acq headers ... ", array_data.acq_headers.shape)

        for slc in range(0,SLC):
            for s in range(0,S):
                for phs in range(0,PHS):
                    print("send out image %d-%d-%d" % (phs, s, slc))
                    a = array_data.data[:,:,:,:,phs,s,slc]
                    # print(a.shape)
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
