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

class ImageArrayReconRTCinePlotting(Gadget):

    def process_config(self, cfg):
        print("ImageArrayReconRTCinePlotting, process config ... ")
        self.header = ismrmrd.xsd.CreateFromDocument(cfg)

        # first encoding space
        self.enc = self.header.encoding[0]

        #Parallel imaging factor
        self.acc_factor = self.enc.parallelImaging.accelerationFactor.kspace_encoding_step_1
        self.slc = self.enc.encodingLimits.slice.maximum+1

        self.data = None
        self.headers = None
        self.waveforms = list()
        self.acq_headers = None
        self.meta = list()

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

        print("ImageArrayReconRTCinePlotting, maximal number of slice ", self.slc)
        print("ImageArrayReconRTCinePlotting, find debug folder ", self.debug_folder)

    def process(self, array_data):

        ndim = array_data.data.shape

        RO = ndim[0]
        E1 = ndim[1]
        E2 = ndim[2]
        CHA = ndim[3]
        PHS = ndim[4]
        S = ndim[5]
        SLC = ndim[6]

        curr_slc = array_data.headers[0,0,0].slice
        print("\ImageArrayReconRTCinePlotting, receiving image array, ", ndim, " for slice ", curr_slc)

        print("ImageArrayReconRTCinePlotting, parse meta ... ")
        mt = list()
        for x in array_data.meta:
            curr_meta = ismrmrd.Meta.deserialize(x)
            mt.append(curr_meta)

        print("ImageArrayReconRTCinePlotting, convert %d meta containers ... " % len(mt))

        buffer_for_plotting = True
        if array_data.waveform is not None:
            print("ImageArrayReconRTCinePlotting, receive %d waveforms ... " % len(array_data.waveform))
        else:
            buffer_for_plotting = False

        if array_data.acq_headers is not None:
            print("ImageArrayReconRTCinePlotting, receive acq headers ... ", array_data.acq_headers.shape)
        else:
            buffer_for_plotting = False
            print("Will not buffer for plotting ... ")

        if buffer_for_plotting:
            if (self.num_processed_==0):
                # allocate data buffer
                self.data = np.zeros([RO, E1, E2, CHA, PHS, S, self.slc])
                print("Create data buffer", self.data.shape)
                self.headers = array_data.headers
                self.acq_headers = array_data.acq_headers
                print("Create headers buffer", self.headers.shape)
                print("Create acq_headers buffer", self.acq_headers.shape)

            self.data[:,:,:,:,:,:,curr_slc] = array_data.data.reshape([RO, E1, E2, CHA, PHS, S])
            self.waveforms.extend(array_data.waveform)
            self.meta.extend(array_data.meta)

            if (self.num_processed_>0):
                print("Incoming headers ", array_data.headers.shape)
                self.headers = np.append(self.headers, array_data.headers, axis=2)
                print("Incoming acq headers ", array_data.acq_headers.shape)
                self.acq_headers = np.append(self.acq_headers, array_data.acq_headers, axis=3)

            if (self.debug_folder is not None):
                save_file = os.path.join(self.debug_folder, "image_array"+str(self.num_processed_)+".dat")
                with open(save_file, "wb") as f:
                    pickle.dump(array_data, f, pickle.HIGHEST_PROTOCOL)
                    print("Save incoming array data to %s" % save_file)

        # send out images to next gadget
        for slc in range(0,SLC):
            for s in range(0,S):
                for phs in range(0,PHS):
                    # print("send out image %d-%d-%d" % (phs, s, slc))
                    a = array_data.data[:,:,:,:,phs,s,slc]
                    # print(a.shape)
                    self.put_next(array_data.headers[phs,s,slc], a, array_data.meta[phs+s*PHS + slc*S*PHS])

        if buffer_for_plotting:
            self.num_processed_+=1

            if(self.num_processed_==self.slc):
                print("ImageArrayReconRTCinePlotting, receive all slices ... ")

                # save data
                save_file = os.path.join(self.debug_folder, "image_array_all_data.dat")
                with open(save_file, "wb") as f:
                    pickle.dump(self.data, f)
                    print("Save incoming array data to %s" % save_file)
                    f.close()

                save_file = os.path.join(self.debug_folder, "image_array_all_headers.dat")
                with open(save_file, "wb") as f:
                    pickle.dump(self.headers, f)
                    print("Save incoming array headers to %s" % save_file)
                    f.close()

                save_file = os.path.join(self.debug_folder, "image_array_all_waveforms.dat")
                with open(save_file, "wb") as f:
                    pickle.dump(self.waveforms, f)
                    print("Save incoming array waveforms to %s" % save_file)
                    f.close()

                save_file = os.path.join(self.debug_folder, "image_array_all_acq_headers.dat")
                with open(save_file, "wb") as f:
                    pickle.dump(self.acq_headers, f)
                    print("Save incoming array acq_headers to %s" % save_file)
                    f.close()
                save_file = os.path.join(self.debug_folder, "image_array_all_acq_meta.dat")
                with open(save_file, "wb") as f:
                    pickle.dump(self.meta, f)
                    print("Save incoming array meta to %s" % save_file)
                    f.close()

        return 0
