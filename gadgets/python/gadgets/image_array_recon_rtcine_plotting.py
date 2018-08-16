import os
import sys
import pickle
import ismrmrd
import ismrmrd.xsd
import numpy as np
import platform
# from ismrmrdtools import transform, coils, grappa
from gadgetron import Gadget,IsmrmrdImageArray
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

        curr_slc = array_data.headers.slice
        print("\ImageArrayReconRTCinePlotting, receiving image array, ", ndim, " for slice ", curr_slc)

        print("ImageArrayReconRTCinePlotting, parse meta ... ")
        mt = list()
        for x in array_data.meta:
            curr_meta = ismrmrd.Meta.deserialize(x)
            mt.append(curr_meta)

        print("ImageArrayReconRTCinePlotting, convert %d meta containers ... " % len(mt))

        if array_data.waveform is not None:
            print("ImageArrayReconRTCinePlotting, receive %d waveforms ... " % len(array_data.waveform))
            print("ImageArrayReconRTCinePlotting, receive waveforms ", array_data.waveform)

        if array_data.acq_headers is not None:
            print("ImageArrayReconRTCinePlotting, receive acq headers ... ", array_data.acq_headers.shape)

        if (self.num_processed_==0):
            # allocate data buffer
            self.data = np.zeros([RO, E1, E2, CHA, PHS, S, self.slc])
            self.header = np.empty([PHS, S, self.slc], dtype=ismrmrd.ImageHeader)
            self.acq_headers = np.empty([E2, 1, PHS, S, self.slc], dtype=ismrmrd.AcquisitionHeader)

        self.data[:,:,:,:,:,:,curr_slc] = array_data.data
        self.headers[:,:,curr_slc] = array_data.headers
        self.waveforms.append(array_data.waveform)
        self.acq_headers[:,:, :,:,curr_slc] = array_data.acq_headers
        self.meta.append(array_data.meta)

        if (self.debug_folder is not None):
            save_file = os.path.join(self.debug_folder, "image_array"+str(self.num_processed_)+".dat")
            with open(save_file, "wb") as f:
                pickle.dump(array_data, f, pickle.HIGHEST_PROTOCOL)
                print("Save incoming array data to %s" % save_file)

        # plot the ecg

        # ecg_plot = plot_ecg(array_data.acq_headers, array_data.waveform)

        # make the image header for ecg_plot

        # set the meta fields for ecg_plot

        # send out images to next gadget
        for slc in range(0,SLC):
            for s in range(0,S):
                for phs in range(0,PHS):
                    print("send out image %d-%d-%d" % (phs, s, slc))
                    a = array_data.data[:,:,:,:,phs,s,slc]
                    # print(a.shape)
                    self.put_next(array_data.headers[phs,s,slc], a, array_data.meta[phs+s*PHS + slc*S*PHS])

        self.num_processed_+=1

        if(self.num_processed_==self.slc):
            print("ImageArrayReconRTCinePlotting, receive all slices ... ")

            # save data
            save_file = os.path.join(self.debug_folder, "image_array_all_data.dat")
            with open(save_file, "wb") as f:
                pickle.dump(self.data, f)
                print("Save incoming array data to %s" % save_file)

            save_file = os.path.join(self.debug_folder, "image_array_all_headers.dat")
            with open(save_file, "wb") as f:
                pickle.dump(self.headers, f)
                print("Save incoming array headers to %s" % save_file)

            save_file = os.path.join(self.debug_folder, "image_array_all_waveforms.dat")
            with open(save_file, "wb") as f:
                pickle.dump(self.waveforms, f)
                print("Save incoming array waveforms to %s" % save_file)

            save_file = os.path.join(self.debug_folder, "image_array_all_acq_headers.dat")
            with open(save_file, "wb") as f:
                pickle.dump(self.acq_headers, f)
                print("Save incoming array acq_headers to %s" % save_file)

            save_file = os.path.join(self.debug_folder, "image_array_all_acq_meta.dat")
            with open(save_file, "wb") as f:
                pickle.dump(self.meta, f)
                print("Save incoming array meta to %s" % save_file)

        return 0
