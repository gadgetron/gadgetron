import os
import sys
import pickle
import platform

import ismrmrd
import ismrmrd.xsd
import numpy as np

from gadgetron import Gadget, IsmrmrdImageArray

class ImageArrayReconRTCinePlotting(Gadget):

    def process_config(self, cfg):
        print("ImageArrayReconRTCinePlotting, process config ... ", file=sys.stderr)
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

        print(f"ImageArrayReconRTCinePlotting, maximal number of slice {self.slc}", file=sys.stderr)
        print(f"ImageArrayReconRTCinePlotting, find debug folder {self.debug_folder}", file=sys.stderr)

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
        print(f"\ImageArrayReconRTCinePlotting, receiving image array, {ndim} for slice {curr_slc}", file=sys.stderr)

        print("ImageArrayReconRTCinePlotting, parse meta ... ", file=sys.stderr)
        mt = list()
        for x in array_data.meta:
            curr_meta = ismrmrd.Meta.deserialize(x)
            mt.append(curr_meta)

        print(f"ImageArrayReconRTCinePlotting, convert {len(mt)} meta containers ... ", file=sys.stderr)

        buffer_for_plotting = True
        if array_data.waveform is not None:
            print(f"ImageArrayReconRTCinePlotting, receive {len(array_data.waveform)} waveforms ... ", file=sys.stderr)
        else:
            buffer_for_plotting = False

        if array_data.acq_headers is not None:
            print(f"ImageArrayReconRTCinePlotting, receive acq headers ... {array_data.acq_headers.shape}", file=sys.stderr)
        else:
            buffer_for_plotting = False
            print("Will not buffer for plotting ... ", file=sys.stderr)

        if buffer_for_plotting:
            if (self.num_processed_==0):
                # allocate data buffer
                self.data = np.zeros([RO, E1, E2, CHA, PHS, S, self.slc])
                print(f"Create data buffer {self.data.shape}", file=sys.stderr)
                self.headers = array_data.headers
                self.acq_headers = array_data.acq_headers
                print(f"Create headers buffer {self.headers.shape}", file=sys.stderr)
                print(f"Create acq_headers buffer {self.acq_headers.shape}", file=sys.stderr)

            self.data[:,:,:,:,:,:,curr_slc] = array_data.data.reshape([RO, E1, E2, CHA, PHS, S])
            self.waveforms.extend(array_data.waveform)
            self.meta.extend(array_data.meta)

            if (self.num_processed_>0):
                print(f"Incoming headers  {array_data.headers.shape}", file=sys.stderr)
                self.headers = np.append(self.headers, array_data.headers, axis=2)
                print(f"Incoming acq headers  {array_data.acq_headers.shape}", file=sys.stderr)
                self.acq_headers = np.append(self.acq_headers, array_data.acq_headers, axis=3)

            if (self.debug_folder is not None):
                save_file = os.path.join(self.debug_folder, "image_array"+str(self.num_processed_)+".dat")
                with open(save_file, "wb") as f:
                    pickle.dump(array_data, f, pickle.HIGHEST_PROTOCOL)
                    print(f"Save incoming array data to {save_file}", file=sys.stderr)

        # send out images to next gadget
        for slc in range(0,SLC):
            for s in range(0,S):
                for phs in range(0,PHS):
                    print("send out image %d-%d-%d" % (phs, s, slc), file=sys.stderr)
                    a = array_data.data[:,:,:,:,phs,s,slc]
                    print(a.shape, file=sys.stderr)
                    self.put_next(array_data.headers[phs,s,slc], a, array_data.meta[phs+s*PHS + slc*S*PHS])

        if buffer_for_plotting:
            self.num_processed_+=1

            if(self.num_processed_==self.slc):
                print("ImageArrayReconRTCinePlotting, receive all slices ... ", file=sys.stderr)

                # save data
                save_file = os.path.join(self.debug_folder, "image_array_all_data.dat")
                with open(save_file, "wb") as f:
                    pickle.dump(self.data, f)
                    print(f"Save incoming array data to {save_file}", file=sys.stderr)
                    f.close()

                save_file = os.path.join(self.debug_folder, "image_array_all_headers.dat")
                with open(save_file, "wb") as f:
                    pickle.dump(self.headers, f)
                    print(f"Save incoming array headers to {save_file}", file=sys.stderr)
                    f.close()

                save_file = os.path.join(self.debug_folder, "image_array_all_waveforms.dat")
                with open(save_file, "wb") as f:
                    pickle.dump(self.waveforms, f)
                    print(f"Save incoming array waveforms to {save_file}", file=sys.stderr)
                    f.close()

                save_file = os.path.join(self.debug_folder, "image_array_all_acq_headers.dat")
                with open(save_file, "wb") as f:
                    pickle.dump(self.acq_headers, f)
                    print(f"Save incoming array acq_headers to {save_file}", file=sys.stderr)
                    f.close()
                save_file = os.path.join(self.debug_folder, "image_array_all_acq_meta.dat")
                with open(save_file, "wb") as f:
                    pickle.dump(self.meta, f)
                    print(f"Save incoming array meta to {save_file}", file=sys.stderr)
                    f.close()

        return 0
