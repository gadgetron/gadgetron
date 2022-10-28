import os
import pickle
import platform
import sys

import ismrmrd
import ismrmrd.xsd

from gadgetron import Gadget, IsmrmrdImageArray

class ImageArrayRecon(Gadget):

    def process_config(self, cfg):
        print("ImageArrayRecon, process config ... ", file=sys.stderr)
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

        print(f"ImageArrayRecon, maximal number of slice {self.slc}", file=sys.stderr)
        print(f"ImageArrayRecon, find debug folder {self.debug_folder}", file=sys.stderr)

    def process(self, array_data):

        ndim = array_data.data.shape

        RO = ndim[0]
        E1 = ndim[1]
        E2 = ndim[2]
        CHA = ndim[3]
        PHS = ndim[4]
        S = ndim[5]
        SLC = ndim[6]

        print(f"\nImageArrayRecon, receiving image array, {ndim}", file=sys.stderr)

        if (self.debug_folder is not None):
            save_file = os.path.join(self.debug_folder, "image_array"+str(self.num_processed_)+".dat")
            with open(save_file, "wb") as f:
                pickle.dump(array_data, f, pickle.HIGHEST_PROTOCOL)
                print(f"Save incoming array data to {save_file}", file=sys.stderr)

        # self.put_next(array_data)

        print("ImageArrayRecon, parse meta ... ", file=sys.stderr)
        mt = list()
        for x in array_data.meta:
            curr_meta = ismrmrd.Meta.deserialize(x)
            mt.append(curr_meta)

        print(f"ImageArrayRecon, convert {len(mt)} meta containers ... ", file=sys.stderr)

        if array_data.waveform is not None:
            print(f"ImageArrayRecon, receive %d waveforms ... {len(array_data.waveform)}", file=sys.stderr)
            print(f"ImageArrayRecon, receive waveforms {array_data.waveform}", file=sys.stderr)

        if array_data.acq_headers is not None:
            print(f"ImageArrayRecon, receive acq headers ... {array_data.acq_headers.shape}", file=sys.stderr)

        for slc in range(0,SLC):
            for s in range(0,S):
                for phs in range(0,PHS):
                    print(f"send out image {phs}-{s}-{slc}", file=sys.stderr)
                    a = array_data.data[:,:,:,:,phs,s,slc]
                    self.put_next(array_data.headers[phs,s,slc], a)
        
        return 0
