from __future__ import print_function
import numpy as np
from gadgetron import Gadget

class ArrayImage(Gadget):

    def process_config(self, cfg):
        print("RMS Coil Combine, Config ignored")
        self.images = []
        self.headers = []
        self.metas = []


    def process(self, header, image, metadata=None):
        self.images.append(image)
        self.headers.append(header)
        if metadata is not None:
            self.metas.append(metadata)



        if (header.flags & (1<<7)): #Is this the last scan in slice
            combined_image = np.stack(self.iamges,axis=-1)

            #Create segmentation here.


            print("RMS coil",image.shape,combined_image.shape)

            #Send the combined image and the last header and metadata
            if metadata is not None:
                self.put_next(header,combined_image,metadata)
            else:
                self.put_next(header,combined_image,metadata)
            self.images = []
            self.headers = []
            self.metas = []
        return 0