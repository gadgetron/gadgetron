from __future__ import print_function
import numpy as np
from gadgetron import Gadget

class ArrayImage(Gadget):
    def __init__(self,next_gadget=None):
        super(ArrayImage,self).__init__(next_gadget=next_gadget)
        self.counter = []

    def process_config(self, cfg):
        self.images = []
        self.headers = []
        self.metas = []
        self.counter = 0

    def process(self, header, image, metadata=None):


        self.images.append(image)
        self.headers.append(header)
        if metadata is not None:
            self.metas.append(metadata)
            print(metadata)

        self.counter += 1
        if self.counter == 10:
            print("Last scan in slice")
            combined_image = np.stack(self.images,axis=2) #Stack images along the E2 direction


            header.matrix_size[2] = combined_image.shape[2] #Adjust header to new dimension
            #Create segmentation here.





            #Send the combined image and the modified header and metadata
            if metadata is not None:
                self.put_next(header,combined_image,metadata)
            else:
                self.put_next(header,combined_image)
            self.images = []
            self.headers = []
            self.metas = []
            self.counter = 0
        return 0