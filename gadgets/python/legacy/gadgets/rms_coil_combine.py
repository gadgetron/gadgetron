import numpy as np
from gadgetron import Gadget


class RMSCoilCombine(Gadget):

    def process(self, h, im):
        combined_image = np.sqrt(np.sum(np.square(np.abs(im)), axis=len(im.shape) - 1))
        print("RMS coil", im.shape, combined_image.shape)
        h.channels = 1
        self.put_next(h, combined_image.astype('complex64'))
        return 0
