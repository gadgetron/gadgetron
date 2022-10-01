import sys

from gadgetron import Gadget



class ArrayImagePassThrough(Gadget):
    def __init__(self,next_gadget=None):
        super(ArrayImagePassThrough,self).__init__(next_gadget=next_gadget)
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
            print(metadata, file=sys.stderr)

        #Send the combined image and the modified header and metadata
        if metadata is not None:
            self.put_next(header, image, metadata)
        else:
            self.put_next(header, image)

        return 0
