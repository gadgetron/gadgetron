
from gadgetron import Gadget
from gadgetron.util.cfft import cfftn, cifftn


class Remove2xOversampling(Gadget):

    def process(self, acq, data):
        orig_size = list(data.shape)
        data2 = data.reshape([data.shape[0], int(data.size/data.shape[0])])
        new_length = data2.shape[0] // 2
        data2 = cfftn(cifftn(data2, axes=[0])[(0+(new_length // 2)):(new_length+(new_length // 2)), :], axes=[0])
        orig_size[0] = new_length
        data2.reshape(tuple(orig_size))
        acq.samples = new_length

        self.put_next(acq, data2)
        return 0
