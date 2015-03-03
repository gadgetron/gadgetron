try:
    import GadgetronPythonMRI
    import ismrmrd
except ImportError:
    pass


class Gadget(object):
    def __init__(self, next_gadget=None):
        self.next_gadget = next_gadget
        self.params = dict()
        self.results = []

    def set_parameter(self, name, value):
        self.params[name] = value

    def get_parameter(self, name):
        return self.params.get(name, None)

    def __call__(self, *args):
        self.process(args)
        return self.get_results()

    def set_next_gadget(self, gadget):
        self.next_gadget = gadget

    def process_config(self, conf):
        pass

    def process(self, header, *args):
        # do work here
        self.put_next(results)

    def put_next(self, *args):
        if self.next_gadget is not None:
            if isinstance(self.next_gadget, Gadget):
                self.next_gadget.process(args)
            elif isinstance(self.next_gadget, GadgetronPythonMRI.GadgetReference):
                if len(args) != 2:
                    raise("Only two return arguments are currently supported when returning to Gadgetron framework")
                if isinstance(args[0], ismrmrd.AcquisitionHeader):
                    self.next_gadget.return_acquisition(args[0],args[1].astype('complex64'))
                elif isinstance(args[0], ismrmrd.ImageHeader):
                    self.next_gadget.return_image(args[0],args[1].astype('complex64'))
                else:
                    raise("Unsupported types when returning to Gadgetron framework")
            else:
                raise("next_gadget is set to unsupported type")
        else:
            self.results.append(list(args))

    def get_results(self):
        results = self.results
        self.results = []
        return results

class FunctionGadget(Gadget):
    """A Gadget with a configurable `process` function.

    Params:
        fn: `process` function
    """
    def __init__(self, fn, next_gadget=None):
        super(FunctionGadget, self).__init__(next_gadget)
        self.process = fn
