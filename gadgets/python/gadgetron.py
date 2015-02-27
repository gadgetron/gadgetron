try:
    import GadgetronPythonMRI
except ImportError:
    pass


class Gadget(object):
    def __init__(self, gadget_reference, next_gadget=None):
        self._gadget_reference = gadget_reference
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
            self.next_gadget.process(args)
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
    def __init__(self, fn, gadget_reference, next_gadget=None):
        super(FunctionGadget, self).__init__(gadget_reference, next_gadget)
        self.process = fn
