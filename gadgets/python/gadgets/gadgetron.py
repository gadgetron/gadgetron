try:
    import GadgetronPythonMRI
    import ismrmrd
except ImportError:
    pass

import time
import numpy as np

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
        self.put_next(header,*args)

    def wait(self):
        pass
    
    def put_next(self, *args):
        if self.next_gadget is not None:
            if isinstance(self.next_gadget, Gadget):
                if len(args) == 3 and not isinstance(args[2],ismrmrd.Meta): #Data with meta data we assume
                    meta = ismrmrd.Meta()
                    meta = ismrmrd.Meta.deserialize(args[2])
                    new_args = (args[0], args[1], meta)
                    self.next_gadget.process(*new_args)
                else:
                    self.next_gadget.process(*args)
            elif isinstance(self.next_gadget, GadgetronPythonMRI.GadgetReference):
                if len(args) > 3:
                    raise Exception("Only two or 3 return arguments are currently supported when returning to Gadgetron framework")
                if isinstance(args[0], ismrmrd.AcquisitionHeader):
                    self.next_gadget.return_acquisition(args[0],args[1].astype('complex64'))
                elif isinstance(args[0], ismrmrd.ImageHeader):
                    header = args[0]
                    if (args[1].dtype == np.uint16):
                        if len(args) == 3:
                            self.next_gadget.return_image_ushort_attr(header,args[1], args[2].serialize())
                        else:
                            self.next_gadget.return_image_ushort(header,args[1])
                    elif (args[1].dtype == np.float32):
                        if len(args) == 3:
                            self.next_gadget.return_image_float_attr(header, args[1], args[2].serialize())
                        else:
                            self.next_gadget.return_image_float(header,args[1])
                    else:
                        if len(args) == 3:
                            self.next_gadget.return_image_cplx_attr(header, args[1].astype('complex64'), args[2].serialize())
                        else:
                            self.next_gadget.return_image_cplx(header,args[1].astype('complex64'))
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
        
class WrapperGadget(Gadget):
    
    def __init__(self, dllname, classname, gadgetname=None, next_gadget=None):
        if gadgetname == None:
            gadgetname = classname
        Gadget.__init__(self, next_gadget)
        self.controller_ = GadgetronPythonMRI.GadgetInstrumentationStreamController()
        self.controller_.prepend_gadget(gadgetname,dllname,classname)
        self.controller_.set_python_gadget(self)
    
    def prepend_gadget(self,dllname, classname, gadgetname=None):
        self.controller_.prepend_gadget(gadgetname,dllname,classname)

    def wait(self):
       self.controller_.close()

    def process_config(self, conf):
        self.controller_.put_config(conf)
        return 0

    def process(self, header, *args):
        if len(args) > 2:
            raise("Only two or three arguments are currently supported when sending data to Gadgetron framework")
        if isinstance(header, ismrmrd.AcquisitionHeader):
            self.controller_.put_acquisition(header,args[0].astype('complex64'))
        elif isinstance(header, ismrmrd.ImageHeader):
            if (args[0].dtype == np.uint16):
                if len(args) == 2:
                    self.controller_.put_image_ushort_attr(header,args[0], args[1].serialize())
                else:
                    self.controller_.put_image_ushort(header,args[0])
            elif (args[0].dtype == np.float32):
                if len(args) == 2:
                    self.controller_.put_image_float_attr(header, args[0], args[1].serialize())
                else:
                    self.controller_.put_image_float(header,args[0])
            else:   
                if len(args) == 2:
                    self.controller_.put_image_cplx_attr(header, args[0].astype('complex64'), args[1].serialize())
                else:
                    self.controller_.put_image_cplx(header,args[0].astype('complex64'))
        else:
            raise("Unsupported types when sending data to Gadgetron framework")
        return 0

    def set_parameter(self,gadgetname,parameter,value):
        self.controller_.set_parameter(gadgetname,parameter,value)

class FunctionGadget(Gadget):
    """A Gadget with a configurable `process` function.

    Params:
        fn: `process` function
    """
    def __init__(self, fn, next_gadget=None):
        super(FunctionGadget, self).__init__(next_gadget)
        self.process = fn


def gadget_chain_wait(first_gadget):
    g = first_gadget;
    while (g):
        g.wait()
        g = g.next_gadget

def gadget_chain_config(first_gadget, conf):
    g = first_gadget;
    while (g):
        g.process_config(conf)
        g = g.next_gadget
        
def get_last_gadget(first_gadget):
    g = first_gadget;
    while (True):
        if g.next_gadget:
            g = g.next_gadget
        else:
            break
    return g
