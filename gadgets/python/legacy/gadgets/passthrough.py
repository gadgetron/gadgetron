from gadgetron import Gadget
class Passthrough(Gadget):
    def __init__(self, next_gadget=None):
        Gadget.__init__(self,next_gadget)

    def process_config(self, conf):
        pass 
    def process(self, recondata,*args):
         #Return image to Gadgetron
        self.put_next(recondata)
        return 0    
        #print "Returning to Gadgetron"
 
