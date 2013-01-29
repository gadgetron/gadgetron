#include <octave/oct.h>

#include "OctaveCommunicator.h"
     
DEFUN_DLD (GadgetronReturnIsmrmrd, args, nargout,
	   "GadgetronReturnIsmrmrd will eventually return ISMRMRD data to the gadgetron")
{
  int nargin = args.length ();

  octave_value retval;
     
  if (nargin != 1) {
    print_usage(); 
  } else {
    std::string id(args(0).string_value());
    OctaveCommunicator::instance()->message_gadget(id);    
  }
  return octave_value_list ();
}
