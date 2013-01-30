#include <octave/oct.h>
#include <octave/ov-struct.h>

#include "OctaveCommunicator.h"
     
DEFUN_DLD (GadgetronReturnIsmrmrdImage, args, nargout,
	   "GadgetronReturnIsmrmrdImage return Image to the Gadgetron")
{
  int nargin = args.length ();

  octave_value retval;
     
  if (nargin != 1) {
    print_usage(); 
  } else {
    std::string id(args(0).string_value());

    ACE_Message_Block* m = 0;
    OctaveCommunicator::instance()->message_gadget(id, m);
  }
  return octave_value_list ();
}
