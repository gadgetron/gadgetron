#include "GPUCGConfigurator.h"

#include "Gadgetron.h"
#include "GPUCGGoldenRadial.h"
#include "GPUCGGoldenSpiral.h"
#include "ImageFinishGadget.h"
#include "ImageWriterGadget.h"

GPUCGConfigurator::GPUCGConfigurator(char* config, ACE_UINT16 config_len,GadgetStreamController* controller)
  : GadgetStreamConfigurator(config,config_len,controller)
{

  

}

int GPUCGConfigurator::ConfigureStream(ACE_Stream<ACE_MT_SYNCH>* stream)
{

  ConfigParser cp;
  cp.parse(config_);

  int slices = cp.getIntVal(std::string("encoding"), std::string("slices"));
  std::string trajectory = cp.getStringVal(std::string("encoding"), 
					   std::string("trajectory"));

  

  ACE_Module<ACE_MT_SYNCH> *head = 0;
  ACE_Module<ACE_MT_SYNCH> *tail = 0;

  if (tail == 0) {
    ACE_NEW_RETURN(tail, 
		   ACE_Module<ACE_MT_SYNCH>( ACE_TEXT("EndGadget"), 
					     new EndGadget() ),
		   -1);
    stream->open(0,head,tail);
  }

  

  /*
  ACE_Module<ACE_MT_SYNCH> *gpucg0 = 0;
  ACE_NEW_RETURN (gpucg0,
		  ACE_Module<ACE_MT_SYNCH> (ACE_TEXT ("GPUCG0"),
					    new GPUCGGoldenRadialGadget (true, 0)),
		  -1);

  ACE_Module<ACE_MT_SYNCH> *gpucg1 = 0;
  ACE_NEW_RETURN (gpucg1,
		  ACE_Module<ACE_MT_SYNCH> (ACE_TEXT ("GPUCG1"),
					    new GPUCGGoldenRadialGadget (true, 1)),
		  -1);

  ACE_Module<ACE_MT_SYNCH> *gpucg2 = 0;
  ACE_NEW_RETURN (gpucg2,
		  ACE_Module<ACE_MT_SYNCH> (ACE_TEXT ("GPUCG2"),
					    new GPUCGGoldenRadialGadget (true, 2)),
		  -1);
  */

  /*
  ACE_Module<ACE_MT_SYNCH> *imwriter = 0;
  ACE_NEW_RETURN (imwriter,
		  ACE_Module<ACE_MT_SYNCH> (ACE_TEXT ("ImageWriter"),
					    new ImageWriterGadget ()),
		  -1);
  
  */

  ACE_Module<ACE_MT_SYNCH> *imaFinish = 0;
  ACE_NEW_RETURN (imaFinish,
		  ACE_Module<ACE_MT_SYNCH> (ACE_TEXT ("ImageFinish"),
			  new ImageFinishGadget (controller_)),
		  -1);

  if (stream->push (imaFinish) == -1)
    ACE_ERROR_RETURN ((LM_ERROR,
		       ACE_TEXT ("Failed to push %p\n"),
		       ACE_TEXT ("ImageFinish")),
		      -1);

  /*
  if (stream->push (imwriter) == -1)
    ACE_ERROR_RETURN ((LM_ERROR,
		       ACE_TEXT ("Failed to push %p\n"),
		       ACE_TEXT ("ImageWriter")),
		      -1);

  */

  for (int s = slices; s > 0; s--) {
      ACE_Module<ACE_MT_SYNCH> *gpucg = 0;
      if (trajectory.compare(std::string("spiral")) == 0) {
	ACE_NEW_RETURN (gpucg,
			ACE_Module<ACE_MT_SYNCH> (ACE_TEXT ("GPUCG0"),
						  new GPUCGGoldenSpiralGadget (true, s-1)),
			GADGET_FAIL);

      } else if (trajectory.compare(std::string("radial")) == 0) {
	ACE_NEW_RETURN (gpucg,
			ACE_Module<ACE_MT_SYNCH> (ACE_TEXT ("GPUCG0"),
						  new GPUCGGoldenRadialGadget (true, s-1)),
			GADGET_FAIL);
	
      } else {
	GADGET_DEBUG2("Unknown trajectory: %s", trajectory.c_str());
	return GADGET_FAIL;
      }

      if (stream->push (gpucg) == -1) {
	GADGET_DEBUG1("Failed to push GPUCG Gadget\n");
	return GADGET_FAIL;	
      }

  }

  /*

  if (stream->push (gpucg2) == -1) {
    GADGET_DEBUG1("Failed to push GPUCG Gadget\n");
    return -1;
    
  }

  if (stream->push (gpucg1) == -1) {
    GADGET_DEBUG1("Failed to push GPUCG Gadget\n");
    return -1;
    
  }

  if (stream->push (gpucg0) == -1) {
    GADGET_DEBUG1("Failed to push GPUCG Gadget\n");
    return -1;
    
  }
  */
  
  return 0;
}
