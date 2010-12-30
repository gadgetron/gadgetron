#include "GPUCGConfigurator.h"

#include "Gadgetron.h"
#include "GPUCGGadget.h"
#include "ImageFinishGadget.h"

GPUCGConfigurator::GPUCGConfigurator(char* config, ACE_UINT16 config_len,GadgetStreamController* controller)
  : GadgetStreamConfigurator(config,config_len,controller)
{

  

}

int GPUCGConfigurator::ConfigureStream(ACE_Stream<ACE_MT_SYNCH>* stream)
{

  ACE_Module<ACE_MT_SYNCH> *head = 0;
  ACE_Module<ACE_MT_SYNCH> *tail = 0;

  if (tail == 0) {
    ACE_NEW_RETURN(tail, 
		   ACE_Module<ACE_MT_SYNCH>( ACE_TEXT("EndGadget"), 
					     new EndGadget() ),
		   -1);
    stream->open(0,head,tail);
  }


  ACE_Module<ACE_MT_SYNCH> *gpucg = 0;
  ACE_NEW_RETURN (gpucg,
		  ACE_Module<ACE_MT_SYNCH> (ACE_TEXT ("GPUCG"),
			  new GPUCGGadget ()),
		  -1);
  
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


  if (stream->push (gpucg) == -1) {
    GADGET_DEBUG1("Failed to push GPUCG Gadget\n");
    return -1;
    
  }
  
  return 0;
}
