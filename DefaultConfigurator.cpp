#include "DefaultConfigurator.h"
#include "Gadget.h"
#include "AcquisitionPassthroughGadget.h"

DefaultConfigurator::DefaultConfigurator(char* config, ACE_UINT16 config_len)
  : GadgetStreamConfigurator(config,config_len)
{

  

}

int DefaultConfigurator::ConfigureStream(ACE_Stream<ACE_MT_SYNCH>* stream)
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


  ACE_Module<ACE_MT_SYNCH> *passThrough = 0;
  ACE_NEW_RETURN (passThrough,
		  ACE_Module<ACE_MT_SYNCH> (ACE_TEXT ("PassThrough"),
			  new AcquisitionPassthroughGadget ()),
		  -1);
  
  if (stream->push (passThrough) == -1)
    ACE_ERROR_RETURN ((LM_ERROR,
		       ACE_TEXT ("Failed to push %p\n"),
		       ACE_TEXT ("PassThrough")),
		      -1);

  
  return 0;
}
