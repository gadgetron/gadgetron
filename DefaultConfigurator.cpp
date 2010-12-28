#include "DefaultConfigurator.h"
#include "Gadget.h"
#include "AcquisitionPassthroughGadget.h"
#include "AcquisitionFinishGadget.h"
#include "ImageFinishGadget.h"
#include "AccumulatorGadget.h"
#include "FFTGadget.h"
#include "CropAndCombineGadget.h"
#include "ImageWriterGadget.h"

DefaultConfigurator::DefaultConfigurator(char* config, ACE_UINT16 config_len,GadgetStreamController* controller)
  : GadgetStreamConfigurator(config,config_len,controller)
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
  
  ACE_Module<ACE_MT_SYNCH> *accum = 0;
  ACE_NEW_RETURN (accum,
		  ACE_Module<ACE_MT_SYNCH> (ACE_TEXT ("Accu"),
					    new AccumulatorGadget ()),
		  -1);

  ACE_Module<ACE_MT_SYNCH> *fft = 0;
  ACE_NEW_RETURN (fft,
		  ACE_Module<ACE_MT_SYNCH> (ACE_TEXT ("FFT"),
					    new FFTGadget ()),
		  -1);

  ACE_Module<ACE_MT_SYNCH> *cropcombine = 0;
  ACE_NEW_RETURN (cropcombine,
		  ACE_Module<ACE_MT_SYNCH> (ACE_TEXT ("CropCombine"),
					    new CropAndCombineGadget ()),
		  -1);

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


  if (stream->push (cropcombine) == -1)
    ACE_ERROR_RETURN ((LM_ERROR,
		       ACE_TEXT ("Failed to push %p\n"),
		       ACE_TEXT ("CropCombine")),
		      -1);

  if (stream->push (fft) == -1)
    ACE_ERROR_RETURN ((LM_ERROR,
		       ACE_TEXT ("Failed to push %p\n"),
		       ACE_TEXT ("FFT")),
		      -1);

  /*
  if (stream->push (imwriter) == -1)
    ACE_ERROR_RETURN ((LM_ERROR,
		       ACE_TEXT ("Failed to push %p\n"),
		       ACE_TEXT ("ImageWriter")),
		      -1);
  */

  if (stream->push (accum) == -1)
    ACE_ERROR_RETURN ((LM_ERROR,
		       ACE_TEXT ("Failed to push %p\n"),
		       ACE_TEXT ("Accumulator")),
		      -1);

  if (stream->push (passThrough) == -1)
    ACE_ERROR_RETURN ((LM_ERROR,
		       ACE_TEXT ("Failed to push %p\n"),
		       ACE_TEXT ("PassThrough")),
		      -1);

  
  return 0;
}
