#ifndef WEBINTERFACEGADGET_H
#define WEBINTERFACEGADGET_H

#include "gadgetronwebinterface_export.h"
#include "Gadget.h"
#include "GadgetStreamController.h"

#include <complex>

class EXPORTGADGETSWEBINTERFACE WebInterfaceGadget : public Gadget
{
 public:
  GADGET_DECLARE(WebInterfaceGadget);
  
  //We have to overwrite close in this gadget to make sure we shutdown the http server.
  virtual int close(unsigned long flags);

 protected:
  virtual int process(ACE_Message_Block * m);
  virtual int process_config(ACE_Message_Block * m);

  struct mg_context *mongoose_ctx_;
};


#endif //WEBINTERFACEGADGET_H
