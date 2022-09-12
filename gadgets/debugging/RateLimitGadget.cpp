#include "RateLimitGadget.h"
#include "ismrmrd/xml.h"
#include <thread>
#include <chrono>
namespace Gadgetron{
RateLimitGadget::RateLimitGadget()
{

}
 
RateLimitGadget::~RateLimitGadget()
{
}

/**
 *   Expects ISMRMRD XML configuration
 *
 */
int RateLimitGadget::process_config(ACE_Message_Block* mb)
{
  this->sleep_time = std::chrono::milliseconds(this->sleep_time_.value());
  return GADGET_OK;
}

int RateLimitGadget::process(ACE_Message_Block* mb)
{

  std::this_thread::sleep_for(this->sleep_time);
  this->next()->putq(mb);

  return GADGET_OK;

}

GADGET_FACTORY_DECLARE(RateLimitGadget)
}
