/** 
    SerializeObject is the base class for serializable objects
*/

#pragma once

#include "GadgetronCommon.h"
#include "GadgetronException.h"
#include "cpucore_export.h"

#include <complex>
#include <iostream>

namespace Gadgetron
{
  class SerializableObject
  {
  public:
    
    SerializableObject() {}
    virtual ~SerializableObject() {}
    
    // serialize and deserialize to/from the buffer
    virtual void serialize(char*& buf, size_t& len) const = 0;
    virtual void deserialize(char* buf, size_t& len) = 0;
  };  
}
