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
    virtual bool serialize(char*& buf, size_t& len) const = 0; // Should be a void function
    virtual bool deserialize(char* buf, size_t& len) = 0; // Should be a void function
  };  
}
