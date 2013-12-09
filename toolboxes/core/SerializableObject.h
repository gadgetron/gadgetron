/********************************************************************
    created:    2013/09/30
    created:    30:9:2013   16:55
    author:     Hui Xue

    purpose:    provide the interface for serializable objects
*********************************************************************/
#pragma once

#include <complex>
#include <iostream>
#include "cpucore_export.h"
#include "GadgetronCommon.h"
#include "GadgetronException.h"

namespace Gadgetron
{

class SerializableObject
{
public:

    SerializableObject() {}
    virtual ~SerializableObject() {}

    // serialize and deserialize to/from the buffer
    virtual bool serialize(char*& buf, unsigned long long& len) const = 0;
    virtual bool deserialize(char* buf, unsigned long long& len) = 0;
};

}
