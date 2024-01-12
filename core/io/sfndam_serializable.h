#pragma once

#include "sfndam.h"

namespace Gadgetron::Core::IO {

template <class T> struct SfndamSerializable {
    virtual ~SfndamSerializable() {}
    virtual void SerializeToSfndam(std::ostream& stream) const = 0;
    static T DeserializeSfndam(std::istream& stream);
};

} // namespace Gadgetron::Core::IO