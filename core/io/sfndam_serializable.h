#pragma once

#include "sfndam.h"

namespace Gadgetron::Core::IO {

template <class T> struct SfndamSerializable {
    virtual ~SfndamSerializable() {}
    virtual void SerializeToSfndam(std::ostream& stream) const = 0;
    static T DeserializeSfndam(std::istream& stream);
};

template <class T>
std::enable_if_t<std::is_base_of_v<Gadgetron::Core::IO::SfndamSerializable<T>, T>> write(std::ostream &stream, const T &t);

} // namespace Gadgetron::Core::IO
