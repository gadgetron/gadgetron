#pragma once

#include <istream>

namespace Gadgetron::Core::Readers {

    template<class T>
    void read_into(std::istream &stream, T &t) {
        stream.read(reinterpret_cast<char *>(&t), sizeof(t));
    }

    template<class T>
    T read_t(std::istream &stream) {
        T t;
        read_into(stream, t);
        return t;
    }

    template<class T>
    std::string read_string_from_stream(std::istream &stream) {

        auto n = read_t<T>(stream);

        std::string str(n, 0);
        stream.read(str.data(), n);

        return str;
    }
}
