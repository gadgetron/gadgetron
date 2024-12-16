#pragma once

#include <set>
#include <vector>
#include <iostream>
#include <type_traits>

#include "hoNDArray.h"
#include "hoNDImage.h"
#include <boost/hana/adapt_struct.hpp>

namespace Gadgetron::Core::IO {

    template<class T>
    std::enable_if_t<std::is_trivially_copyable_v<T>> read(std::istream &stream, T &t);

    template<class T>
    std::enable_if_t<std::is_trivially_copyable_v<T>> read(std::istream &stream, T *data, size_t elements);

    template<class T>
    std::enable_if_t<!std::is_trivially_copyable_v<T>> read(std::istream &stream, T *data, size_t elements);

    template<class T>
    void read(std::istream &stream, std::optional<T> &opt);

    template<class T>
    void read(std::istream &stream, std::vector<T> &vec);

    template<class T>
    void read(std::istream &stream, std::set<T> &vec);

    template<class T>
    void read(std::istream &stream, Gadgetron::hoNDArray<T> &array);

    template<class T, unsigned int D>
    void read(std::istream &stream, Gadgetron::hoNDImage<T, D> &image);

    template<class T, unsigned int D>
    void read(std::istream &stream, Gadgetron::hoNDArray< Gadgetron::hoNDImage<T, D> > &array);

    template<class... ARGS>
    void read(std::istream& stream, std::tuple<ARGS...>& tup);

    template<class T>
    std::enable_if_t<boost::hana::Struct<T>::value> read(std::istream &istream, T &x);



    template<class T>
    T read(std::istream &stream) {
        T val;
        IO::read(stream, val);
        return val;
    }

    template<class TObjectType>
    void read_many(std::istream &stream, Gadgetron::hoNDArray<TObjectType> &array);

    template<class T = uint64_t >
    std::string read_string_from_stream(std::istream &stream);

    template<class T, class V = std::enable_if_t<std::is_trivially_copyable_v<T>>>
    void write(std::ostream &stream, const T &value);

    template<class T>
    void write(std::ostream &ostream, const std::optional<T> &val);

    template<class T>
    void write(std::ostream &ostream, const std::vector<T> &val);

    template<class T>
    void write(std::ostream &ostream, const std::set<T> &val);

    template<class... ARGS>
    void write(std::ostream& ostream, const std::tuple<ARGS...>& tup);

    template<class T>
    std::enable_if_t<boost::hana::Struct<T>::value, void> write(std::ostream &ostream, const T &x);

    template<class T>
    std::enable_if_t<std::is_trivially_copyable_v<T>>
    write(std::ostream &stream, const T *data, size_t number_of_elements);

    template<class T>
    std::enable_if_t<!std::is_trivially_copyable_v<T>>
    write(std::ostream &stream, const T *data, size_t number_of_elements);

    template<class iterator_type>
    using enable_if_forward_iterator = std::enable_if_t<std::is_base_of<std::forward_iterator_tag, typename std::iterator_traits<iterator_type>::iterator_category>::value>;

    template<class iterator_type, class = enable_if_forward_iterator<iterator_type>>
    void write(std::ostream &stream, iterator_type begin, iterator_type end) {
        for (; begin!=end; begin++) {
            IO::write(stream, *begin);
        }
    }

    template<class T>
    void write(std::ostream &stream, const Gadgetron::hoNDArray<T> &array);

    template<class T, unsigned int D>
    void write(std::ostream &stream, const Gadgetron::hoNDImage<T, D> &image);

    template<class T, unsigned int D>
    void write(std::ostream &stream, const Gadgetron::hoNDArray< Gadgetron::hoNDImage<T, D> > &array);

    template<class TObjectType>
    void write_many(std::ostream &stream, const Gadgetron::hoNDArray<TObjectType> &array);

    template<class T = uint64_t>
    void write_string_to_stream(std::ostream &stream, const std::string &str);
}

#include "primitives.hpp"