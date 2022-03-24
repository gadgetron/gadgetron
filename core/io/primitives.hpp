#pragma once

#include <boost/hana/for_each.hpp>
#include <boost/hana/keys.hpp>
#include <boost/hana/at_key.hpp>

template<class T>
void Gadgetron::Core::IO::write(std::ostream &ostream, const Core::optional<T> &val) {
    IO::write(ostream, bool(val));
    if (val) write(ostream, *val);
}

template<class T>
void Gadgetron::Core::IO::write(std::ostream &ostream, const std::vector<T> &val) {
    IO::write(ostream, val.size());
    IO::write(ostream, val.data(), val.size());
}

template<class T>
void Gadgetron::Core::IO::write(std::ostream &ostream, const std::set<T> &val) {
    IO::write(ostream, val.size());
    IO::write(ostream, val.begin(), val.end());
}

template<class T>
std::enable_if_t<boost::hana::Struct<T>::value, void> Gadgetron::Core::IO::write(std::ostream &ostream, const T &x) {
    namespace hana = boost::hana;
    hana::for_each(hana::keys(x), [&](auto name) {
        IO::write(ostream, hana::at_key(x, name));
    });
}

template<class T, class V>
void Gadgetron::Core::IO::write(std::ostream &stream, const T &value) {
    stream.write(reinterpret_cast<const char *>(&value), sizeof(value));
}

template<class T>
std::enable_if_t<Gadgetron::Core::is_trivially_copyable_v<T>>
Gadgetron::Core::IO::write(std::ostream &stream, const T *data, size_t number_of_elements) {
    stream.write(reinterpret_cast<const char *>(data), number_of_elements * sizeof(T));
}

template<class T>
std::enable_if_t<!Gadgetron::Core::is_trivially_copyable_v<T>>
Gadgetron::Core::IO::write(std::ostream &stream, const T *data, size_t number_of_elements) {
    for (size_t i = 0; i < number_of_elements; i++) {
        write(stream, data[i]);
    }
}

template<class T>
void Gadgetron::Core::IO::write(std::ostream &stream, const hoNDArray <T> &array) {
    write(stream, *array.get_dimensions());
    write(stream, array.get_data_ptr(), array.get_number_of_elements());
}

template<class... ARGS>
void Gadgetron::Core::IO::write(std::ostream &stream, const Core::tuple<ARGS...>& tup) {
    Core::apply([&](const auto&... elements){(write(stream,elements),...);},tup);
}

template<class TObjectType>
void Gadgetron::Core::IO::write_objects_to_stream(std::ostream &stream, const Gadgetron::hoNDArray<TObjectType> &array)
{
    std::vector<size_t> dimensions;
    array.get_dimensions(dimensions);

    Gadgetron::Core::IO::write(stream, dimensions);

    size_t N = array.get_number_of_elements();
    for (size_t i = 0; i < N; i++) {
        Gadgetron::Core::IO::write(stream, array[i]);
    }
}

template<class T>
void Gadgetron::Core::IO::write_string_to_stream(std::ostream &stream, const std::string &str) {
    auto string_length = static_cast<T>(str.size());
    write(stream, string_length);
    stream.write(str.data(), string_length);
}

template<class T, unsigned int D>
void Gadgetron::Core::IO::write(std::ostream &stream, const Gadgetron::hoNDImage<T, D> &array) {

    typedef typename Gadgetron::hoNDImage<T, D>::coord_type coord_type;
    typedef typename Gadgetron::hoNDImage<T, D>::axis_type axis_type;
    
    std::vector<size_t> dimensions;
    std::vector<coord_type> pixelSize;
    std::vector<coord_type> origin;
    axis_type axis;

    array.get_dimensions(dimensions);
    array.get_pixel_size(pixelSize);
    array.get_origin(origin);
    array.get_axis(axis);

    std::vector<coord_type> axis_values(D*D);
    for (auto d=0; d<D; d++) {
        for (auto a=0; a<D; a++) {
            axis_values[a + d*D] = axis[d][a];
        }
    }

    Gadgetron::Core::IO::write(stream, dimensions);
    Gadgetron::Core::IO::write(stream, pixelSize);
    Gadgetron::Core::IO::write(stream, origin);
    Gadgetron::Core::IO::write(stream, axis_values);

    Gadgetron::Core::IO::write(stream, array.get_data_ptr(), array.get_number_of_elements());
}

template<class T, unsigned int D>
void Gadgetron::Core::IO::write(std::ostream &stream, const Gadgetron::hoNDArray< Gadgetron::hoNDImage<T, D> > &array) {
    Gadgetron::Core::IO::write_objects_to_stream(stream, array);
}

template<class T>
std::enable_if_t<Gadgetron::Core::is_trivially_copyable_v<T>> Gadgetron::Core::IO::read(std::istream &stream, T *data, size_t elements) {
    stream.read(reinterpret_cast<char *>(data), elements*sizeof(T));
}

template<class T>
std::enable_if_t<!Gadgetron::Core::is_trivially_copyable_v<T>> Gadgetron::Core::IO::read(std::istream &stream, T *data, size_t elements) {
    for (size_t i = 0; i < elements; i++) read(stream,data[i]);
}

template<class T>
std::string Gadgetron::Core::IO::read_string_from_stream(std::istream &stream) {
    auto n = IO::read<T>(stream);
    std::vector<char> data(n);
    stream.read(data.data(), n);

    return std::string(data.data(),data.size());
}

template<class T>
std::enable_if_t<Gadgetron::Core::is_trivially_copyable_v<T>> Gadgetron::Core::IO::read(std::istream &stream, T &value) {
    stream.read(reinterpret_cast<char *>(&value), sizeof(value));
}

template<class T>
void Gadgetron::Core::IO::read(std::istream &stream, Gadgetron::Core::optional<T> &opt) {
    auto is_set = IO::read<bool>(stream);
    if (!is_set) {
        opt = Core::none;
        return;
    }
    opt = IO::read<T>(stream);
}

template<class T>
void Gadgetron::Core::IO::read(std::istream &stream, std::vector<T> &vec) {
    auto length = IO::read<size_t>(stream);
    vec.resize(length);
    IO::read(stream,vec.data(),length);

}

template<class T>
void Gadgetron::Core::IO::read(std::istream &stream, std::set<T> &setT) {
    auto length = IO::read<size_t>(stream);
    for (size_t i = 0; i < length; i++){
        setT.insert(IO::read<T>(stream));
    }

}

template<class T>
void Gadgetron::Core::IO::read(std::istream &stream, Gadgetron::hoNDArray<T> &array) {
    auto dimensions = IO::read<std::vector<size_t>>(stream);
    array = hoNDArray<T>(dimensions);
    IO::read(stream,array.data(),array.size());

}

template<class T, unsigned int D>
void Gadgetron::Core::IO::read(std::istream &stream, Gadgetron::hoNDImage<T, D> &array) {

    typedef typename Gadgetron::hoNDImage<T, D>::coord_type coord_type;
    typedef typename Gadgetron::hoNDImage<T, D>::axis_type axis_type;

    std::vector<size_t> dimensions;
    std::vector<coord_type> pixelSize, origin, axis_values;

    Gadgetron::Core::IO::read(stream, dimensions);
    Gadgetron::Core::IO::read(stream, pixelSize);
    Gadgetron::Core::IO::read(stream, origin);
    Gadgetron::Core::IO::read(stream, axis_values);

    axis_type axis(D);
    for (auto d=0; d<D; d++) {
        for (auto a=0; a<D; a++) {
            axis[d][a] = axis_values[a + d*D]; 
        }
    }

    array = Gadgetron::hoNDImage<T, D>(dimensions, pixelSize, origin, axis);
    Gadgetron::Core::IO::read(stream, array.data(), array.size());
}

template<class TObjectType>
void Gadgetron::Core::IO::read_objects_from_stream(std::istream &stream, Gadgetron::hoNDArray<TObjectType> &array)
{
    std::vector<size_t> dimensions;
    Gadgetron::Core::IO::read(stream, dimensions);

    array.create(dimensions);

    size_t N = array.get_number_of_elements();
    for (size_t i = 0; i < N; i++) {
        Gadgetron::Core::IO::read(stream, array[i]);
    }
}
    
template<class T, unsigned int D>
void Gadgetron::Core::IO::read(std::istream &stream, Gadgetron::hoNDArray< Gadgetron::hoNDImage<T, D> > &array) {
    Gadgetron::Core::IO::read_objects_from_stream(stream, array);
}

template<class... ARGS>
void Gadgetron::Core::IO::read(std::istream &stream,  Core::tuple<ARGS...>& tup) {
    Core::apply([&](auto&&... elements){(read(stream,elements),...);},tup);
}
template<class T>
std::enable_if_t<boost::hana::Struct<T>::value> Gadgetron::Core::IO::read(std::istream &istream, T &x) {
  namespace hana = boost::hana;
    hana::for_each(hana::keys(x), [&](auto name) {
        IO::read(istream, hana::at_key(x, name));
    });
}
