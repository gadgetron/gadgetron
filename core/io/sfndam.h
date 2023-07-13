#pragma once
#include <complex>
#include <cstdint>
#include <istream>
#include <numeric>
#include <ostream>
#include <vector>

/**
 * Implementation of the Simple NDArray with Metadata (SFNDAM) format for exhanging
 * multi-dimensional data with some metadata (a property bag).
 *
 * Overall the format consists of an NDArray of type T and a JSON document with metadata.
 * The binary layout is as follows:
 *
 *    1. Magic bytes 'SFNDAM' followed by a uint32 with the version number.
 *    2. uint32 indicating the data type
 *    3. uint32 indicating the number of dimensions in the array, R
 *    4. Vector of R uint64 values specifying the length of each dimension.
 *       Indexes and sizes are "row first", e.g. dimensions = [z, y, x]
 *       Total number of elements (N) are prod(dimensions).
 *    5. uint64 indicating the length (M) of the metadata (string characters).
 *    6. vector of N*size(data_type) bytes for the array values. Row major storage.
 *    7. M characters of metadata.
 */

namespace sfndam {

const std::uint32_t version = 1;
const auto version_bytes = (char*)(&version);
const char magic_bytes[] = {'S', 'F', 'N', 'D', 'A', 'M', version_bytes[0], version_bytes[1], version_bytes[2], version_bytes[3]};

enum SFNDAM_DataTypes {
  USHORT = 1,  // uint16_t
  SHORT = 2,   // int16_t
  UINT = 3,    // uint32_t
  INT = 4,     // int32_t
  FLOAT = 5,   // float
  DOUBLE = 6,  // double
  CXFLOAT = 7, // complex float
  CXDOUBLE = 8 // complex double
};

template <typename T>
uint32_t typeid_for_type();

template <>
inline uint32_t typeid_for_type<uint16_t>() { return SFNDAM_DataTypes::USHORT; }
template <>
inline uint32_t typeid_for_type<int16_t>() { return SFNDAM_DataTypes::SHORT; }
template <>
inline uint32_t typeid_for_type<uint32_t>() { return SFNDAM_DataTypes::UINT; }
template <>
inline uint32_t typeid_for_type<int32_t>() { return SFNDAM_DataTypes::INT; }
template <>
inline uint32_t typeid_for_type<float>() { return SFNDAM_DataTypes::FLOAT; }
template <>
inline uint32_t typeid_for_type<double>() { return SFNDAM_DataTypes::DOUBLE; }
template <>
inline uint32_t typeid_for_type<std::complex<float>>() { return SFNDAM_DataTypes::CXFLOAT; }
template <>
inline uint32_t typeid_for_type<std::complex<double>>() { return SFNDAM_DataTypes::CXDOUBLE; }

template <typename T>
struct sfndam {
  std::vector<uint64_t> array_dimensions;
  std::vector<T> data;
  std::string meta;
};

template <typename T>
void serialize(const sfndam<T>& a, std::ostream& out) {
  // Check that length of data matches the dimensions
  if (a.data.size() != static_cast<size_t>(std::accumulate(a.array_dimensions.begin(), a.array_dimensions.end(), 1, std::multiplies<uint32_t>()))) {
    throw std::runtime_error("Length of data does not match dimensions");
  }
  // Write magic bytes
  out.write(&magic_bytes[0], sizeof(magic_bytes));

  // Write data type
  auto type_id = typeid_for_type<T>();
  out.write((char*)(&type_id), sizeof(uint32_t));
  // Write number of dimensions
  auto num_dims = static_cast<uint32_t>(a.array_dimensions.size());
  out.write((char*)(&num_dims), sizeof(uint32_t));
  // Write dimensions
  for (auto d : a.array_dimensions) {
    out.write((char*)(&d), sizeof(uint64_t));
  }
  // Write metadata length
  auto meta_len = static_cast<uint64_t>(a.meta.size());
  out.write((char*)(&meta_len), sizeof(uint64_t));
  // Write data
  out.write((char*)(&a.data[0]), a.data.size() * sizeof(T));
  // Write metadata
  out.write(a.meta.c_str(), a.meta.size());

  if (!out.good()) {
    throw std::runtime_error("SFNDAM: Error writing to stream");
  }
}

template <typename T>
sfndam<T> deserialize(std::istream& in) {
  sfndam<T> a;
  // Read magic bytes
  char magic[sizeof(magic_bytes)];
  in.read(magic, sizeof(magic));
  if (memcmp(magic, magic_bytes, sizeof(magic_bytes)) != 0) {
    throw std::runtime_error("Invalid magic bytes");
  }
  // Read data type
  uint32_t type_id;
  in.read((char*)(&type_id), sizeof(uint32_t));
  if (type_id != typeid_for_type<T>()) {
    throw std::runtime_error("Invalid data type");
  }
  // Read number of dimensions
  uint32_t num_dims;
  in.read((char*)(&num_dims), sizeof(uint32_t));
  // Read dimensions
  a.array_dimensions.resize(num_dims);
  for (auto& d : a.array_dimensions) {
    in.read((char*)(&d), sizeof(uint64_t));
  }
  // Read metadata length
  uint64_t meta_len;
  in.read((char*)(&meta_len), sizeof(uint64_t));
  // Read data
  a.data.resize(std::accumulate(a.array_dimensions.begin(), a.array_dimensions.end(), 1, std::multiplies<uint64_t>()));
  in.read((char*)(&a.data[0]), a.data.size() * sizeof(T));
  // Read metadata
  a.meta.resize(meta_len);
  in.read(&a.meta[0], meta_len);

  if (in.fail() || in.bad()) {
    throw std::runtime_error("SFNDAM: Error reading file.");
  }
  return a;
}

}
