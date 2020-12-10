#include "cuNDArray.h"
#include "vector_td.h"
#include <sstream>

namespace Gadgetron{


template void cuNDArray_permute<>(cuNDArray<int>* in,
				 cuNDArray<int>* out,
				 std::vector<unsigned int> *order,
				 int shift_mode);

template void cuNDArray_permute<>(cuNDArray<int2>* in,
				 cuNDArray<int2>* out,
				 std::vector<unsigned int> *order,
				 int shift_mode);

template void cuNDArray_permute<>(cuNDArray<int3>* in,
				 cuNDArray<int3>* out,
				 std::vector<unsigned int> *order,
				 int shift_mode);

template void cuNDArray_permute<>(cuNDArray<int4>* in,
				 cuNDArray<int4>* out,
				 std::vector<unsigned int> *order,
				 int shift_mode);

template void cuNDArray_permute<>(cuNDArray<unsigned int>* in,
				 cuNDArray<unsigned int>* out,
				 std::vector<unsigned int> *order,
				 int shift_mode);

template void cuNDArray_permute<>(cuNDArray<uint2>* in,
				 cuNDArray<uint2>* out,
				 std::vector<unsigned int> *order,
				 int shift_mode);

template void cuNDArray_permute<>(cuNDArray<uint3>* in,
				 cuNDArray<uint3>* out,
				 std::vector<unsigned int> *order,
				 int shift_mode);

template void cuNDArray_permute<>(cuNDArray<uint4>* in,
				 cuNDArray<uint4>* out,
				 std::vector<unsigned int> *order,
				 int shift_mode);

template void cuNDArray_permute<>(cuNDArray<float>* in,
				 cuNDArray<float>* out,
				 std::vector<unsigned int> *order,
				 int shift_mode);

template void cuNDArray_permute<>(cuNDArray<float2>* in,
				 cuNDArray<float2>* out,
				 std::vector<unsigned int> *order,
				 int shift_mode);

template void cuNDArray_permute<>(cuNDArray<float3>* in,
				 cuNDArray<float3>* out,
				 std::vector<unsigned int> *order,
				 int shift_mode);

template void cuNDArray_permute<>(cuNDArray<float4>* in,
				 cuNDArray<float4>* out,
				 std::vector<unsigned int> *order,
				 int shift_mode);

template void cuNDArray_permute<>(cuNDArray<double>* in,
				 cuNDArray<double>* out,
				 std::vector<unsigned int> *order,
				 int shift_mode);

template void cuNDArray_permute<>(cuNDArray<double2>* in,
				 cuNDArray<double2>* out,
				 std::vector<unsigned int> *order,
				 int shift_mode);

template void cuNDArray_permute<>(cuNDArray<double3>* in,
				 cuNDArray<double3>* out,
				 std::vector<unsigned int> *order,
				 int shift_mode);

template void cuNDArray_permute<>(cuNDArray<double4>* in,
				 cuNDArray<double4>* out,
				 std::vector<unsigned int> *order,
				 int shift_mode);



template void cuNDArray_permute<>(cuNDArray<intd1>* in,
				 cuNDArray<intd1>* out,
				 std::vector<unsigned int> *order,
				 int shift_mode);

template void cuNDArray_permute<>(cuNDArray<intd2>* in,
				 cuNDArray<intd2>* out,
				 std::vector<unsigned int> *order,
				 int shift_mode);

template void cuNDArray_permute<>(cuNDArray<intd3>* in,
				 cuNDArray<intd3>* out,
				 std::vector<unsigned int> *order,
				 int shift_mode);

template void cuNDArray_permute<>(cuNDArray<intd4>* in,
				 cuNDArray<intd4>* out,
				 std::vector<unsigned int> *order,
				 int shift_mode);

template void cuNDArray_permute<>(cuNDArray<uint64d1>* in,
				 cuNDArray<uint64d1>* out,
				 std::vector<unsigned int> *order,
				 int shift_mode);

template void cuNDArray_permute<>(cuNDArray<uint64d2>* in,
				 cuNDArray<uint64d2>* out,
				 std::vector<unsigned int> *order,
				 int shift_mode);

template void cuNDArray_permute<>(cuNDArray<uint64d3>* in,
				 cuNDArray<uint64d3>* out,
				 std::vector<unsigned int> *order,
				 int shift_mode);

template void cuNDArray_permute<>(cuNDArray<uint64d4>* in,
				 cuNDArray<uint64d4>* out,
				 std::vector<unsigned int> *order,
				 int shift_mode);

template void cuNDArray_permute<>(cuNDArray<floatd1>* in,
				 cuNDArray<floatd1>* out,
				 std::vector<unsigned int> *order,
				 int shift_mode);

template void cuNDArray_permute<>(cuNDArray<floatd2>* in,
				 cuNDArray<floatd2>* out,
				 std::vector<unsigned int> *order,
				 int shift_mode);

template void cuNDArray_permute<>(cuNDArray<floatd3>* in,
				 cuNDArray<floatd3>* out,
				 std::vector<unsigned int> *order,
				 int shift_mode);

template void cuNDArray_permute<>(cuNDArray<floatd4>* in,
				 cuNDArray<floatd4>* out,
				 std::vector<unsigned int> *order,
				 int shift_mode);

template void cuNDArray_permute<>(cuNDArray<doubled1>* in,
				 cuNDArray<doubled1>* out,
				 std::vector<unsigned int> *order,
				 int shift_mode);

template void cuNDArray_permute<>(cuNDArray<doubled2>* in,
				 cuNDArray<doubled2>* out,
				 std::vector<unsigned int> *order,
				 int shift_mode);

template void cuNDArray_permute<>(cuNDArray<doubled3>* in,
				 cuNDArray<doubled3>* out,
				 std::vector<unsigned int> *order,
				 int shift_mode);

template void cuNDArray_permute<>(cuNDArray<doubled4>* in,
				 cuNDArray<doubled4>* out,
				 std::vector<unsigned int> *order,
				 int shift_mode);
				   
template void cuNDArray_permute<>(cuNDArray<float_complext>* in,
				 cuNDArray<float_complext>* out,
				 std::vector<unsigned int> *order,
				 int shift_mode);

template void cuNDArray_permute<>(cuNDArray<double_complext>* in,
				 cuNDArray<double_complext>* out,
				 std::vector<unsigned int> *order,
				 int shift_mode);
}
