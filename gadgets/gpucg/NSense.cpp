
#include "NSense.hpp"
#include "preprocess.hpp"

#include <vector_functions.h>

// Preprocessing

template< class UINTd, class FLOATd, char TYPE > mr_recon::NFFT_iteration_plan<UINTd, FLOATd, TYPE>*
mr_recon::preprocess_NSense( UINTd matrix_size, UINTd matrix_size_os, UINTd fixed_dims, UINTd domain_size_grid, unsigned int domain_size_samples, unsigned int domain_size_coils, float W, RealFloatArray *trajectory )
{
  return preprocess_NFFT<UINTd, FLOATd, TYPE>( matrix_size, matrix_size_os, fixed_dims, domain_size_grid, domain_size_samples, domain_size_coils, W, trajectory );
}
 
// Instantiation
template mr_recon::NFFT_iteration_plan< uint2, float2, 0>* mr_recon::preprocess_NSense< uint2,float2, 0>( uint2, uint2, uint2, uint2, unsigned int, unsigned int, float, mr_recon::RealFloatArray* );
template mr_recon::NFFT_iteration_plan< uint3, float3, 0>* mr_recon::preprocess_NSense< uint3,float3, 0>( uint3, uint3, uint3, uint3, unsigned int, unsigned int, float, mr_recon::RealFloatArray* );
template mr_recon::NFFT_iteration_plan< uint4, float4, 0>* mr_recon::preprocess_NSense< uint4,float4, 0>( uint4, uint4, uint4, uint4, unsigned int, unsigned int, float, mr_recon::RealFloatArray* );

template mr_recon::NFFT_iteration_plan< uint2, float2, 1>* mr_recon::preprocess_NSense< uint2,float2, 1>( uint2, uint2, uint2, uint2, unsigned int, unsigned int, float, mr_recon::RealFloatArray* );
template mr_recon::NFFT_iteration_plan< uint3, float3, 1>* mr_recon::preprocess_NSense< uint3,float3, 1>( uint3, uint3, uint3, uint3, unsigned int, unsigned int, float, mr_recon::RealFloatArray* );
template mr_recon::NFFT_iteration_plan< uint4, float4, 1>* mr_recon::preprocess_NSense< uint4,float4, 1>( uint4, uint4, uint4, uint4, unsigned int, unsigned int, float, mr_recon::RealFloatArray* );
