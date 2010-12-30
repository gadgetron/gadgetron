/*

	-------------------------
	NFFT template instantion.
	-------------------------
	
	This file is included from NFFT.cu. It shall not be compiled individually.
*/

template bool NFFT_initialize( mr_recon::NFFT_plan< uint2, float2, 0>* );
template bool NFFT_initialize( mr_recon::NFFT_plan< uint3, float3, 0>* );
template bool NFFT_initialize( mr_recon::NFFT_plan< uint4, float4, 0>* );

template bool NFFT_initialize( mr_recon::NFFT_plan< uint2, float2, 1>* );
template bool NFFT_initialize( mr_recon::NFFT_plan< uint3, float3, 1>* );
template bool NFFT_initialize( mr_recon::NFFT_plan< uint4, float4, 1>* );

template bool NFFT_initialize( mr_recon::NFFT_plan< uint2, float2, 2>* );
template bool NFFT_initialize( mr_recon::NFFT_plan< uint3, float3, 2>* );
template bool NFFT_initialize( mr_recon::NFFT_plan< uint4, float4, 2>* );

template bool NFFT_initialize( mr_recon::NFFT_H_plan< uint2, float2, 0>* );
template bool NFFT_initialize( mr_recon::NFFT_H_plan< uint3, float3, 0>* );
template bool NFFT_initialize( mr_recon::NFFT_H_plan< uint4, float4, 0>* );

template bool NFFT_initialize( mr_recon::NFFT_H_plan< uint2, float2, 1>* );
template bool NFFT_initialize( mr_recon::NFFT_H_plan< uint3, float3, 1>* );
template bool NFFT_initialize( mr_recon::NFFT_H_plan< uint4, float4, 1>* );

template bool NFFT_initialize( mr_recon::NFFT_H_plan< uint2, float2, 2>* );
template bool NFFT_initialize( mr_recon::NFFT_H_plan< uint3, float3, 2>* );
template bool NFFT_initialize( mr_recon::NFFT_H_plan< uint4, float4, 2>* );

template bool NFFT_initialize( mr_recon::NFFT_iteration_plan< uint2, float2, 0>* );
template bool NFFT_initialize( mr_recon::NFFT_iteration_plan< uint3, float3, 0>* );
template bool NFFT_initialize( mr_recon::NFFT_iteration_plan< uint4, float4, 0>* );

template bool NFFT_initialize( mr_recon::NFFT_iteration_plan< uint2, float2, 1>* );
template bool NFFT_initialize( mr_recon::NFFT_iteration_plan< uint3, float3, 1>* );
template bool NFFT_initialize( mr_recon::NFFT_iteration_plan< uint4, float4, 1>* );

template bool NFFT_initialize( mr_recon::NFFT_iteration_plan< uint2, float2, 2>* );
template bool NFFT_initialize( mr_recon::NFFT_iteration_plan< uint3, float3, 2>* );
template bool NFFT_initialize( mr_recon::NFFT_iteration_plan< uint4, float4, 2>* );

template bool NFFT_compute( mr_recon::NFFT_plan< uint2, float2, 0>*, cuFloatComplex*, bool, cuFloatComplex*, bool );
template bool NFFT_compute( mr_recon::NFFT_plan< uint3, float3, 0>*, cuFloatComplex*, bool, cuFloatComplex*, bool );
template bool NFFT_compute( mr_recon::NFFT_plan< uint4, float4, 0>*, cuFloatComplex*, bool, cuFloatComplex*, bool );

template bool NFFT_compute( mr_recon::NFFT_plan< uint2, float2, 1>*, cuFloatComplex*, bool, cuFloatComplex*, bool );
template bool NFFT_compute( mr_recon::NFFT_plan< uint3, float3, 1>*, cuFloatComplex*, bool, cuFloatComplex*, bool );
template bool NFFT_compute( mr_recon::NFFT_plan< uint4, float4, 1>*, cuFloatComplex*, bool, cuFloatComplex*, bool );

template bool NFFT_compute( mr_recon::NFFT_plan< uint2, float2, 2>*, cuFloatComplex*, bool, cuFloatComplex*, bool );
template bool NFFT_compute( mr_recon::NFFT_plan< uint3, float3, 2>*, cuFloatComplex*, bool, cuFloatComplex*, bool );
template bool NFFT_compute( mr_recon::NFFT_plan< uint4, float4, 2>*, cuFloatComplex*, bool, cuFloatComplex*, bool );

template bool NFFT_compute( mr_recon::NFFT_H_plan< uint2, float2, 0>*, cuFloatComplex*, bool, cuFloatComplex*, bool );
template bool NFFT_compute( mr_recon::NFFT_H_plan< uint3, float3, 0>*, cuFloatComplex*, bool, cuFloatComplex*, bool );
template bool NFFT_compute( mr_recon::NFFT_H_plan< uint4, float4, 0>*, cuFloatComplex*, bool, cuFloatComplex*, bool );

template bool NFFT_compute( mr_recon::NFFT_H_plan< uint2, float2, 1>*, cuFloatComplex*, bool, cuFloatComplex*, bool );
template bool NFFT_compute( mr_recon::NFFT_H_plan< uint3, float3, 1>*, cuFloatComplex*, bool, cuFloatComplex*, bool );
template bool NFFT_compute( mr_recon::NFFT_H_plan< uint4, float4, 1>*, cuFloatComplex*, bool, cuFloatComplex*, bool );

template bool NFFT_compute( mr_recon::NFFT_H_plan< uint2, float2, 2>*, cuFloatComplex*, bool, cuFloatComplex*, bool );
template bool NFFT_compute( mr_recon::NFFT_H_plan< uint3, float3, 2>*, cuFloatComplex*, bool, cuFloatComplex*, bool );
template bool NFFT_compute( mr_recon::NFFT_H_plan< uint4, float4, 2>*, cuFloatComplex*, bool, cuFloatComplex*, bool );

template bool NFFT_compute( mr_recon::NFFT_iteration_plan< uint2, float2, 0>*, cuFloatComplex*, bool, cuFloatComplex*, bool );
template bool NFFT_compute( mr_recon::NFFT_iteration_plan< uint3, float3, 0>*, cuFloatComplex*, bool, cuFloatComplex*, bool );
template bool NFFT_compute( mr_recon::NFFT_iteration_plan< uint4, float4, 0>*, cuFloatComplex*, bool, cuFloatComplex*, bool );

template bool NFFT_compute( mr_recon::NFFT_iteration_plan< uint2, float2, 1>*, cuFloatComplex*, bool, cuFloatComplex*, bool );
template bool NFFT_compute( mr_recon::NFFT_iteration_plan< uint3, float3, 1>*, cuFloatComplex*, bool, cuFloatComplex*, bool );
template bool NFFT_compute( mr_recon::NFFT_iteration_plan< uint4, float4, 1>*, cuFloatComplex*, bool, cuFloatComplex*, bool );

template bool NFFT_compute( mr_recon::NFFT_iteration_plan< uint2, float2, 2>*, cuFloatComplex*, bool, cuFloatComplex*, bool );
template bool NFFT_compute( mr_recon::NFFT_iteration_plan< uint3, float3, 2>*, cuFloatComplex*, bool, cuFloatComplex*, bool );
template bool NFFT_compute( mr_recon::NFFT_iteration_plan< uint4, float4, 2>*, cuFloatComplex*, bool, cuFloatComplex*, bool );

template bool NFFT_cleanup( mr_recon::NFFT_plan< uint2, float2, 0>** );
template bool NFFT_cleanup( mr_recon::NFFT_plan< uint3, float3, 0>** );
template bool NFFT_cleanup( mr_recon::NFFT_plan< uint4, float4, 0>** );

template bool NFFT_cleanup( mr_recon::NFFT_plan< uint2, float2, 1>** );
template bool NFFT_cleanup( mr_recon::NFFT_plan< uint3, float3, 1>** );
template bool NFFT_cleanup( mr_recon::NFFT_plan< uint4, float4, 1>** );

template bool NFFT_cleanup( mr_recon::NFFT_plan< uint2, float2, 2>** );
template bool NFFT_cleanup( mr_recon::NFFT_plan< uint3, float3, 2>** );
template bool NFFT_cleanup( mr_recon::NFFT_plan< uint4, float4, 2>** );

template bool NFFT_cleanup( mr_recon::NFFT_H_plan< uint2, float2, 0>** );
template bool NFFT_cleanup( mr_recon::NFFT_H_plan< uint3, float3, 0>** );
template bool NFFT_cleanup( mr_recon::NFFT_H_plan< uint4, float4, 0>** );

template bool NFFT_cleanup( mr_recon::NFFT_H_plan< uint2, float2, 1>** );
template bool NFFT_cleanup( mr_recon::NFFT_H_plan< uint3, float3, 1>** );
template bool NFFT_cleanup( mr_recon::NFFT_H_plan< uint4, float4, 1>** );

template bool NFFT_cleanup( mr_recon::NFFT_H_plan< uint2, float2, 2>** );
template bool NFFT_cleanup( mr_recon::NFFT_H_plan< uint3, float3, 2>** );
template bool NFFT_cleanup( mr_recon::NFFT_H_plan< uint4, float4, 2>** );

template bool NFFT_cleanup( mr_recon::NFFT_iteration_plan< uint2, float2, 0>** );
template bool NFFT_cleanup( mr_recon::NFFT_iteration_plan< uint3, float3, 0>** );
template bool NFFT_cleanup( mr_recon::NFFT_iteration_plan< uint4, float4, 0>** );

template bool NFFT_cleanup( mr_recon::NFFT_iteration_plan< uint2, float2, 1>** );
template bool NFFT_cleanup( mr_recon::NFFT_iteration_plan< uint3, float3, 1>** );
template bool NFFT_cleanup( mr_recon::NFFT_iteration_plan< uint4, float4, 1>** );

template bool NFFT_cleanup( mr_recon::NFFT_iteration_plan< uint2, float2, 2>** );
template bool NFFT_cleanup( mr_recon::NFFT_iteration_plan< uint3, float3, 2>** );
template bool NFFT_cleanup( mr_recon::NFFT_iteration_plan< uint4, float4, 2>** );

template void deapodize( mr_recon::NFFT_plan< uint2, float2, 0>*, bool, cuFloatComplex*, bool );
template void deapodize( mr_recon::NFFT_plan< uint3, float3, 0>*, bool, cuFloatComplex*, bool );
template void deapodize( mr_recon::NFFT_plan< uint4, float4, 0>*, bool, cuFloatComplex*, bool );

template void deapodize( mr_recon::NFFT_plan< uint2, float2, 1>*, bool, cuFloatComplex*, bool );
template void deapodize( mr_recon::NFFT_plan< uint3, float3, 1>*, bool, cuFloatComplex*, bool );
template void deapodize( mr_recon::NFFT_plan< uint4, float4, 1>*, bool, cuFloatComplex*, bool );

template void deapodize( mr_recon::NFFT_plan< uint2, float2, 2>*, bool, cuFloatComplex*, bool );
template void deapodize( mr_recon::NFFT_plan< uint3, float3, 2>*, bool, cuFloatComplex*, bool );
template void deapodize( mr_recon::NFFT_plan< uint4, float4, 2>*, bool, cuFloatComplex*, bool );

template void deapodize( mr_recon::NFFT_H_plan< uint2, float2, 0>*, bool, cuFloatComplex*, bool );
template void deapodize( mr_recon::NFFT_H_plan< uint3, float3, 0>*, bool, cuFloatComplex*, bool );
template void deapodize( mr_recon::NFFT_H_plan< uint4, float4, 0>*, bool, cuFloatComplex*, bool );

template void deapodize( mr_recon::NFFT_H_plan< uint2, float2, 1>*, bool, cuFloatComplex*, bool );
template void deapodize( mr_recon::NFFT_H_plan< uint3, float3, 1>*, bool, cuFloatComplex*, bool );
template void deapodize( mr_recon::NFFT_H_plan< uint4, float4, 1>*, bool, cuFloatComplex*, bool );

template void deapodize( mr_recon::NFFT_H_plan< uint2, float2, 2>*, bool, cuFloatComplex*, bool );
template void deapodize( mr_recon::NFFT_H_plan< uint3, float3, 2>*, bool, cuFloatComplex*, bool );
template void deapodize( mr_recon::NFFT_H_plan< uint4, float4, 2>*, bool, cuFloatComplex*, bool );

template void deapodize( mr_recon::NFFT_iteration_plan< uint2, float2, 0>*, bool, cuFloatComplex*, bool );
template void deapodize( mr_recon::NFFT_iteration_plan< uint3, float3, 0>*, bool, cuFloatComplex*, bool );
template void deapodize( mr_recon::NFFT_iteration_plan< uint4, float4, 0>*, bool, cuFloatComplex*, bool );

template void deapodize( mr_recon::NFFT_iteration_plan< uint2, float2, 1>*, bool, cuFloatComplex*, bool );
template void deapodize( mr_recon::NFFT_iteration_plan< uint3, float3, 1>*, bool, cuFloatComplex*, bool );
template void deapodize( mr_recon::NFFT_iteration_plan< uint4, float4, 1>*, bool, cuFloatComplex*, bool );

template void deapodize( mr_recon::NFFT_iteration_plan< uint2, float2, 2>*, bool, cuFloatComplex*, bool );
template void deapodize( mr_recon::NFFT_iteration_plan< uint3, float3, 2>*, bool, cuFloatComplex*, bool );
template void deapodize( mr_recon::NFFT_iteration_plan< uint4, float4, 2>*, bool, cuFloatComplex*, bool );

template bool NFFT_iteration_convolve_to_image( mr_recon::NFFT_iteration_plan< uint2, float2, 0>*, cuFloatComplex*, cuFloatComplex*, bool );
template bool NFFT_iteration_convolve_to_image( mr_recon::NFFT_iteration_plan< uint3, float3, 0>*, cuFloatComplex*, cuFloatComplex*, bool );
template bool NFFT_iteration_convolve_to_image( mr_recon::NFFT_iteration_plan< uint4, float4, 0>*, cuFloatComplex*, cuFloatComplex*, bool );

template bool NFFT_iteration_convolve_to_image( mr_recon::NFFT_iteration_plan< uint2, float2, 1>*, cuFloatComplex*, cuFloatComplex*, bool );
template bool NFFT_iteration_convolve_to_image( mr_recon::NFFT_iteration_plan< uint3, float3, 1>*, cuFloatComplex*, cuFloatComplex*, bool );
template bool NFFT_iteration_convolve_to_image( mr_recon::NFFT_iteration_plan< uint4, float4, 1>*, cuFloatComplex*, cuFloatComplex*, bool );

template bool NFFT_iteration_convolve_to_image( mr_recon::NFFT_iteration_plan< uint2, float2, 2>*, cuFloatComplex*, cuFloatComplex*, bool );
template bool NFFT_iteration_convolve_to_image( mr_recon::NFFT_iteration_plan< uint3, float3, 2>*, cuFloatComplex*, cuFloatComplex*, bool );
template bool NFFT_iteration_convolve_to_image( mr_recon::NFFT_iteration_plan< uint4, float4, 2>*, cuFloatComplex*, cuFloatComplex*, bool );

template bool NFFT_iteration_convolve_to_samples( mr_recon::NFFT_iteration_plan< uint2, float2, 0>*, cuFloatComplex*, cuFloatComplex* );
template bool NFFT_iteration_convolve_to_samples( mr_recon::NFFT_iteration_plan< uint3, float3, 0>*, cuFloatComplex*, cuFloatComplex* );
template bool NFFT_iteration_convolve_to_samples( mr_recon::NFFT_iteration_plan< uint4, float4, 0>*, cuFloatComplex*, cuFloatComplex* );

template bool NFFT_iteration_convolve_to_samples( mr_recon::NFFT_iteration_plan< uint2, float2, 1>*, cuFloatComplex*, cuFloatComplex* );
template bool NFFT_iteration_convolve_to_samples( mr_recon::NFFT_iteration_plan< uint3, float3, 1>*, cuFloatComplex*, cuFloatComplex* );
template bool NFFT_iteration_convolve_to_samples( mr_recon::NFFT_iteration_plan< uint4, float4, 1>*, cuFloatComplex*, cuFloatComplex* );

template bool NFFT_iteration_convolve_to_samples( mr_recon::NFFT_iteration_plan< uint2, float2, 2>*, cuFloatComplex*, cuFloatComplex* );
template bool NFFT_iteration_convolve_to_samples( mr_recon::NFFT_iteration_plan< uint3, float3, 2>*, cuFloatComplex*, cuFloatComplex* );
template bool NFFT_iteration_convolve_to_samples( mr_recon::NFFT_iteration_plan< uint4, float4, 2>*, cuFloatComplex*, cuFloatComplex* );

template float2* compute_radial_sample_positions( mr_recon::NFFT_iteration_plan<uint2, float2, 1>*, bool, bool );
template float3* compute_radial_sample_positions( mr_recon::NFFT_iteration_plan<uint3, float3, 1>*, bool, bool );
template float4* compute_radial_sample_positions( mr_recon::NFFT_iteration_plan<uint4, float4, 1>*, bool, bool );

template float2* compute_radial_sample_positions( mr_recon::NFFT_iteration_plan<uint2, float2, 2>*, bool, bool );
template float3* compute_radial_sample_positions( mr_recon::NFFT_iteration_plan<uint3, float3, 2>*, bool, bool );
template float4* compute_radial_sample_positions( mr_recon::NFFT_iteration_plan<uint4, float4, 2>*, bool, bool );

template float* get_density_compensation_weights( mr_recon::NFFT_iteration_plan<uint2, float2, 1>* );
template float* get_density_compensation_weights( mr_recon::NFFT_iteration_plan<uint3, float3, 1>* );
template float* get_density_compensation_weights( mr_recon::NFFT_iteration_plan<uint4, float4, 1>* );

template float* get_density_compensation_weights( mr_recon::NFFT_iteration_plan<uint2, float2, 2>* );
template float* get_density_compensation_weights( mr_recon::NFFT_iteration_plan<uint3, float3, 2>* );
template float* get_density_compensation_weights( mr_recon::NFFT_iteration_plan<uint4, float4, 2>* );

template float2* get_sample_trajectories( mr_recon::NFFT_iteration_plan<uint2, float2, 1>* );
template float3* get_sample_trajectories( mr_recon::NFFT_iteration_plan<uint3, float3, 1>* );
template float4* get_sample_trajectories( mr_recon::NFFT_iteration_plan<uint4, float4, 1>* );

template float2* get_sample_trajectories( mr_recon::NFFT_iteration_plan<uint2, float2, 2>* );
template float3* get_sample_trajectories( mr_recon::NFFT_iteration_plan<uint3, float3, 2>* );
template float4* get_sample_trajectories( mr_recon::NFFT_iteration_plan<uint4, float4, 2>* );

template float* estimate_dcw( mr_recon::NFFT_iteration_plan<uint2, float2, 0>*, float2* );
