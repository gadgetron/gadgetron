/*--------------------------------------------------------------------------*\
Copyright (c) 2008-2013, Danny Ruijters. All rights reserved.
http://www.dannyruijters.nl/cubicinterpolation/
This file is part of CUDA Cubic B-Spline Interpolation (CI).

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
*  Redistributions of source code must retain the above copyright
   notice, this list of conditions and the following disclaimer.
*  Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions and the following disclaimer in the
   documentation and/or other materials provided with the distribution.
*  Neither the name of the copyright holders nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.

The views and conclusions contained in the software and documentation are
those of the authors and should not be interpreted as representing official
policies, either expressed or implied.

When using this code in a scientific project, please cite one or all of the
following papers:
*  Daniel Ruijters and Philippe Thévenaz,
   GPU Prefilter for Accurate Cubic B-Spline Interpolation, 
   The Computer Journal, vol. 55, no. 1, pp. 15-20, January 2012.
   http://dannyruijters.nl/docs/cudaPrefilter3.pdf
*  Daniel Ruijters, Bart M. ter Haar Romeny, and Paul Suetens,
   Efficient GPU-Based Texture Interpolation using Uniform B-Splines,
   Journal of Graphics Tools, vol. 13, no. 4, pp. 61-69, 2008.
\*--------------------------------------------------------------------------*/

#ifndef _EXTRA_ARGS
#define _EXTRA_ARGS
#define _PASS_EXTRA_ARGS
#endif

#ifndef _TEX2D
#define _TEX2D tex2D
#define _TEXTYPE 2
#endif

//! Bicubic interpolated texture lookup, using unnormalized coordinates.
//! Fast implementation, using 4 bilinear lookups.
//! @param tex  2D texture
//! @param x  unnormalized x texture coordinate
//! @param y  unnormalized y texture coordinate
template<class floatN, class T, enum cudaTextureReadMode mode>
__device__ floatN CUBICTEX2D(texture<T, _TEXTYPE, mode> tex, float x, float y _EXTRA_ARGS)
{
	// transform the coordinate from [0,extent] to [-0.5, extent-0.5]
	const float2 coord_grid = make_float2(x - 0.5f, y - 0.5f);
	const float2 index = floor(coord_grid);
	const float2 fraction = coord_grid - index;
	float2 w0, w1, w2, w3;
	WEIGHTS(fraction, w0, w1, w2, w3);

	const float2 g0 = w0 + w1;
	const float2 g1 = w2 + w3;
	const float2 h0 = (w1 / g0) - make_float2(0.5f) + index;  //h0 = w1/g0 - 1, move from [-0.5, extent-0.5] to [0, extent]
	const float2 h1 = (w3 / g1) + make_float2(1.5f) + index;  //h1 = w3/g1 + 1, move from [-0.5, extent-0.5] to [0, extent]

	// fetch the four linear interpolations
	floatN tex00 = _TEX2D(tex, h0.x, h0.y _PASS_EXTRA_ARGS);
	floatN tex10 = _TEX2D(tex, h1.x, h0.y _PASS_EXTRA_ARGS);
	floatN tex01 = _TEX2D(tex, h0.x, h1.y _PASS_EXTRA_ARGS);
	floatN tex11 = _TEX2D(tex, h1.x, h1.y _PASS_EXTRA_ARGS);

	// weigh along the y-direction
	tex00 = g0.y * tex00 + g1.y * tex01;
	tex10 = g0.y * tex10 + g1.y * tex11;

	// weigh along the x-direction
	return (g0.x * tex00 + g1.x * tex10);
}


// Specializations

// These specializations fill in the floatN and T class types and therefore
// allow the cubicTex2D/cubicTex2DLayered function to be called without any
// template arguments, thus without any <> brackets.

// 1-dimensional pixels
template<enum cudaTextureReadMode mode> __device__ float CUBICTEX2D(texture<float, _TEXTYPE, mode> tex, float x, float y _EXTRA_ARGS) {return CUBICTEX2D<float, float, mode>(tex, x, y _PASS_EXTRA_ARGS);}
template<enum cudaTextureReadMode mode> __device__ float CUBICTEX2D(texture<uchar, _TEXTYPE, mode> tex, float x, float y _EXTRA_ARGS) {return CUBICTEX2D<float, uchar, mode>(tex, x, y _PASS_EXTRA_ARGS);}
template<enum cudaTextureReadMode mode> __device__ float CUBICTEX2D(texture<char, _TEXTYPE, mode> tex, float x, float y _EXTRA_ARGS) {return CUBICTEX2D<float, char, mode>(tex, x, y _PASS_EXTRA_ARGS);}
template<enum cudaTextureReadMode mode> __device__ float CUBICTEX2D(texture<ushort, _TEXTYPE, mode> tex, float x, float y _EXTRA_ARGS) {return CUBICTEX2D<float, ushort, mode>(tex, x, y _PASS_EXTRA_ARGS);}
template<enum cudaTextureReadMode mode> __device__ float CUBICTEX2D(texture<short, _TEXTYPE, mode> tex, float x, float y _EXTRA_ARGS) {return CUBICTEX2D<float, short, mode>(tex, x, y _PASS_EXTRA_ARGS);}
template<enum cudaTextureReadMode mode> __device__ float CUBICTEX2D(texture<uint, _TEXTYPE, mode> tex, float x, float y _EXTRA_ARGS) {return CUBICTEX2D<float, uint, mode>(tex, x, y _PASS_EXTRA_ARGS);}
template<enum cudaTextureReadMode mode> __device__ float CUBICTEX2D(texture<int, _TEXTYPE, mode> tex, float x, float y _EXTRA_ARGS) {return CUBICTEX2D<float, int, mode>(tex, x, y _PASS_EXTRA_ARGS);}
// 2-dimensional pixels
template<enum cudaTextureReadMode mode> __device__ float2 CUBICTEX2D(texture<float2, _TEXTYPE, mode> tex, float x, float y _EXTRA_ARGS) {return CUBICTEX2D<float2, float2, mode>(tex, x, y _PASS_EXTRA_ARGS);}
template<enum cudaTextureReadMode mode> __device__ float2 CUBICTEX2D(texture<uchar2, _TEXTYPE, mode> tex, float x, float y _EXTRA_ARGS) {return CUBICTEX2D<float2, uchar2, mode>(tex, x, y _PASS_EXTRA_ARGS);}
template<enum cudaTextureReadMode mode> __device__ float2 CUBICTEX2D(texture<char2, _TEXTYPE, mode> tex, float x, float y _EXTRA_ARGS) {return CUBICTEX2D<float2, char2, mode>(tex, x, y _PASS_EXTRA_ARGS);}
template<enum cudaTextureReadMode mode> __device__ float2 CUBICTEX2D(texture<ushort2, _TEXTYPE, mode> tex, float x, float y _EXTRA_ARGS) {return CUBICTEX2D<float2, ushort2, mode>(tex, x, y _PASS_EXTRA_ARGS);}
template<enum cudaTextureReadMode mode> __device__ float2 CUBICTEX2D(texture<short2, _TEXTYPE, mode> tex, float x, float y _EXTRA_ARGS) {return CUBICTEX2D<float2, short2, mode>(tex, x, y _PASS_EXTRA_ARGS);}
template<enum cudaTextureReadMode mode> __device__ float2 CUBICTEX2D(texture<uint2, _TEXTYPE, mode> tex, float x, float y _EXTRA_ARGS) {return CUBICTEX2D<float2, uint2, mode>(tex, x, y _PASS_EXTRA_ARGS);}
template<enum cudaTextureReadMode mode> __device__ float2 CUBICTEX2D(texture<int2, _TEXTYPE, mode> tex, float x, float y _EXTRA_ARGS) {return CUBICTEX2D<float2, int2, mode>(tex, x, y _PASS_EXTRA_ARGS);}
// 3-dimensional pixels
template<enum cudaTextureReadMode mode> __device__ float3 CUBICTEX2D(texture<float3, _TEXTYPE, mode> tex, float x, float y _EXTRA_ARGS) {return CUBICTEX2D<float3, float3, mode>(tex, x, y _PASS_EXTRA_ARGS);}
template<enum cudaTextureReadMode mode> __device__ float3 CUBICTEX2D(texture<uchar3, _TEXTYPE, mode> tex, float x, float y _EXTRA_ARGS) {return CUBICTEX2D<float3, uchar3, mode>(tex, x, y _PASS_EXTRA_ARGS);}
template<enum cudaTextureReadMode mode> __device__ float3 CUBICTEX2D(texture<char3, _TEXTYPE, mode> tex, float x, float y _EXTRA_ARGS) {return CUBICTEX2D<float3, char3, mode>(tex, x, y _PASS_EXTRA_ARGS);}
template<enum cudaTextureReadMode mode> __device__ float3 CUBICTEX2D(texture<ushort3, _TEXTYPE, mode> tex, float x, float y _EXTRA_ARGS) {return CUBICTEX2D<float3, ushort3, mode>(tex, x, y _PASS_EXTRA_ARGS);}
template<enum cudaTextureReadMode mode> __device__ float3 CUBICTEX2D(texture<short3, _TEXTYPE, mode> tex, float x, float y _EXTRA_ARGS) {return CUBICTEX2D<float3, short3, mode>(tex, x, y _PASS_EXTRA_ARGS);}
template<enum cudaTextureReadMode mode> __device__ float3 CUBICTEX2D(texture<uint3, _TEXTYPE, mode> tex, float x, float y _EXTRA_ARGS) {return CUBICTEX2D<float3, uint3, mode>(tex, x, y _PASS_EXTRA_ARGS);}
template<enum cudaTextureReadMode mode> __device__ float3 CUBICTEX2D(texture<int3, _TEXTYPE, mode> tex, float x, float y _EXTRA_ARGS) {return CUBICTEX2D<float3, int3, mode>(tex, x, y _PASS_EXTRA_ARGS);}
// 4-dimensional pixels
template<enum cudaTextureReadMode mode> __device__ float4 CUBICTEX2D(texture<float4, _TEXTYPE, mode> tex, float x, float y _EXTRA_ARGS) {return CUBICTEX2D<float4, float4, mode>(tex, x, y _PASS_EXTRA_ARGS);}
template<enum cudaTextureReadMode mode> __device__ float4 CUBICTEX2D(texture<uchar4, _TEXTYPE, mode> tex, float x, float y _EXTRA_ARGS) {return CUBICTEX2D<float4, uchar4, mode>(tex, x, y _PASS_EXTRA_ARGS);}
template<enum cudaTextureReadMode mode> __device__ float4 CUBICTEX2D(texture<char4, _TEXTYPE, mode> tex, float x, float y _EXTRA_ARGS) {return CUBICTEX2D<float4, char4, mode>(tex, x, y _PASS_EXTRA_ARGS);}
template<enum cudaTextureReadMode mode> __device__ float4 CUBICTEX2D(texture<ushort4, _TEXTYPE, mode> tex, float x, float y _EXTRA_ARGS) {return CUBICTEX2D<float4, ushort4, mode>(tex, x, y _PASS_EXTRA_ARGS);}
template<enum cudaTextureReadMode mode> __device__ float4 CUBICTEX2D(texture<short4, _TEXTYPE, mode> tex, float x, float y _EXTRA_ARGS) {return CUBICTEX2D<float4, short4, mode>(tex, x, y _PASS_EXTRA_ARGS);}
template<enum cudaTextureReadMode mode> __device__ float4 CUBICTEX2D(texture<uint4, _TEXTYPE, mode> tex, float x, float y _EXTRA_ARGS) {return CUBICTEX2D<float4, uint4, mode>(tex, x, y _PASS_EXTRA_ARGS);}
template<enum cudaTextureReadMode mode> __device__ float4 CUBICTEX2D(texture<int4, _TEXTYPE, mode> tex, float x, float y _EXTRA_ARGS) {return CUBICTEX2D<float4, int4, mode>(tex, x, y _PASS_EXTRA_ARGS);}

// Function definitions with (float2 coord) instead of (float x, float y)
// 1-dimensional pixels
template<enum cudaTextureReadMode mode> __device__ float CUBICTEX2D(texture<float, _TEXTYPE, mode> tex, float2 coord _EXTRA_ARGS) {return CUBICTEX2D<float, float, mode>(tex, coord.x, coord.y _PASS_EXTRA_ARGS);}
template<enum cudaTextureReadMode mode> __device__ float CUBICTEX2D(texture<uchar, _TEXTYPE, mode> tex, float2 coord _EXTRA_ARGS) {return CUBICTEX2D<float, uchar, mode>(tex, coord.x, coord.y _PASS_EXTRA_ARGS);}
template<enum cudaTextureReadMode mode> __device__ float CUBICTEX2D(texture<char, _TEXTYPE, mode> tex, float2 coord _EXTRA_ARGS) {return CUBICTEX2D<float, char, mode>(tex, coord.x, coord.y _PASS_EXTRA_ARGS);}
template<enum cudaTextureReadMode mode> __device__ float CUBICTEX2D(texture<ushort, _TEXTYPE, mode> tex, float2 coord _EXTRA_ARGS) {return CUBICTEX2D<float, ushort, mode>(tex, coord.x, coord.y _PASS_EXTRA_ARGS);}
template<enum cudaTextureReadMode mode> __device__ float CUBICTEX2D(texture<short, _TEXTYPE, mode> tex, float2 coord _EXTRA_ARGS) {return CUBICTEX2D<float, short, mode>(tex, coord.x, coord.y _PASS_EXTRA_ARGS);}
template<enum cudaTextureReadMode mode> __device__ float CUBICTEX2D(texture<uint, _TEXTYPE, mode> tex, float2 coord _EXTRA_ARGS) {return CUBICTEX2D<float, uint, mode>(tex, coord.x, coord.y _PASS_EXTRA_ARGS);}
template<enum cudaTextureReadMode mode> __device__ float CUBICTEX2D(texture<int, _TEXTYPE, mode> tex, float2 coord _EXTRA_ARGS) {return CUBICTEX2D<float, int, mode>(tex, coord.x, coord.y _PASS_EXTRA_ARGS);}
// 2-dimensional pixels
template<enum cudaTextureReadMode mode> __device__ float2 CUBICTEX2D(texture<float2, _TEXTYPE, mode> tex, float2 coord _EXTRA_ARGS) {return CUBICTEX2D<float2, float2, mode>(tex, coord.x, coord.y _PASS_EXTRA_ARGS);}
template<enum cudaTextureReadMode mode> __device__ float2 CUBICTEX2D(texture<uchar2, _TEXTYPE, mode> tex, float2 coord _EXTRA_ARGS) {return CUBICTEX2D<float2, uchar2, mode>(tex, coord.x, coord.y _PASS_EXTRA_ARGS);}
template<enum cudaTextureReadMode mode> __device__ float2 CUBICTEX2D(texture<char2, _TEXTYPE, mode> tex, float2 coord _EXTRA_ARGS) {return CUBICTEX2D<float2, char2, mode>(tex, coord.x, coord.y _PASS_EXTRA_ARGS);}
template<enum cudaTextureReadMode mode> __device__ float2 CUBICTEX2D(texture<ushort2, _TEXTYPE, mode> tex, float2 coord _EXTRA_ARGS) {return CUBICTEX2D<float2, ushort2, mode>(tex, coord.x, coord.y _PASS_EXTRA_ARGS);}
template<enum cudaTextureReadMode mode> __device__ float2 CUBICTEX2D(texture<short2, _TEXTYPE, mode> tex, float2 coord _EXTRA_ARGS) {return CUBICTEX2D<float2, short2, mode>(tex, coord.x, coord.y _PASS_EXTRA_ARGS);}
template<enum cudaTextureReadMode mode> __device__ float2 CUBICTEX2D(texture<uint2, _TEXTYPE, mode> tex, float2 coord _EXTRA_ARGS) {return CUBICTEX2D<float2, uint2, mode>(tex, coord.x, coord.y _PASS_EXTRA_ARGS);}
template<enum cudaTextureReadMode mode> __device__ float2 CUBICTEX2D(texture<int2, _TEXTYPE, mode> tex, float2 coord _EXTRA_ARGS) {return CUBICTEX2D<float2, int2, mode>(tex, coord.x, coord.y _PASS_EXTRA_ARGS);}
// 3-dimensional pixels
template<enum cudaTextureReadMode mode> __device__ float3 CUBICTEX2D(texture<float3, _TEXTYPE, mode> tex, float2 coord _EXTRA_ARGS) {return CUBICTEX2D<float3, float3, mode>(tex, coord.x, coord.y _PASS_EXTRA_ARGS);}
template<enum cudaTextureReadMode mode> __device__ float3 CUBICTEX2D(texture<uchar3, _TEXTYPE, mode> tex, float2 coord _EXTRA_ARGS) {return CUBICTEX2D<float3, uchar3, mode>(tex, coord.x, coord.y _PASS_EXTRA_ARGS);}
template<enum cudaTextureReadMode mode> __device__ float3 CUBICTEX2D(texture<char3, _TEXTYPE, mode> tex, float2 coord _EXTRA_ARGS) {return CUBICTEX2D<float3, char3, mode>(tex, coord.x, coord.y _PASS_EXTRA_ARGS);}
template<enum cudaTextureReadMode mode> __device__ float3 CUBICTEX2D(texture<ushort3, _TEXTYPE, mode> tex, float2 coord _EXTRA_ARGS) {return CUBICTEX2D<float3, ushort3, mode>(tex, coord.x, coord.y _PASS_EXTRA_ARGS);}
template<enum cudaTextureReadMode mode> __device__ float3 CUBICTEX2D(texture<short3, _TEXTYPE, mode> tex, float2 coord _EXTRA_ARGS) {return CUBICTEX2D<float3, short3, mode>(tex, coord.x, coord.y _PASS_EXTRA_ARGS);}
template<enum cudaTextureReadMode mode> __device__ float3 CUBICTEX2D(texture<uint3, _TEXTYPE, mode> tex, float2 coord _EXTRA_ARGS) {return CUBICTEX2D<float3, uint3, mode>(tex, coord.x, coord.y _PASS_EXTRA_ARGS);}
template<enum cudaTextureReadMode mode> __device__ float3 CUBICTEX2D(texture<int3, _TEXTYPE, mode> tex, float2 coord _EXTRA_ARGS) {return CUBICTEX2D<float3, int3, mode>(tex, coord.x, coord.y _PASS_EXTRA_ARGS);}
// 4-dimensional pixels
template<enum cudaTextureReadMode mode> __device__ float4 CUBICTEX2D(texture<float4, _TEXTYPE, mode> tex, float2 coord _EXTRA_ARGS) {return CUBICTEX2D<float4, float4, mode>(tex, coord.x, coord.y _PASS_EXTRA_ARGS);}
template<enum cudaTextureReadMode mode> __device__ float4 CUBICTEX2D(texture<uchar4, _TEXTYPE, mode> tex, float2 coord _EXTRA_ARGS) {return CUBICTEX2D<float4, uchar4, mode>(tex, coord.x, coord.y _PASS_EXTRA_ARGS);}
template<enum cudaTextureReadMode mode> __device__ float4 CUBICTEX2D(texture<char4, _TEXTYPE, mode> tex, float2 coord _EXTRA_ARGS) {return CUBICTEX2D<float4, char4, mode>(tex, coord.x, coord.y _PASS_EXTRA_ARGS);}
template<enum cudaTextureReadMode mode> __device__ float4 CUBICTEX2D(texture<ushort4, _TEXTYPE, mode> tex, float2 coord _EXTRA_ARGS) {return CUBICTEX2D<float4, ushort4, mode>(tex, coord.x, coord.y _PASS_EXTRA_ARGS);}
template<enum cudaTextureReadMode mode> __device__ float4 CUBICTEX2D(texture<short4, _TEXTYPE, mode> tex, float2 coord _EXTRA_ARGS) {return CUBICTEX2D<float4, short4, mode>(tex, coord.x, coord.y _PASS_EXTRA_ARGS);}
template<enum cudaTextureReadMode mode> __device__ float4 CUBICTEX2D(texture<uint4, _TEXTYPE, mode> tex, float2 coord _EXTRA_ARGS) {return CUBICTEX2D<float4, uint4, mode>(tex, coord.x, coord.y _PASS_EXTRA_ARGS);}
template<enum cudaTextureReadMode mode> __device__ float4 CUBICTEX2D(texture<int4, _TEXTYPE, mode> tex, float2 coord _EXTRA_ARGS) {return CUBICTEX2D<float4, int4, mode>(tex, coord.x, coord.y _PASS_EXTRA_ARGS);}
