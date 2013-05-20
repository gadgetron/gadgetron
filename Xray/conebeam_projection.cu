#include "conebeam_projection.h"

//#include "float_util.hcu"
#include "check_CUDA.h"

#include <cuda_runtime_api.h>
#include <math_constants.h>
#include <assert.h>
#include <float.h>

#include "hoNDArray.h"
#include "cuNDArray_elemwise.h"
#include "cuNDArray_operators.h"
#include "hoNDArray_fileio.h"
#include "cuNFFT.h"
#include "GPUTimer.h"
#include "float3x3.h"
#include "cuNDArray_utils.h"
//#include "../../libs/device_utils/float3x3.h"
//#include <cuda_utilities.h>
#include <iostream>
#include <fstream>

using namespace Gadgetron;


//const bool printTiming = false;
const bool debugInfo = false;

//#define PI CUDART_PI_F
#define PI (4.0f*atan(1.0f))

//#define zAxisDir -1.0f;
#define zAxisDir 1.0f;

//#define PS_ORIGON_CENTERING
//#define IS_ORIGON_CENTERING

// Read the projection/image data respectively as a texture (taking advantage of the cache and interpolation)
#define NORMALIZED_TC 0
texture<float,  3, cudaReadModeElementType> image_tex( NORMALIZED_TC, cudaFilterModeLinear, cudaAddressModeBorder );
texture<float,  cudaTextureType2DLayered, cudaReadModeElementType> projection_tex( NORMALIZED_TC, cudaFilterModeLinear, cudaAddressModeBorder );

inline __host__ __device__
float degrees2radians(float degree) {
    return degree * (PI/180.0f);
}

inline __host__ __device__
float calcSagValue(float angle_in_degrees, floatd3 sag_parameters) {
    float A = sag_parameters[0];
    float B = sag_parameters[1];
    float C = degrees2radians(angle_in_degrees + sag_parameters[3]);
    return A + B * cos( C );
}

inline __host__ __device__
floatd2 calcSagVector(float ioe_angle_in_degrees, floatd3 sag_parameters_x, floatd3 sag_parameters_y) {
    float d_xp = calcSagValue(ioe_angle_in_degrees, sag_parameters_x);
    float d_yp = calcSagValue(ioe_angle_in_degrees, sag_parameters_y);
    return floatd2(d_xp, d_yp);
}

__global__ void
conebeam_forwards_projection_cb_kernel(float* projections, float* angles,
                                       floatd3 is_dims_in_pixels, floatd3 is_spacing_in_mm,
                                       uint2 ps_dims_in_pixels_uint, floatd2 ps_dims_in_mm,
                                       floatd3 sag_parameters_x, floatd3 sag_parameters_y,
                                       float SDD, float SAD,
                                       const unsigned int numSamplesPerRay, bool use_circular_cutoff) {

    const int PROJ_IDX_X = threadIdx.x + blockIdx.x * blockDim.x;
    const int PROJ_IDX_Y = threadIdx.y + blockIdx.y * blockDim.y;

    const int PROJECTION = blockIdx.z;

    const int PROJECTION_RES_X = ps_dims_in_pixels_uint.x;
    const int PROJECTION_RES_Y = ps_dims_in_pixels_uint.y;
    if ( PROJ_IDX_X < PROJECTION_RES_X && PROJ_IDX_Y < PROJECTION_RES_Y) {

    floatd2 ps_dims_in_pixels = floatd2(PROJECTION_RES_X, PROJECTION_RES_Y);

    floatd2 ps_spacing = ps_dims_in_mm / ps_dims_in_pixels;

    // Determine projection angle and rotation matrix
    const float angle = angles[PROJECTION];
    const float3x3 rotation = calcRotationMatrixAroundZ(degrees2radians(angle));
    floatd2 sag = calcSagVector(angle, sag_parameters_x, sag_parameters_y);

    // Find start and end point for the line integral
    floatd3 startPoint = floatd3(0.0f, -SAD, 0.0f);
    startPoint = mul(rotation, startPoint);

    floatd2 ps_pc = floatd2(PROJ_IDX_X, PROJ_IDX_Y); // pixel coordinates
    /*
#ifdef PS_ORIGIN_CENTERING
    // old origin, between pixels (negative and positive sides have equal number of pixels)
    ps_pc += floatd2(0.5);
#endif
    */
    floatd2 ps_nc01 = ps_pc / ps_dims_in_pixels; // [0; 1[

    //floatd2 ps_nc01_offset = floatd2(0.5);
    floatd2 ps_nc01_offset = floor(ps_dims_in_pixels/2.0) / ps_dims_in_pixels;
    floatd2 ps_nc = ps_nc01 - ps_nc01_offset; // [-0.5; 0.5[
    floatd2 proj_coords = ps_nc * ps_dims_in_mm; // ps in mm
    proj_coords += sag; // this sag is in mm.
    //proj_coords += floatd2(148.0f, 0.0f); // Half-fan

    //floatd2 ps_origin_offset = floatd2(0.0f);
    //proj_coords += ps_origin_offset; // this is in mm.

    const float ADD = SDD - SAD; // in mm.
    floatd3 endPoint = floatd3(proj_coords[0], ADD, proj_coords[1]);
    //floatd3 endPoint = floatd3(proj_coords.y, ADD, proj_coords.x);
    endPoint = mul(rotation, endPoint);

    // Find direction of the line integral
    floatd3 dir = endPoint-startPoint;

    // Scale start and end onto radius
    floatd3 is_dims_in_mm = is_spacing_in_mm * is_dims_in_pixels;
    float radius = norm( (is_dims_in_mm + is_spacing_in_mm) / floatd3(2.0f) );
    floatd3 unitDir = dir / SDD;
    startPoint = startPoint + (unitDir * (SAD-radius));
    endPoint = endPoint - (unitDir * (ADD-radius));
    float rayLength = norm(startPoint-endPoint);

    const int projectionSpaceId =
        PROJ_IDX_X +
        PROJ_IDX_Y * PROJECTION_RES_X +
        PROJECTION * PROJECTION_RES_X * PROJECTION_RES_Y;

    float result = 0.0f;

        for ( int sampleIndex = 0; sampleIndex<numSamplesPerRay; sampleIndex++) {
            float i = float(sampleIndex)/float(numSamplesPerRay); // [0; 1[
            floatd3 samplePoint = startPoint+dir*i;
            /*
            if (use_circular_cutoff) {
                float inplane_inner_circle_radius = fmin(is_dims_in_mm.x, is_dims_in_mm.y) / 2.0f;
                float length_of_p2D = length(floatd2( floor(fabsf(samplePoint.x)),
                                                          floor(fabsf(samplePoint.y)) ));
                if ( length_of_p2D > inplane_inner_circle_radius )
                    continue;
            }
            */

            floatd3 is_nc = samplePoint / is_dims_in_mm; // normalized coordinates [-0.5, 0.5]
            is_nc[2] *= zAxisDir; // swap z coordinate in image space


            floatd3 is_nc01_offset = floor(is_dims_in_pixels/2.0) / is_dims_in_pixels;
            floatd3 is_nc01 = is_nc + is_nc01_offset; // normalized coordinate space offset [0.0, 1.0]
            floatd3 is_pc = is_nc01 * is_dims_in_pixels; // pixel coordinates
            /*
#ifndef IS_ORIGON_CENTERING
            is_pc += floatd3(0.5) * zAxisDir;
#endif
            */

            is_pc += floatd3(0.5f); // into OpenGL texture coordinates

            // Accumulate result
            float value = tex3D( image_tex, is_pc[0], is_pc[1], is_pc[2] );
            result += value * (rayLength / float(numSamplesPerRay));
        } // ends sampleIndex loop

    projections[projectionSpaceId] += result;
    }
}

// Forwards projection of the 3D volume of a given phase onto projections
template <class TYPE>
bool Gadgetron::conebeam_forwards_projection( hoNDArray<TYPE>* projections, hoNDArray<TYPE>* x, unsigned int bin,
                                   std::vector<unsigned int>& binningdata,
                                   std::vector<float>& angles,
                                   unsigned int orig_ppb, floatd3 is_spacing_in_mm, floatd2 ps_dims_in_mm,
                                   floatd3 sag_parameters_x, floatd3 sag_parameters_y,
                                   float SDD, float SAD,
                                   unsigned int numSamplesPerRay, bool use_circular_cutoff,
                                   bool accumulate) {
    assert( projections->get_number_of_dimensions() == 3 );
    assert( (x->get_number_of_dimensions() == 3) || (x->get_number_of_dimensions() == 4) );

    if (accumulate)
        std::cout << "accumulate seams to be broken!" << std::endl;

    const unsigned int projection_res_x = projections->get_size(0);
    const unsigned int projection_res_y = projections->get_size(1);
    const unsigned int num_projections  = binningdata.size();

    const unsigned int matrix_size_x = x->get_size(0);
    const unsigned int matrix_size_y = x->get_size(1);
    const unsigned int matrix_size_z = x->get_size(2);
    floatd3 is_dims_in_pixels = floatd3(matrix_size_x, matrix_size_y, matrix_size_z);

    unsigned int offs = bin * matrix_size_x * matrix_size_y * matrix_size_z;

    // Build array for input texture
    cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
    cudaExtent extent;
    extent.width = matrix_size_x;
    extent.height = matrix_size_y;
    extent.depth = matrix_size_z;

    cudaMemcpy3DParms cpy_params = {0};
    cpy_params.kind = cudaMemcpyHostToDevice;
    cpy_params.extent = extent;

    cudaArray *image_array;
    cudaMalloc3DArray(&image_array, &channelDesc, extent);
    CHECK_FOR_CUDA_ERROR();

    cpy_params.dstArray = image_array;
    cpy_params.srcPtr =
        make_cudaPitchedPtr((void*)(x->get_data_ptr()+offs), extent.width*sizeof(float),
                            extent.width, extent.height);
    cudaMemcpy3D(&cpy_params);
    CHECK_FOR_CUDA_ERROR();

    cudaBindTextureToArray(image_tex, image_array, channelDesc);
    CHECK_FOR_CUDA_ERROR();

    unsigned int num_batches = std::ceil(num_projections / float(orig_ppb));

    if (debugInfo) {
        printf("number of batches: %i\n", num_batches);
        printf("original projections per batch: %i\n", orig_ppb);
        printf("total number then: %i\n", orig_ppb*num_batches);
    }

    for (unsigned int batch = 0; batch < num_batches; batch++) {
        unsigned int from_projection = batch * orig_ppb;
        unsigned int to_projection = (batch+1) * orig_ppb;
        if (from_projection > num_projections)
            from_projection = num_projections;
        if (to_projection > num_projections)
            to_projection = num_projections;

        unsigned int ppb = to_projection-from_projection;

        if (debugInfo) {
            printf("batch: %03i, handling projections: %03i - %03i, angle %.2f - %.2f, ppb %i \n",
                   batch, from_projection, to_projection-1,
                   angles[from_projection], angles[to_projection-1], ppb);
        }

        float* angles_DevPtr;
        cudaMalloc( (void**) &angles_DevPtr, ppb*sizeof(float));
        CHECK_FOR_CUDA_ERROR();

        // Allocate temporary GPU memory for the forward projections
        float* projections_DevPtr;
        unsigned int projSize = projection_res_x*projection_res_y;
        cudaMalloc( (void**) &projections_DevPtr, projSize*ppb*sizeof(float));
        CHECK_FOR_CUDA_ERROR();

        for (unsigned int p=from_projection, i=0; p<to_projection; p++, i++) {
            unsigned int from_id = binningdata[p];
            cudaMemcpy(angles_DevPtr+i, &angles[from_id],
                       sizeof(float), cudaMemcpyHostToDevice);
            CHECK_FOR_CUDA_ERROR();

            if (accumulate) {
                cudaMemcpy(projections_DevPtr+i*projSize, projections->get_data_ptr()+from_id*projSize,
                           projSize*sizeof(float), cudaMemcpyHostToDevice);
                CHECK_FOR_CUDA_ERROR();
            } else {
                cudaMemset(projections_DevPtr+i*projSize, 0, projSize*sizeof(float));
                CHECK_FOR_CUDA_ERROR();
            }

        }



        dim3 dimBlock(256, 1);
		dim3 dimGrid((unsigned int) std::ceil((double)projection_res_x/(double)dimBlock.x),
                     (unsigned int) std::ceil((double)projection_res_y/(double)dimBlock.y), ppb);
        uint2 ps_dims_in_pixels = make_uint2(projection_res_x, projection_res_y);

        cudaFuncSetCacheConfig(conebeam_forwards_projection_cb_kernel, cudaFuncCachePreferL1);

        conebeam_forwards_projection_cb_kernel<<< dimGrid, dimBlock >>>
            (projections_DevPtr, angles_DevPtr,
             is_dims_in_pixels, is_spacing_in_mm, ps_dims_in_pixels, ps_dims_in_mm,
             sag_parameters_x, sag_parameters_y, SDD, SAD, numSamplesPerRay, use_circular_cutoff);
        CHECK_FOR_CUDA_ERROR();

        // Copy result from device to host
        for (unsigned int p=from_projection, i=0; p<to_projection; p++, i++) {
            unsigned int to_id = binningdata[p];
            CUDA_CALL(cudaMemcpy(projections->get_data_ptr()+projSize*to_id,
                       projections_DevPtr+i*projSize, projSize*sizeof(float), cudaMemcpyDeviceToHost));
            CHECK_FOR_CUDA_ERROR();
        }

        cudaFree(projections_DevPtr);
        cudaFree(angles_DevPtr);
    }

    cudaUnbindTexture(image_tex);
    cudaFreeArray(image_array);

    return true;
}

__global__ void
conebeam_backwards_projection_cb_kernel( floatd3 is_dims_in_pixels, floatd2 ps_dims_in_pixels,
                                unsigned int num_projections,
                                float total_num_projections, float *x, float *angles,
                                unsigned int from_projection,
                                floatd2 ps_dims_in_mm, floatd3 is_spacing_in_mm,
                                floatd3 sag_parameters_x, floatd3 sag_parameters_y,
                                float SDD, float SAD, bool use_circular_cutoff,
                                bool use_fbp) {
    const unsigned int IMAGE_Y = blockIdx.x;
    const unsigned int IMAGE_Z = blockIdx.y;
    const unsigned int IMAGE_X = threadIdx.x;

    //const floatd3 is_spacing = is_dims_in_mm / is_dims_in_pixels;
    floatd3 is_dims_in_mm = is_spacing_in_mm * is_dims_in_pixels;
    const floatd2 ps_spacing = ps_dims_in_mm / ps_dims_in_pixels;

    floatd3 is_pc = floatd3( IMAGE_X, IMAGE_Y, IMAGE_Z );


    // Normailzed coordinater ]0, 1[
    floatd3 is_nc01 = is_pc / is_dims_in_pixels;

    // Normailzed coordinater ]-0.5, 0.5[
    floatd3 is_nc01_offset = floor(is_dims_in_pixels/2.0) / is_dims_in_pixels;
    floatd3 is_nc = is_nc01 - is_nc01_offset;
    is_nc[2] *= zAxisDir; // swap Z Axis in image space

    // real world coordinates in mm.
    const floatd3 pos_traj_original = is_nc * is_dims_in_mm;

    // Backproject all projection data
    float result = 0.0f;

    for (int projection = 0; projection < num_projections; projection++ ) {
        const float angle = angles[projection];
        float3x3 inverseRotation = calcRotationMatrixAroundZ(degrees2radians(-angle));
        floatd3 pos_traj = mul(inverseRotation, pos_traj_original);

        floatd3 startPoint = floatd3(0.0f, -SAD, 0.0f);
        floatd3 dir = pos_traj - startPoint;
        dir = dir / dir[1];
        floatd3 endPoint = startPoint + dir * SDD;

        floatd2 endPoint2d = floatd2(endPoint[0], endPoint[2]);

        floatd2 sag = calcSagVector(angle, sag_parameters_x, sag_parameters_y);
        endPoint2d -= sag; // this sag is in mm.

        floatd2 ps_nc = (endPoint2d / ps_dims_in_mm); // [-0.5; 0.5]
        floatd2 ps_nc01_offset = floor(ps_dims_in_pixels/2.0) / ps_dims_in_pixels;
        floatd2 ps_nc01 = ps_nc + ps_nc01_offset; // [0; 1]
        floatd2 ps_pc = ps_nc01 * ps_dims_in_pixels;

        ps_pc += floatd2(0.5f); // into OpenGL texture coordinates

        float weight = 1.0;
        if (use_fbp) {
            const float xx = pos_traj_original[0];
            const float yy = pos_traj_original[1];
            //const float zz = pos_traj_original.z;

            const float beta = degrees2radians(angle);
            const float r = hypot(xx,yy);
            const float phi = atan2(yy,xx);
            // equation 3.59, page 96 and equation 10.2, page 386
            // in Computed Tomography 2nd edition, Jiang Hsieh
            const float D = SAD;
            const float ym = r*sin(beta-phi);
            const float U = (D+ym)/D;
            weight = 1.0f/(U*U);
        }

	// Enforce a boundary condition in which we repeat the first and last slice indefinitely --in the 'z' direction--
	// For the x and y axis it is ok to use the border mode (0 return value) if the FOV is set correctly.

	//if( ps_pc.y < 0.5f ) ps_pc.y = 0.5f;
	//if( ps_pc.y > (ps_dims_in_pixels.y-0.5f) ) ps_pc.y = ps_dims_in_pixels.y-0.5f;

        // Read the projection data (bilinear interpolation enabled) and accumulate
        float val = weight * tex2DLayered( projection_tex, ps_pc[0], ps_pc[1], projection);
        result += val / float(total_num_projections);

    } // end projection loop

    int id =
        IMAGE_X +
        IMAGE_Y * ((unsigned int)is_dims_in_pixels[0]) +
        IMAGE_Z * ((unsigned int)is_dims_in_pixels[0]) * ((unsigned int)is_dims_in_pixels[1]);
    x[ id ] += result;
}

boost::shared_ptr< cuNDArray<float> > get_conebeam_cosinus_weights(uintd2 ps_dims_in_pixels, floatd2 ps_dims_in_mm,
                                                          float D0d, float Ds0) {
    std::vector<unsigned int> dims;
    dims.push_back(ps_dims_in_pixels[0]);
    dims.push_back(ps_dims_in_pixels[1]);

    hoNDArray<float> weights;
    weights.create(&dims);

    const float Dsd = D0d + Ds0;

    float* data = weights.get_data_ptr();
    for (unsigned int x=0; x<ps_dims_in_pixels[0]; x++) {
        for (unsigned int y=0; y<ps_dims_in_pixels[1]; y++) {
            double xx = (( double(x) / double(ps_dims_in_pixels[0])) - 0.5) * ps_dims_in_mm[0];
            double yy = (( double(y) / double(ps_dims_in_pixels[1])) - 0.5) * ps_dims_in_mm[1];
	    double s = Ds0 * xx/Dsd;
	    double t = Ds0 * yy/Dsd;
        // equation 10.1, page 386 in Computed Tomography 2nd edition, Jiang Hsieh
	    double value = Ds0 / sqrt( Ds0*Ds0 + s*s + t*t );
            data[x+y*ps_dims_in_pixels[0]] = value;
        }
    }
    return boost::shared_ptr< cuNDArray<float> >( new cuNDArray<float>(&weights) );
}

boost::shared_ptr< cuNDArray<float> > get_conebeam_ramp(unsigned int dim, float delta) {
    std::vector<unsigned int> dims, dims_os;
    dims.push_back(dim);
    dims_os.push_back(dim<<1);

    hoNDArray<float> weights;

    cuNDArray<float> *res = new cuNDArray<float>();

   weights.create(&dims);


    res->create(&dims);

    float* data = weights.get_data_ptr();
    for (int i=0; i<dim; i++) {
      // equation 3.29, page 73 in Computed Tomography 2nd edition, Jiang Hsieh
      int n = i-dim/2;
      float value;
      if(n==0) {
          value = 1.0f/(4.0f*delta*delta);
      } else if ((i%2)==0) { // even
          value = 0.0f;
      } else { // odd
          float tmp = n*PI*delta;
          value = -1.0f/(tmp*tmp);
      }
      data[i] = value;
    }

    cuNDArray<float> tmp_real(&weights);
    boost::shared_ptr< cuNDArray<float> > weights_os = pad<float,1>( from_std_vector<unsigned int,1>(dims_os), &tmp_real );

    boost::shared_ptr< cuNDArray<float_complext> > tmp_cplx = real_to_complex<float_complext>(weights_os.get());
    cuNFFT_plan<float,1>().fft(tmp_cplx.get(),cuNFFT_plan<float,1>::NFFT_FORWARDS);
    boost::shared_ptr< cuNDArray<float> >  tmp_res = real(tmp_cplx.get());

    //zero_fill_border<float,1>( vector_to_uintd<1>(*weights.get_dimensions())>>1, res.get());

    crop<float,1>( from_std_vector<unsigned int,1>(dims)>>1, tmp_res.get(), res );

    return boost::shared_ptr< cuNDArray<float> >(res);
}

// Backproject all projections for a given phase.
template <class TYPE>
bool Gadgetron::conebeam_backwards_projection( hoNDArray<TYPE>* projections, hoNDArray<TYPE>* x, unsigned int bin,
                                    std::vector<unsigned int>& binningdata,
                                    std::vector<float>& angles,
                                    unsigned int orig_ppb, floatd3 is_spacing_in_mm, floatd2 ps_dims_in_mm,
                                    floatd3 sag_parameters_x, floatd3 sag_parameters_y,
                                    float SDD, float SAD, bool use_circular_cutoff,
                                    bool use_fbp, bool accumulate) {

    /*hoNDArray<TYPE> projections;
    std::vector<unsigned int> image_dims = *(projections_in.get_dimensions().get());
    projections->create(&image_dims);*/

    assert( projections->get_number_of_dimensions() == 3 );
    assert( (x->get_number_of_dimensions() == 3) || (x->get_number_of_dimensions() == 4) );

    unsigned int matrix_size_x = x->get_size(0);
    unsigned int matrix_size_y = x->get_size(1);
    unsigned int matrix_size_z = x->get_size(2);
    floatd3 is_dims = floatd3(matrix_size_x, matrix_size_y, matrix_size_z);
    unsigned int num_image_elements = matrix_size_x*matrix_size_y*matrix_size_z;

    unsigned int projection_res_x = projections->get_size(0);
    unsigned int projection_res_y = projections->get_size(1);
    floatd2 ps_dims_in_pixels = floatd2(projection_res_x, projection_res_y);
    unsigned int num_projections = binningdata.size();

    // BACKWARDS PROJECT
    float *x_DevPtr;
    cudaMalloc( (void**) &x_DevPtr, num_image_elements*sizeof(float) );

    if (accumulate) {
        unsigned int offs = bin * matrix_size_x * matrix_size_y * matrix_size_z;
        cudaMemcpy(x_DevPtr, x->get_data_ptr()+offs,
                   num_image_elements*sizeof(float), cudaMemcpyHostToDevice);
        CHECK_FOR_CUDA_ERROR();
    } else {
        cudaMemset(x_DevPtr, 0, num_image_elements*sizeof(float));
        CHECK_FOR_CUDA_ERROR();
    }

    unsigned int num_batches = std::ceil(num_projections / float(orig_ppb));

    if (debugInfo) {
        printf("number of batches: %i\n", num_batches);
        printf("original projections per batch: %i\n", orig_ppb);
        printf("total number then: %i\n", orig_ppb*num_batches);
    }

    for (unsigned int batch = 0; batch < num_batches; batch++) {
        unsigned int from_projection = batch * orig_ppb;
        unsigned int to_projection = (batch+1) * orig_ppb;
        if (from_projection > num_projections)
            from_projection = num_projections;
        if (to_projection > num_projections)
            to_projection = num_projections;

        unsigned int ppb = to_projection-from_projection;

        if (debugInfo) {
            printf("batch: %03i, handling projections: %03i - %03i, angle %.2f - %.2f\n", batch,
                   from_projection, to_projection-1, angles[from_projection], angles[to_projection-1]);
        }

        CHECK_FOR_CUDA_ERROR();
        float* angles_DevPtr;
        cudaMalloc( (void**) &angles_DevPtr, ppb*sizeof(float));
        CHECK_FOR_CUDA_ERROR();

        // Build array for input texture
        cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
        cudaExtent extent;
        extent.width = projection_res_x;
        extent.height = projection_res_y;
        extent.depth = ppb;
        cudaArray *projections_array;
        cudaMalloc3DArray( &projections_array, &channelDesc, extent,cudaArrayLayered );
        CHECK_FOR_CUDA_ERROR();

        unsigned int projSize = projection_res_x*projection_res_y;
        unsigned int pepb = projSize*ppb; // projection elements per batch

        // Allocate device memory
        float *projections_DevPtr;
        cudaMalloc( (void**) &projections_DevPtr, pepb*sizeof(float) );
        CHECK_FOR_CUDA_ERROR();


        //std::cout << "ppb: " << ppb << std::endl;
        for (unsigned int p=from_projection, i=0; p<to_projection; p++, i++) {
            unsigned int from_id = binningdata[p];
            //std::cout << "from: " << from_id << ", to: " << i << std::endl;
            cudaMemcpy(angles_DevPtr+i, &angles[from_id],
                       sizeof(float), cudaMemcpyHostToDevice);
            CHECK_FOR_CUDA_ERROR();

            // Upload projections to device
            cudaMemcpy(projections_DevPtr+i*projSize, projections->get_data_ptr()+from_id*projSize,
                       projSize*sizeof(float), cudaMemcpyHostToDevice);
            CHECK_FOR_CUDA_ERROR();


            if (use_fbp) {
                const float ADD = SDD - SAD; // in mm.
                /*
                  unsigned int num_projection_elements = projection_res_x * projection_res_y;
                */
                  std::vector<unsigned int> dims;
                  dims.push_back(projection_res_x);
                  dims.push_back(projection_res_y);

                // copy to device
                  cuNDArray<TYPE> device_projection_original;
                  device_projection_original.create( &dims, projections_DevPtr+i*projSize );

                // 1. Cosine weighting "SDD / sqrt( SDD*SDD + u*u + v*v)"
                //     elementwize multiply each pixel with weights
                static boost::shared_ptr< cuNDArray<TYPE> > cos_weights;
                static bool initialized = false;
                if (!initialized) {
                    cos_weights = get_conebeam_cosinus_weights(uintd2(projection_res_x, projection_res_y),
                                                      ps_dims_in_mm, ADD, SAD);
                    initialized = true;
                    //write_nd_array<TYPE>( cos_weights->to_host().get(), "cos_weights.real");
                }
                device_projection_original *= *cos_weights;

                cuNDArray<TYPE> device_projection;
                //std::vector<unsigned int> dims = *(host_projection.get_dimensions().get());
                std::vector<unsigned int> double_dims;
                for (unsigned int i=0; i<dims.size(); i++)
                    double_dims.push_back(dims[i] << 1);
                device_projection.create(&double_dims);
                pad<float,2>( &device_projection_original, &device_projection );

                // 2. Ramp filter using FFT
                static boost::shared_ptr< cuNDArray<TYPE> > ramp;
                static bool ramp_initialized = false;
                if (!ramp_initialized) {
                    ramp = get_conebeam_ramp(projection_res_x*2, ps_dims_in_mm[0]/(projection_res_x*2));
                    ramp_initialized = true;
                    //write_nd_array<TYPE>( ramp->to_host().get(), "ramp.real");
                }

                boost::shared_ptr< cuNDArray<complext<TYPE> > > complex_projection =
                    real_to_complex< complext<TYPE> >(&device_projection);
                cuNFFT_plan<float,1>().fft(complex_projection.get(),cuNFFT_plan<float,1>::NFFT_FORWARDS);
                *complex_projection *= *ramp;
                cuNFFT_plan<float,1>().fft(complex_projection.get(),cuNFFT_plan<float,1>::NFFT_BACKWARDS);

                // copy cuNDArray back to hoNDArray
                boost::shared_ptr< cuNDArray<TYPE> > result_original =
                    real(complex_projection.get());

                /*
                boost::shared_ptr< cuNDArray<TYPE> > result(new cuNDArray<TYPE>());
                result->create(&dims);
                */
                uintd<2>::Type offset;
                offset.vec[0]= dims[0]>>1;
                offset.vec[1]= dims[1]>>1;
                crop<float,2>( offset, result_original.get(), &device_projection_original );
            }
            //}
        }

        cudaMemcpy3DParms cpy_params = {0};
        cpy_params.extent = extent;
        cpy_params.dstArray = projections_array;
        cpy_params.kind = cudaMemcpyDeviceToDevice;
        cpy_params.srcPtr =
            make_cudaPitchedPtr( (void*)projections_DevPtr, projection_res_x*sizeof(float),
                                 projection_res_x, projection_res_y );
        cudaMemcpy3D( &cpy_params );
        CHECK_FOR_CUDA_ERROR();

        cudaBindTextureToArray( projection_tex, projections_array, channelDesc );
        CHECK_FOR_CUDA_ERROR();

        // Define dimensions of grid/blocks.
        dim3 dimBlock( matrix_size_x );
        dim3 dimGrid( matrix_size_y, matrix_size_z );

        cudaFuncSetCacheConfig(conebeam_backwards_projection_cb_kernel, cudaFuncCachePreferL1);
        // Invoke kernel
        conebeam_backwards_projection_cb_kernel<<< dimGrid, dimBlock >>>
            (is_dims, ps_dims_in_pixels, ppb, num_projections, x_DevPtr, angles_DevPtr,
             from_projection, ps_dims_in_mm, is_spacing_in_mm,
             sag_parameters_x, sag_parameters_y, SDD, SAD, use_circular_cutoff, use_fbp);
        CHECK_FOR_CUDA_ERROR();

        // Cleanup
        cudaUnbindTexture( projection_tex );
        CHECK_FOR_CUDA_ERROR();
        cudaFreeArray( projections_array );
        CHECK_FOR_CUDA_ERROR();
        cudaFree( projections_DevPtr );
        CHECK_FOR_CUDA_ERROR();
        cudaFree( angles_DevPtr );
        CHECK_FOR_CUDA_ERROR();
    }

    // Copy result from device to host
    unsigned int offs = bin * matrix_size_x * matrix_size_y * matrix_size_z;
    cudaMemcpy( x->get_data_ptr()+offs, x_DevPtr, num_image_elements*sizeof(float), cudaMemcpyDeviceToHost );
    CHECK_FOR_CUDA_ERROR();

    cudaFree( x_DevPtr );
    CHECK_FOR_CUDA_ERROR();

    return true;
}

template bool Gadgetron::conebeam_forwards_projection(hoNDArray<float>*, hoNDArray<float>*, unsigned int bin,
                                           std::vector<unsigned int>& binningdata,
                                           std::vector<float>& angles,
                                           unsigned int orig_ppb,
                                           floatd3 is_spacing_in_mm, floatd2 ps_dims_in_mm,
                                           floatd3 sag_parameters_x, floatd3 sag_parameters_y,
                                           float SDD, float SAD,
                                           unsigned int numSamplesPerRay, bool use_circular_cutoff,
                                           bool accumulate);

template bool Gadgetron::conebeam_backwards_projection(hoNDArray<float>*, hoNDArray<float>*, unsigned int bin,
																						std::vector<unsigned int>& binningdata,
																						std::vector<float>& angles,
                                            unsigned int orig_ppb,
                                            floatd3 is_spacing_in_mm, floatd2 ps_dims_in_mm,
                                            floatd3 sag_parameters_x, floatd3 sag_parameters_y,
                                            float SDD, float SAD, bool use_circular_cutoff,
                                            bool use_fbp, bool accumulate);
