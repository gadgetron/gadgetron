#include "cuNonCartesianMOCOOperator.h"
#include "cubicPrefilter3D.cu"
#include "cubicTex3D.cu"
#include "vector_td_utilities.h"
#include <thrust/extrema.h>
#include "stdio.h"
using namespace Gadgetron;

template <class REAL, unsigned int D>
cuNonCartesianMOCOOperator<REAL, D>::cuNonCartesianMOCOOperator(ConvolutionType conv) : cuSenseOperator<REAL, D>() {

    convolutionType = conv;
    is_preprocessed_ = false;
}

template <class REAL, unsigned int D>
void cuNonCartesianMOCOOperator<REAL, D>::mult_M(cuNDArray<complext<REAL>>* in, cuNDArray<complext<REAL>>* out,
                                                 bool accumulate) {

    if (!in || !out) {
        throw std::runtime_error("cuNonCartesianSenseOperator::mult_M : 0x0 input/output not accepted");
    }
    if (!in->dimensions_equal(&this->domain_dims_) || !out->dimensions_equal(&this->codomain_dims_)) {
        throw std::runtime_error(
            "cuNonCartesianSenseOperator::mult_H: input/output arrays do not match specified domain/codomains");
    }
    // Cart -> noncart
    std::vector<size_t> full_dimensions = *this->get_domain_dimensions();   // cart
    std::vector<size_t> data_dimensions = *this->get_codomain_dimensions(); // Non-cart
    data_dimensions.pop_back();                                             // remove coil dimension from tmp_data;

    // auto timeD = full_dimensions[full_dimensions.size() - 1];
    // full_dimensions.pop_back();
    full_dimensions.push_back(this->ncoils_);
    // full_dimensions.push_back(timeD);

    // std::iter_swap(full_dimensions.end(), full_dimensions.end() - 1); // swap the coil dimension and time

    // full_dimensions.pop_back(); // remove time dimension

    std::vector<size_t> slice_dimensions = *this->get_domain_dimensions();
    // slice_dimensions.pop_back(); // remove time
    auto stride = std::accumulate(slice_dimensions.begin(), slice_dimensions.end(), 1,
                                  std::multiplies<size_t>()); // product of X,Y,and Z

    std::vector<size_t> tmp_dims = *this->get_codomain_dimensions();
    auto stride_data = std::accumulate(tmp_dims.begin(), tmp_dims.end() - 1, 1, std::multiplies<size_t>());

    auto input = cuNDArray<complext<REAL>>(slice_dimensions, in->data());

    slice_dimensions.push_back(shots_per_time_.size()); // add time dimenstion
    cuNDArray<complext<REAL>> moving_images(slice_dimensions);
    slice_dimensions.pop_back();
    for (size_t it = 0; it < shots_per_time_.size(); it++) {

        auto inter_acc = std::accumulate(shots_per_time_.begin(), shots_per_time_.begin() + it, size_t(0)) *
                         tmp_dims[0]; // sum of cum sum shots per time

        // auto slice_view_in = cuNDArray<complext<REAL>>(slice_dimensions, moving_images.data() + stride * it);

        // Move the image to moving image
        auto slice_view_in = input;
        // auto adj_deformation = forward_deformation_[it];
        // adj_deformation *= (REAL)-1.0;
        // applyDeformationbSpline(&slice_view_in, backward_deformation_[it]);
        deform_image(&slice_view_in, backward_deformation_[it]);
        cuNDArray<complext<REAL>> tmp(&full_dimensions);
        this->mult_csm(&slice_view_in, &tmp);

        data_dimensions.pop_back();                     // remove interleave
        data_dimensions.push_back(shots_per_time_[it]); // insert correct interleave
        // data_dimensions.push_back(this->ncoils_);       // insert coils again

        // cuNDArray<complext<REAL>> tmp_data(&data_dimensions);

        full_dimensions.pop_back(); // remove ch
        for (size_t iCHA = 0; iCHA < this->ncoils_; iCHA++) {
            auto tmp_view = cuNDArray<complext<REAL>>(full_dimensions, tmp.data() + stride * iCHA);
            auto slice_view_out =
                cuNDArray<complext<REAL>>(data_dimensions, out->data() + inter_acc + stride_data * iCHA);

            if (accumulate) {
                cuNDArray<complext<REAL>> tmp_out(&full_dimensions);
                plan_[it]->compute(tmp_view, tmp_out, &dcw_[it], NFFT_comp_mode::FORWARDS_C2NC);
                slice_view_out += tmp_out;
            } else
                plan_[it]->compute(tmp_view, slice_view_out, &dcw_[it], NFFT_comp_mode::FORWARDS_C2NC);
        }
        full_dimensions.push_back(this->ncoils_);
        // size_t inter_acc = 0;
        // if (it > 0)

        // This is not correct yet ! -- AJ
        // for (size_t iCHA = 0; iCHA < this->ncoils_; iCHA++)
        //     cudaMemcpy(out->get_data_ptr() + inter_acc + stride_data * iCHA,
        //                tmp_data.get_data_ptr() + tmp_data.get_size(0) * tmp_data.get_size(1) * iCHA,
        //                tmp_data.get_size(0) * tmp_data.get_size(1) * sizeof(complext<REAL>), cudaMemcpyDefault);
    }
}

template <class REAL, unsigned int D>
void cuNonCartesianMOCOOperator<REAL, D>::mult_MH(cuNDArray<complext<REAL>>* in, cuNDArray<complext<REAL>>* out,
                                                  bool accumulate) {

    if (!in || !out) {
        throw std::runtime_error("cuNonCartesianSenseOperator::mult_MH : 0x0 input/output not accepted");
    }

    if (!in->dimensions_equal(&this->codomain_dims_) || !out->dimensions_equal(&this->domain_dims_)) {
        throw std::runtime_error(
            "cuNonCartesianSenseOperator::mult_MH: input/output arrays do not match specified domain/codomains");
    }
    std::vector<size_t> out_dimensions = *this->get_domain_dimensions();
    std::vector<size_t> in_dimensions = *this->get_codomain_dimensions();

    auto RO = in->get_size(0);
    auto E1E2 = in->get_size(1);
    auto CHA = in->get_size(2);

    in_dimensions.pop_back(); // Remove CH dimension

    // out_dimensions.pop_back();               // Remove the timeDimension
    out_dimensions.push_back(this->ncoils_); // add coil dimension
    cuNDArray<complext<REAL>> tmp(&out_dimensions);
    out_dimensions.pop_back(); // rm coil dimension
    // cuNDArray<complext<REAL>> tmp_coilCmb(&out_dimensions);

    auto stride_ch = std::accumulate(in_dimensions.begin(), in_dimensions.end(), 1,
                                     std::multiplies<size_t>()); // product of X,Y,and Z

    auto stride_out = std::accumulate(out_dimensions.begin(), out_dimensions.end(), 1,
                                      std::multiplies<size_t>()); // product of X,Y,and Z
    if (!accumulate) {
        clear(out);
    }
    out_dimensions.push_back(shots_per_time_.size()); // add time dimension
    cuNDArray<complext<REAL>> moving_images(out_dimensions);
    out_dimensions.pop_back(); // rm time

    auto output = cuNDArray<complext<REAL>>(out_dimensions, out->data());
    fill<complext<REAL>>(&output, complext<REAL>((REAL)0, (REAL)0));
    for (size_t it = 0; it < shots_per_time_.size(); it++) {

        size_t inter_acc = std::accumulate(shots_per_time_.begin(), shots_per_time_.begin() + it, 0) * in_dimensions[0];
        in_dimensions.pop_back(); // Remove INT dimension
        in_dimensions.push_back(shots_per_time_[it]);

        for (size_t ich = 0; ich < CHA; ich++) {

            auto slice_view = cuNDArray<complext<REAL>>(in_dimensions, in->data() + stride_ch * ich + inter_acc);
            auto out_view_ch = cuNDArray<complext<REAL>>(out_dimensions, tmp.data() + stride_out * ich);

            plan_[it]->compute(slice_view, out_view_ch, &dcw_[it], NFFT_comp_mode::BACKWARDS_NC2C);
        }

        auto slice_view_output = cuNDArray<complext<REAL>>(out_dimensions, moving_images.data() + stride_out * it);

        this->mult_csm_conj_sum(&tmp, &slice_view_output);
        // applyDeformationbSpline(&slice_view_output, forward_deformation_[it]);
        deform_image(&slice_view_output, forward_deformation_[it]);
        output += slice_view_output;
    }
    output /= complext<REAL>((REAL)shots_per_time_.size(), (REAL)0);
}

template <class REAL, unsigned int D>
void cuNonCartesianMOCOOperator<REAL, D>::setup(_uint64d matrix_size, _uint64d matrix_size_os, REAL W) {
    for (auto ii = 0; ii < shots_per_time_.size(); ii++)
        plan_.push_back(NFFT<cuNDArray, REAL, D>::make_plan(matrix_size, matrix_size_os, W, convolutionType));
}

template <class REAL, unsigned int D>
void cuNonCartesianMOCOOperator<REAL, D>::preprocess(std::vector<cuNDArray<_reald>>& trajectory) {
    if (&(*trajectory.begin()) == 0x0) {
        throw std::runtime_error("cuNonCartesianSenseOperator: cannot preprocess 0x0 trajectory.");
    }

    boost::shared_ptr<std::vector<size_t>> domain_dims = this->get_domain_dimensions();
    if (domain_dims.get() == 0x0 || domain_dims->size() == 0) {
        throw std::runtime_error("cuNonCartesianSenseOperator::preprocess : operator domain dimensions not set");
    }
    for (auto ii = 0; ii < shots_per_time_.size(); ii++)
        plan_[ii]->preprocess(trajectory[ii], NFFT_prep_mode::ALL);
    is_preprocessed_ = true;
}

template <class REAL, unsigned int D>
void cuNonCartesianMOCOOperator<REAL, D>::set_dcw(std::vector<cuNDArray<REAL>> dcw) {
    dcw_ = dcw;
}
template <class REAL, unsigned int D>
void cuNonCartesianMOCOOperator<REAL, D>::applyDeformation(cuNDArray<complext<REAL>>* moving_image,
                                                           cuNDArray<REAL> transformation) {

    // Setup solver
    boost::shared_ptr<cuLinearResampleOperator<REAL, D>> R(new cuLinearResampleOperator<REAL, D>());
    cuCKOpticalFlowSolver<REAL, D> CK;
    CK.set_interpolator(R);

    boost::shared_ptr<cuNDArray<REAL>> result = boost::make_shared<cuNDArray<REAL>>(transformation);

    auto mir = *real(moving_image);
    auto mii = *imag(moving_image);

    auto deformed_movingr = *CK.deform(&mir, result);
    auto deformed_movingi = *CK.deform(&mii, result);

    // Gadgetron::cuNDArray<complext<REAL>> deformed_image =
    //     *cureal_imag_to_complex<complext<REAL>>(&deformed_movingr, &deformed_movingi);
    moving_image = cureal_imag_to_complex<complext<REAL>>(&deformed_movingr, &deformed_movingi).get();
    // return deformed_image;
}
template <class REAL, unsigned int D>
void cuNonCartesianMOCOOperator<REAL, D>::applyDeformationbSpline(cuNDArray<complext<REAL>>* moving_image,
                                                                  cuNDArray<REAL> transformation) {
    auto mir = *real(moving_image);
    auto mii = *imag(moving_image);

    auto ho_mir = hoNDArray<float>(*mir.to_host());
    auto ho_mii = hoNDArray<float>(*mii.to_host());

    auto ho_transformation = *transformation.to_host();
    auto dims = *ho_transformation.get_dimensions();
    //        dims.pop_back();
    dims.erase(dims.begin());

    auto transform = Gadgetron::hoNDArray<Gadgetron::vector_td<float, 3>>(dims);
    GDEBUG("ho_transformation.size(): %d \n", ho_transformation.size());
    std::copy_n(ho_transformation.data(), ho_transformation.size(), (float*)transform.data());

    auto defr =
        cuNDArray<REAL>(hoNDArray<REAL>(Gadgetron::Registration::deform_image_bspline<float, 3>(ho_mir, transform)));
    auto defi =
        cuNDArray<REAL>(hoNDArray<REAL>(Gadgetron::Registration::deform_image_bspline<float, 3>(ho_mii, transform)));
    moving_image = cureal_imag_to_complex<complext<REAL>>(&defr, &defi).get();
}
// Simple transformation kernel
template <typename REAL>
__global__ static void deform_imageKernel(REAL* output, const REAL* vector_field, cudaTextureObject_t texObj, int width,
                                          int height, int depth) {

    const int ixo = blockDim.x * blockIdx.x + threadIdx.x;
    const int iyo = blockDim.y * blockIdx.y + threadIdx.y;
    const int izo = blockDim.z * blockIdx.z + threadIdx.z;

    if (ixo < width && iyo < height && izo < depth) {

        const int idx = ixo + iyo * width + izo * width * height;
        const int elements = width * height * depth;
        REAL ux = vector_field[idx] + (0.5f + ixo);
        REAL uy = vector_field[idx + elements] + (0.5f + iyo);
        REAL uz = vector_field[idx + 2 * elements] + (0.5f + izo);
       // printf("ux: %0.2f, vector_field[%d]: %0.2f, ixo: %d \n",ux,idx,vector_field[idx],ixo);
       // printf("output: %d, input: %d",tex3D<REAL>(texObj, ux, uy, uz),
       // tex3D<REAL>(texObj, (0.5f + ixo), (REAL)(0.5f + iyo), (REAL)(0.5f + izo)));
        output[idx] = tex3D<REAL>(texObj, ux,uy,uz);
    }
}
// Simple transformation kernel
template <typename REAL>
__global__ static void deform_imageKernel(REAL* output, const REAL* vector_field,
                                          texture<REAL, 3, cudaReadModeElementType>& texObj, int width, int height,
                                          int depth) {

    const int ixo = blockDim.x * blockIdx.x + threadIdx.x;
    const int iyo = blockDim.y * blockIdx.y + threadIdx.y;
    const int izo = blockDim.z * blockIdx.z + threadIdx.z;

    if (ixo < width && iyo < height && izo < depth) {

        const int idx = ixo + iyo * width + izo * width * height;
        const int elements = width * height * depth;
        REAL ux = vector_field[idx] + (0.5f + ixo);
        REAL uy = vector_field[idx + elements] + (0.5f + iyo);
        REAL uz = vector_field[idx + 2 * elements] + (0.5f + izo);

        output[idx] = cubicTex3D(texObj, (REAL)ux, (REAL)uy, (REAL)uz);
    }
}
//--------------------------------------------------------------------------
// Declare the typecast CUDA kernels
//--------------------------------------------------------------------------
template <class T> __device__ float Multiplier() { return 1.0f; }
template <> __device__ float Multiplier<uchar>() { return 255.0f; }
template <> __device__ float Multiplier<schar>() { return 127.0f; }
template <> __device__ float Multiplier<ushort>() { return 65535.0f; }
template <> __device__ float Multiplier<short>() { return 32767.0f; }
template <class T> __global__ void CopyCast(uchar* destination, const T* source, uint pitch, uint width) {
    uint2 index =
        make_uint2(__umul24(blockIdx.x, blockDim.x) + threadIdx.x, __umul24(blockIdx.y, blockDim.y) + threadIdx.y);

    float* dest = (float*)(destination + index.y * pitch) + index.x;
    *dest = (1.0f / Multiplier<T>()) * (float)(source[index.y * width + index.x]);
}
//--------------------------------------------------------------------------
// Declare the typecast templated function
// This function can be called directly in C++ programs
//--------------------------------------------------------------------------

//! Allocate GPU memory and copy a voxel volume from CPU to GPU memory
//! and cast it to the normalized floating point format
//! @return the pointer to the GPU copy of the voxel volume
//! @param host  pointer to the voxel volume in CPU (host) memory
//! @param width   volume width in number of voxels
//! @param height  volume height in number of voxels
//! @param depth   volume depth in number of voxels
template <class T> extern cudaPitchedPtr CastVolumeHostToDevice(const T* host, uint width, uint height, uint depth) {
    cudaPitchedPtr device = {0};
    const cudaExtent extent = make_cudaExtent(width * sizeof(float), height, depth);
    cudaMalloc3D(&device, extent);
    const size_t pitchedBytesPerSlice = device.pitch * device.ysize;

    T* temp = 0;
    const uint voxelsPerSlice = width * height;
    const size_t nrOfBytesTemp = voxelsPerSlice * sizeof(T);
    cudaMalloc((void**)&temp, nrOfBytesTemp);

    uint dimX = min(PowTwoDivider(width), 64);
    dim3 dimBlock(dimX, min(PowTwoDivider(height), 512 / dimX));
    dim3 dimGrid(width / dimBlock.x, height / dimBlock.y);
    size_t offsetHost = 0;
    size_t offsetDevice = 0;

    for (uint slice = 0; slice < depth; slice++) {
        cudaMemcpy(temp, host + offsetHost, nrOfBytesTemp, cudaMemcpyHostToDevice);
        CopyCast<T><<<dimGrid, dimBlock>>>((uchar*)device.ptr + offsetDevice, temp, (uint)device.pitch, width);
        offsetHost += voxelsPerSlice;
        offsetDevice += pitchedBytesPerSlice;
    }

    cudaFree(temp); // free the temp GPU volume
    return device;
}
//! Copy a voxel volume from a pitched pointer to a texture
//! @param tex      [output]  pointer to the texture
//! @param texArray [output]  pointer to the texArray
//! @param volume   [input]   pointer to the the pitched voxel volume
//! @param extent   [input]   size (width, height, depth) of the voxel volume
//! @param onDevice [input]   boolean to indicate whether the voxel volume resides in GPU (true) or CPU (false) memory
//! @note When the texArray is not yet allocated, this function will allocate it
template <class T, enum cudaTextureReadMode mode>
void CreateTextureFromVolume(texture<T, 3, mode>* tex, cudaArray** texArray, const cudaPitchedPtr volume,
                             cudaExtent extent, bool onDevice) {

    cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<T>();
    if (*texArray == 0)
        cudaMalloc3DArray(texArray, &channelDesc, extent);
    // copy data to 3D array
    cudaMemcpy3DParms p = {0};
    p.extent = extent;
    p.srcPtr = volume;
    p.dstArray = *texArray;
    p.kind = onDevice ? cudaMemcpyDeviceToDevice : cudaMemcpyHostToDevice;
    cudaMemcpy3D(&p);
    // bind array to 3D texture
    cudaBindTextureToArray(*tex, *texArray, channelDesc);
    tex->normalized = false; // access with absolute texture coordinates
    tex->filterMode = cudaFilterModeLinear;
}
template <class REAL, unsigned int D>
cuNDArray<complext<REAL>> cuNonCartesianMOCOOperator<REAL, D>::deform_image(cuNDArray<complext<REAL>>* image, cuNDArray<REAL> vector_field) {

    cuNDArray<REAL> mir = *real(image);
    cuNDArray<REAL> mii = *imag(image);

    // clear(&output);

    cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<REAL>();
    cudaExtent extent;
    extent.width = image->get_size(0);
    extent.height = image->get_size(1);
    extent.depth = image->get_size(2);

    cudaMemcpy3DParms cpy_params_r = {0};
    cpy_params_r.kind = cudaMemcpyDeviceToDevice;
    cpy_params_r.extent = extent;

    cudaArray* image_array_r;
    cudaMalloc3DArray(&image_array_r, &channelDesc, extent);
    cpy_params_r.dstArray = image_array_r;
    cpy_params_r.srcPtr =
        make_cudaPitchedPtr((void*)mir.data(), extent.width * sizeof(float), extent.width, extent.height);

    cudaMemcpy3DParms cpy_params_i = {0};
    channelDesc = cudaCreateChannelDesc<REAL>();

    cudaArray* image_array_i;
    cudaMalloc3DArray(&image_array_i, &channelDesc, extent);
    cpy_params_i.kind = cudaMemcpyDeviceToDevice;
    cpy_params_i.extent = extent;
    cpy_params_i.dstArray = image_array_i;
    cpy_params_i.srcPtr =
        make_cudaPitchedPtr((void*)mii.data(), extent.width * sizeof(float), extent.width, extent.height);

    // CubicBSplinePrefilter3DTimer((float*)cpy_params_r.srcPtr.ptr, (uint)cpy_params_r.srcPtr.pitch, extent.width,
    //                              extent.height, extent.depth);
    // CubicBSplinePrefilter3DTimer((float*)cpy_params_i.srcPtr.ptr, (uint)cpy_params_i.srcPtr.pitch, extent.width,
    //                              extent.height, extent.depth);
    cudaMemcpy3D(&cpy_params_r);
    cudaMemcpy3D(&cpy_params_i);

    // create the b-spline coefficients texture
    // CreateTextureFromVolume(&coeffs, &coeffArray, cpy_params.srcPtr, volumeExtent, true);
    // cudaDestroyTextureObject(texObj);

    struct cudaResourceDesc resDesc;
    memset(&resDesc, 0, sizeof(resDesc));
    resDesc.resType = cudaResourceTypeArray;
    resDesc.res.array.array = image_array_r;

    struct cudaTextureDesc texDesc;
    memset(&texDesc, 0, sizeof(texDesc));
    texDesc.addressMode[0] = cudaAddressModeClamp;
    texDesc.addressMode[1] = cudaAddressModeClamp;
    texDesc.addressMode[2] = cudaAddressModeClamp;
    texDesc.filterMode = cudaFilterModeLinear;
    texDesc.readMode = cudaReadModeElementType;
    texDesc.normalizedCoords = 0;
    cudaTextureObject_t texObj_r = 0;
    cudaTextureObject_t texObj_i = 0;
    cudaCreateTextureObject(&texObj_r, &resDesc, &texDesc, NULL);
    resDesc.res.array.array = image_array_i;
    cudaCreateTextureObject(&texObj_i, &resDesc, &texDesc, NULL);

    // cudaPitchedPtr bsplineCoeffs = make_cudaPitchedPtr((void*)mii.data(), extent.width * sizeof(float), extent.width,
    // extent.height);

    // create the b-spline coefficients texture
    // cudaArray* coeffArray = 0;
    // texture<float, 3, cudaReadModeElementType> coeffs; // 3D texture
    // cudaExtent volumeExtent = make_cudaExtent(extent.width, extent.height, extent.depth);
    // CreateTextureFromVolume(&coeffs, &coeffArray, cpy_params.srcPtr, volumeExtent, true);
    // cudaFree(bsplineCoeffs.ptr);  //they are now in the coeffs texture, we do not need this anymore

    dim3 threads(8, 8, 8);

    dim3 grid((extent.width + threads.x - 1) / threads.x, (extent.height + threads.y - 1) / threads.y,
              (extent.depth + threads.z - 1) / threads.z);

    cuNDArray<REAL> outputr(image->get_dimensions());
    cuNDArray<REAL> outputi(image->get_dimensions());

    deform_imageKernel<REAL>
        <<<grid, threads>>>(outputr.data(), vector_field.data(), texObj_r, extent.width, extent.height, extent.depth);

    // deform_imageKernel<REAL>
    //     <<<grid, threads>>>(outputr.data(), vector_field.data(), coeffs, extent.width, extent.height, extent.depth);

    // cudaFreeArray(image_array_r);

    cudaDeviceSynchronize();

    deform_imageKernel<REAL>
        <<<grid, threads>>>(outputi.data(), vector_field.data(), texObj_i, extent.width, extent.height, extent.depth);

    cudaDeviceSynchronize();
    // cudaFree(&coeffs);
    // Free device memory

    auto output = *cureal_imag_to_complex<float_complext>(&outputr, &outputi);
    return output;
}

template <class REAL, unsigned int D>
void cuNonCartesianMOCOOperator<REAL, D>::set_shots_per_time(std::vector<size_t> shots_per_time) {
    shots_per_time_ = shots_per_time;
}
template <class REAL, unsigned int D>
void cuNonCartesianMOCOOperator<REAL, D>::set_forward_deformation(std::vector<cuNDArray<REAL>> forward_deformation) {
    forward_deformation_ = forward_deformation;
}
template <class REAL, unsigned int D>
void cuNonCartesianMOCOOperator<REAL, D>::set_backward_deformation(std::vector<cuNDArray<REAL>> backward_deformation) {
    backward_deformation_ = backward_deformation;
}

template class EXPORTGPUPMRI cuNonCartesianMOCOOperator<float, 1>;
template class EXPORTGPUPMRI cuNonCartesianMOCOOperator<float, 2>;
template class EXPORTGPUPMRI cuNonCartesianMOCOOperator<float, 3>;
template class EXPORTGPUPMRI cuNonCartesianMOCOOperator<float, 4>;

// template class EXPORTGPUPMRI cuNonCartesianMOCOOperator<double, 1>;
// template class EXPORTGPUPMRI cuNonCartesianMOCOOperator<double, 2>;
// template class EXPORTGPUPMRI cuNonCartesianMOCOOperator<double, 3>;
// template class EXPORTGPUPMRI cuNonCartesianMOCOOperator<double, 4>;
