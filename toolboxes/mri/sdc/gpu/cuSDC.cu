
#include "cuSDC.h"

#include "cuGriddingConvolution.h"
#include "cuNDArray_elemwise.h"

#include "mri_sdc_export.h"

#include "SDC.hpp"

namespace Gadgetron
{
    template<class REAL>
    struct updates<cuNDArray, REAL>
    {
        void operator()(const cuNDArray<REAL>& src, cuNDArray<REAL>& dst)
        {
            thrust::transform(dst.begin(),
                              dst.end(),
                              src.begin(),
                              dst.begin(),
                              safe_divides<REAL>());
        }
    };

    template<class REAL, unsigned int D>
    struct validates<cuNDArray, REAL, D>
    {
        vector_td<size_t, D> operator()(const vector_td<size_t, D>& size)
        {
            // Get warp size of current device.
            int device;
            if (cudaGetDevice(&device) != cudaSuccess)
                throw cuda_error("Could not retrieve the active device.");
            vector_td<size_t, D> warp_size(static_cast<size_t>(
                cudaDeviceManager::Instance()->warp_size(device)));
            
            // Enforce that matrix size is a multiple of warp size.
            return ((size + warp_size - size_t(1)) / warp_size) * warp_size;
        }
    };
}


template EXPORTSDC std::shared_ptr<Gadgetron::cuNDArray<float>>
Gadgetron::estimate_dcw<Gadgetron::cuNDArray, float, 2>(
    const Gadgetron::cuNDArray<Gadgetron::vector_td<float, 2>>& traj,
    const Gadgetron::vector_td<size_t, 2>& matrix_size,
    float os_factor,
    unsigned int num_iterations);

template EXPORTSDC std::shared_ptr<Gadgetron::cuNDArray<float>>
Gadgetron::estimate_dcw<Gadgetron::cuNDArray, float, 3>(
    const Gadgetron::cuNDArray<Gadgetron::vector_td<float, 3>>& traj,
    const Gadgetron::vector_td<size_t, 3>& matrix_size,
    float os_factor,
    unsigned int num_iterations);

template EXPORTSDC std::shared_ptr<Gadgetron::cuNDArray<float>>
Gadgetron::estimate_dcw<Gadgetron::cuNDArray, float, 2>(
    const Gadgetron::cuNDArray<Gadgetron::vector_td<float, 2>>& traj,
    const Gadgetron::cuNDArray<float>& initial_dcw,
    const Gadgetron::vector_td<size_t, 2>& matrix_size,
    float os_factor,
    unsigned int num_iterations);

template EXPORTSDC std::shared_ptr<Gadgetron::cuNDArray<float>>
Gadgetron::estimate_dcw<Gadgetron::cuNDArray, float, 3>(
    const Gadgetron::cuNDArray<Gadgetron::vector_td<float, 3>>& traj,
    const Gadgetron::cuNDArray<float>& initial_dcw,
    const Gadgetron::vector_td<size_t, 3>& matrix_size,
    float os_factor,
    unsigned int num_iterations);
    