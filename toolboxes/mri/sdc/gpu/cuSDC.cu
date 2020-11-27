
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
    template <class REAL, unsigned int D>
std::shared_ptr<cuNDArray<REAL>> estimate_dcw(const cuNDArray<vector_td<REAL, D>>& traj,
                                              const vector_td<size_t, D>& matrix_size, REAL os_factor,
                                              unsigned int num_iterations, ConvolutionType convtype) {
    // Initialize weights to 1.
    cuNDArray<REAL> dcw(*traj.get_dimensions());
    fill(&dcw, (REAL)1);

    // Compute density compensation weights.
    return estimate_dcw(traj, dcw, matrix_size, os_factor, num_iterations, convtype);
}

template <class REAL, unsigned int D>
std::shared_ptr<cuNDArray<REAL>> estimate_dcw(const cuNDArray<vector_td<REAL, D>>& traj,
                                              const cuNDArray<REAL>& initial_dcw,
                                              const vector_td<size_t, D>& matrix_size, REAL os_factor,
                                              unsigned int num_iterations, ConvolutionType convtype) {
    // Specialized functors.
    auto update_weights = updates<cuNDArray, REAL>();
    auto validate_size = validates<cuNDArray, REAL, D>();

    // Matrix size with oversampling.
    auto matrix_size_os = vector_td<size_t, D>(vector_td<REAL, D>(matrix_size) * os_factor);

    // Validate matrix size.
    auto valid_matrix_size = validate_size(matrix_size);
    auto valid_matrix_size_os = validate_size(matrix_size_os);

    // Convolution kernel.
    auto kernel = JincKernel<REAL, D>(vector_td<unsigned int, D>(valid_matrix_size),
                                      vector_td<unsigned int, D>(valid_matrix_size_os));

    // Prepare gridding convolution.
    auto conv = GriddingConvolution<cuNDArray, REAL, D, JincKernel>::make(valid_matrix_size, valid_matrix_size_os, kernel,convtype);
    // cudaPointerAttributes attributes;
    // cudaPointerGetAttributes(&attributes,traj.get_data_ptr());

    // if(attributes.devicePointer != NULL)
    //     //conv->initialize(ConvolutionType::ATOMIC);

    conv->preprocess(traj,GriddingConvolutionPrepMode::ALL);

    // Working arrays.
    cuNDArray<REAL> dcw(initial_dcw);
    cuNDArray<REAL> grid(to_std_vector(conv->get_matrix_size_os()));
    cuNDArray<REAL> tmp(*dcw.get_dimensions());

    // Iteration loop.
    for (size_t i = 0; i < num_iterations; i++) {
               
        // To intermediate grid.
        conv->compute(dcw, grid, GriddingConvolutionMode::NC2C);
        
         // To original trajectory.
        conv->compute(grid, tmp, GriddingConvolutionMode::C2NC);



        // Update weights.
        update_weights(tmp, dcw);
    }

    return std::make_shared<cuNDArray<REAL>>(dcw);
}
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

template EXPORTSDC std::shared_ptr<Gadgetron::cuNDArray<float>>
Gadgetron::estimate_dcw<float, 3>(
    const Gadgetron::cuNDArray<Gadgetron::vector_td<float, 3>>& traj,
    const Gadgetron::vector_td<size_t, 3>& matrix_size,
    float os_factor,
    unsigned int num_iterations, 
    ConvolutionType convtype);

template EXPORTSDC std::shared_ptr<Gadgetron::cuNDArray<float>>
Gadgetron::estimate_dcw<float, 2>(
    const Gadgetron::cuNDArray<Gadgetron::vector_td<float, 2>>& traj,
    const Gadgetron::vector_td<size_t, 2>& matrix_size,
    float os_factor,
    unsigned int num_iterations, 
    ConvolutionType convtype);

template EXPORTSDC std::shared_ptr<Gadgetron::cuNDArray<float>>
Gadgetron::estimate_dcw<float, 3>(
    const Gadgetron::cuNDArray<Gadgetron::vector_td<float, 3>>& traj,
    const Gadgetron::cuNDArray<float>& initial_dcw,
    const Gadgetron::vector_td<size_t, 3>& matrix_size,
    float os_factor,
    unsigned int num_iterations, 
    ConvolutionType convtype);

template EXPORTSDC std::shared_ptr<Gadgetron::cuNDArray<float>>
Gadgetron::estimate_dcw<float, 2>(
    const Gadgetron::cuNDArray<Gadgetron::vector_td<float, 2>>& traj,
    const Gadgetron::cuNDArray<float>& initial_dcw,
    const Gadgetron::vector_td<size_t, 2>& matrix_size,
    float os_factor,
    unsigned int num_iterations, 
    ConvolutionType convtype);
    