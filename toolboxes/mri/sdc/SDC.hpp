
#include "SDC.h"

#include "ConvolutionKernel.h"
#include "GriddingConvolution.h"

namespace Gadgetron {
template <class T> struct safe_divides {
    __host__ __device__ T operator()(const T& x, const T& y) const { return y == T(0) ? T(0) : x / y; }
};

template <template <class> class ARRAY, class REAL> struct updates {
    void operator()(const ARRAY<REAL>& src, ARRAY<REAL>& dst);
};

template <template <class> class ARRAY, class REAL, unsigned int D> struct validates {
    vector_td<size_t, D> operator()(const vector_td<size_t, D>& size);
};

template <template <class> class ARRAY, class REAL, unsigned int D>
std::shared_ptr<ARRAY<REAL>> estimate_dcw(const ARRAY<vector_td<REAL, D>>& traj,
                                          const vector_td<size_t, D>& matrix_size, REAL os_factor,
                                          unsigned int num_iterations, REAL kernelWidth) {
    // Initialize weights to 1.
    ARRAY<REAL> dcw(traj.get_dimensions());
    fill(&dcw, (REAL)1);

    // Compute density compensation weights.
    return estimate_dcw(traj, dcw, matrix_size, os_factor, num_iterations, kernelWidth);
}

template <template <class> class ARRAY, class REAL, unsigned int D>
std::shared_ptr<ARRAY<REAL>> estimate_dcw(const ARRAY<vector_td<REAL, D>>& traj, const ARRAY<REAL>& initial_dcw,
                                          const vector_td<size_t, D>& matrix_size, REAL os_factor,
                                          unsigned int num_iterations, REAL kernelWidth) {
    // Specialized functors.
    auto update_weights = updates<ARRAY, REAL>();
    auto validate_size = validates<ARRAY, REAL, D>();

    // Matrix size with oversampling.
    auto matrix_size_os = vector_td<size_t, D>(vector_td<REAL, D>(matrix_size) * os_factor);

    // Validate matrix size.
    auto valid_matrix_size = validate_size(matrix_size);
    auto valid_matrix_size_os = validate_size(matrix_size_os);

    // Convolution kernel.
    auto kernel = JincKernel<REAL, D>(kernelWidth);

    // Prepare gridding convolution.
    auto conv = GriddingConvolution<ARRAY, REAL, D, JincKernel>::make(valid_matrix_size, valid_matrix_size_os, kernel);

    // cudaPointerAttributes attributes;
    // cudaPointerGetAttributes(&attributes,traj.get_data_ptr());


    // if(attributes.devicePointer != NULL)
    //     //conv->initialize(ConvolutionType::ATOMIC);

    conv->preprocess(traj);

    // Working arrays.
    ARRAY<REAL> dcw(initial_dcw);
    ARRAY<REAL> grid(to_std_vector(conv->get_matrix_size_os()));
    ARRAY<REAL> tmp(dcw.get_dimensions());

    // Iteration loop.
    for (size_t i = 0; i < num_iterations; i++) {
        // To intermediate grid.
        conv->compute(dcw, grid, GriddingConvolutionMode::NC2C);

        // To original trajectory.
        conv->compute(grid, tmp, GriddingConvolutionMode::C2NC);

        // Update weights.
        update_weights(tmp, dcw);
    }

    return std::make_shared<ARRAY<REAL>>(dcw);
}



} // namespace Gadgetron