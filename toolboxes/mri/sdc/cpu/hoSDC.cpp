
#include "hoSDC.h"

#include "hoGriddingConvolution.h"
#include "hoNDArray_elemwise.h"

#include "SDC.hpp"

namespace Gadgetron
{
    template<class REAL>
    struct updates<hoNDArray, REAL>
    {
        void operator()(const hoNDArray<REAL>& src, hoNDArray<REAL>& dst)
        {
            std::transform(dst.begin(),
                           dst.end(),
                           src.begin(),
                           dst.begin(),
                           safe_divides<REAL>());
        }
    };

    template<class REAL, unsigned int D>
    struct validates<hoNDArray, REAL, D>
    {
        vector_td<size_t, D> operator()(const vector_td<size_t, D>& size)
        {
            return size;
        }
    };
}


template std::shared_ptr<Gadgetron::hoNDArray<float>>
Gadgetron::estimate_dcw<Gadgetron::hoNDArray, float, 2>(
    const Gadgetron::hoNDArray<Gadgetron::vector_td<float, 2>>& traj,
    const Gadgetron::vector_td<size_t, 2>& matrix_size,
    float os_factor,
    unsigned int num_iterations,
    float kernelWidth);

template std::shared_ptr<Gadgetron::hoNDArray<float>>
Gadgetron::estimate_dcw<Gadgetron::hoNDArray, float, 3>(
    const Gadgetron::hoNDArray<Gadgetron::vector_td<float, 3>>& traj,
    const Gadgetron::vector_td<size_t, 3>& matrix_size,
    float os_factor,
    unsigned int num_iterations,
    float kernelWidth);

template std::shared_ptr<Gadgetron::hoNDArray<float>>
Gadgetron::estimate_dcw<Gadgetron::hoNDArray, float, 2>(
    const Gadgetron::hoNDArray<Gadgetron::vector_td<float, 2>>& traj,
    const Gadgetron::hoNDArray<float>& initial_dcw,
    const Gadgetron::vector_td<size_t, 2>& matrix_size,
    float os_factor,
    unsigned int num_iterations,
    float kernelWidth);

template std::shared_ptr<Gadgetron::hoNDArray<float>>
Gadgetron::estimate_dcw<Gadgetron::hoNDArray, float, 3>(
    const Gadgetron::hoNDArray<Gadgetron::vector_td<float, 3>>& traj,
    const Gadgetron::hoNDArray<float>& initial_dcw,
    const Gadgetron::vector_td<size_t, 3>& matrix_size,
    float os_factor,
    unsigned int num_iterations,
    float kernelWidth);
