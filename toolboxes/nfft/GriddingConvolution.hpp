#include "GriddingConvolution.h"

namespace Gadgetron
{
    template<template<class> class ARRAY, class T, unsigned int D, template<class, unsigned int> class K>
    GriddingConvolutionBase<ARRAY, T, D, K>::GriddingConvolutionBase(
        const vector_td<size_t, D>& matrix_size,
        const vector_td<size_t, D>& matrix_size_os,
        const K<REAL, D>& kernel)
      : matrix_size_(matrix_size)
      , matrix_size_os_(matrix_size_os)
      , kernel_(kernel)
    {
        if (matrix_size > matrix_size_os)
            throw std::runtime_error("Oversampled matrix size must be greater "
                                     "than or equal to matrix size.");
    }

    template<template<class> class ARRAY, class T, unsigned int D, template<class, unsigned int> class K>
    GriddingConvolutionBase<ARRAY, T, D, K>::GriddingConvolutionBase(
        const vector_td <size_t, D>& matrix_size,
        const REAL os_factor,
        const K<REAL, D>& kernel)
      : matrix_size_(matrix_size)
      , matrix_size_os_(vector_td<size_t, D>(os_factor * vector_td<REAL, D>(matrix_size)))
      , kernel_(kernel)
    {
        if (os_factor < REAL(1))
            throw std::runtime_error("Oversampling factor must greater than or "
                                     "equal to 1.");
    }

    template<template<class> class ARRAY, class T, unsigned int D, template<class, unsigned int> class K>
    void GriddingConvolutionBase<ARRAY, T, D, K>::preprocess(
        const ARRAY<vector_td<REAL, D>>& trajectory,
        GriddingConvolutionPrepMode prep_mode)
    {
        this->num_samples_ = trajectory.get_size(0);
        this->num_frames_ = trajectory.get_number_of_elements() / num_samples_;
    }

    template<template<class> class ARRAY, class T, unsigned int D, template<class, unsigned int> class K>
    void GriddingConvolutionBase<ARRAY, T, D, K>::compute(
        const ARRAY<T>& input,
        ARRAY<T>& output,
        GriddingConvolutionMode mode,
        bool accumulate)
    {
        switch (mode)
        {
            case GriddingConvolutionMode::C2NC:
            {
                this->compute_C2NC(input, output, accumulate);
                break;
            }
            case GriddingConvolutionMode::NC2C:
            {
                this->compute_NC2C(input, output, accumulate);
                break;
            }
        }
    }

    template<template<class> class ARRAY, class T, unsigned int D, template<class, unsigned int> class K>
    inline const K<realType_t<T>, D>& GriddingConvolutionBase<ARRAY, T, D, K>::get_kernel() const
    {
        return this->kernel_;
    }

    template<template<class> class ARRAY, class T, unsigned int D, template<class, unsigned int> class K>
    inline vector_td<size_t, D> GriddingConvolutionBase<ARRAY, T, D, K>::get_matrix_size() const
    {
        return this->matrix_size_;
    }

    template<template<class> class ARRAY, class T, unsigned int D, template<class, unsigned int> class K>
    inline vector_td<size_t, D> GriddingConvolutionBase<ARRAY, T, D, K>::get_matrix_size_os() const
    {
        return this->matrix_size_os_;
    }

    template<template<class> class ARRAY, class T, unsigned int D, template<class, unsigned int> class K>
    inline size_t GriddingConvolutionBase<ARRAY, T, D, K>::get_num_samples() const
    {
        return this->num_samples_;
    }

    template<template<class> class ARRAY, class T, unsigned int D, template<class, unsigned int> class K>
    inline size_t GriddingConvolutionBase<ARRAY, T, D, K>::get_num_frames() const
    {
        return this->num_frames_;
    }

} // namespace Gadgetron
