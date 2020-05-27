/**
 * \file SDC_kernel.h
 * \brief Kernel for sampling density compensation.
 */

#pragma once

#include "core_defines.h"
#include "hoNDArray.h"
#include "log.h"
#include "vector_td.h"
#include "vector_td_operators.h"
#include "vector_td_utilities.h"

#include <cassert>
#include <cmath>
#include <vector>

namespace Gadgetron
{
    namespace SDC_internal
    {
        template<class REAL, unsigned int D>
        struct KernelProperties
        {
            // polynomial order
            static constexpr size_t POLY_ORDER = 5;

            // length of the table used in polyfit
            static constexpr size_t FIT_LENGTH = 394;

            // length of the spectral domain table
            static constexpr size_t SPECTRAL_LENGTH = 25600;

            // length of the FOV (zeta)
            static constexpr size_t FOV_LENGTH = 63;

            // length of lookup table
            static constexpr size_t LUT_LENGTH = 10000;
        
            // radius FOV product (length to the first truncation point [pixels])
            static constexpr REAL RFP = 0.96960938;
        };
    }

    /**
     * @brief Kernel for sampling density compensation.
     * 
     * @tparam ARRAY Array type.
     * @tparam REAL Floating-point type.
     * @tparam D Number of dimensions.
     * @remark The current implementation does not depend on ARRAY,
     * but the template parameter is intentionally left in place in
     * case future development needs to specialize.
     */
    template<template<class> class ARRAY, class REAL, unsigned int D>
    class SDC_kernel
    {
    public:

        // default constructor
        SDC_kernel();

        // main constructor
        SDC_kernel(const vector_td<size_t, D>& matrix_size, const vector_td<size_t, D>& grid_size);

        // compute kernel given normalized radius
        __host__ __device__ static REAL compute(REAL r);

        // compute kernel using polynomial fit
        REAL compute(const vector_td<REAL, D>& u) const;

        // compute kernel given coordinates and RFP
        __host__ __device__ static REAL compute(const vector_td<REAL, D>& u, REAL rfp);

        // look up kernel value
        REAL lookup(REAL r2) const;

        // look up kernel value
        REAL lookup(const vector_td<REAL, D>& u) const;

        // get grid size
        const vector_td<size_t, D>& get_grid_size() const;

        // get kernel radius FOV product
        REAL get_rfp() const;

    protected:
        
        // compute lookup table
        void compute_lut();

        // coordinates to radius squared
        __host__ __device__ static REAL co_to_r2(const vector_td<REAL, D>& u, REAL rfp);

        // size of intermediate grid - this determines kernel size
        vector_td<size_t, D> grid_size_;

        // oversampling factor
        REAL os_factor_;

        // radius FOV product
        REAL rfp_;

        // lookup table
        std::vector<REAL> lut_;
    };

    using namespace SDC_internal;

    template<template<class> class ARRAY, class REAL, unsigned int D>
    SDC_kernel<ARRAY, REAL, D>::SDC_kernel()
    {
        
    }

    template<template<class> class ARRAY, class REAL, unsigned int D>
    SDC_kernel<ARRAY, REAL, D>::SDC_kernel(const vector_td<size_t, D>& matrix_size, const vector_td<size_t, D>& grid_size)
      : grid_size_(grid_size)
    {
        // oversampling factor
        vector_td<REAL, D> os_factor = vector_td<REAL, D>(grid_size_) / vector_td<REAL, D>(matrix_size);
        os_factor_ = max(os_factor);

        // kernel radius
        rfp_ = KernelProperties<REAL, D>::RFP * os_factor_;

        // compute lookup table
        compute_lut();
    }

    template<template<class> class ARRAY, class REAL, unsigned int D>
    REAL SDC_kernel<ARRAY, REAL, D>::compute(REAL r)
    {
        // not valid for r < zero
        assert(r >= REAL(0));

        // scale index to match table
        REAL x = KernelProperties<REAL, D>::SPECTRAL_LENGTH * r / KernelProperties<REAL, D>::FOV_LENGTH;

        // not valid beyond rfp
        if (x >= KernelProperties<REAL, D>::FIT_LENGTH) return REAL(0);

        // define polynomial coefficients
        // poly[0]*x^(POLY_LEN-1) + ... poly[POLY_LEN-1]*1
        REAL poly[] = {
            -1.1469041640943728E-13, 
            8.5313956268885989E-11, 
            1.3282009203652969E-08, 
            -1.7986635886194154E-05, 
            3.4511129626832091E-05, 
            0.99992359966186584 };

        // get the zeroth order coefficient
        REAL out = poly[KernelProperties<REAL, D>::POLY_ORDER]; // x^0

        // add up polynomial for this point
        for (int i = 1; i <= (int)KernelProperties<REAL, D>::POLY_ORDER; i++)
        {
            out += pow(x, i) * poly[KernelProperties<REAL, D>::POLY_ORDER - i];
        }
            
        // clip negative lobes
        if (out < REAL(0)) out = REAL(0);

        return out;
    }

    template<template<class> class ARRAY, class REAL, unsigned int D>
    REAL SDC_kernel<ARRAY, REAL, D>::compute(const vector_td<REAL, D>& u) const
    {
        return compute(co_to_r2(u, rfp_));
    }

    template<template<class> class ARRAY, class REAL, unsigned int D>
    REAL SDC_kernel<ARRAY, REAL, D>::compute(const vector_td<REAL, D>& u, REAL rfp)
    {
        return compute(sqrt(co_to_r2(u, rfp)));
    }

    template<template<class> class ARRAY, class REAL, unsigned int D>
    REAL SDC_kernel<ARRAY, REAL, D>::lookup(REAL r2) const
    {
        if (r2 > REAL(1.0))
        {
            return REAL(0);
        }
        else
        {
            return lut_[static_cast<size_t>(r2 * (lut_.size() - 1) + 0.5)];
        }
    }

    template<template<class> class ARRAY, class REAL, unsigned int D>
    REAL SDC_kernel<ARRAY, REAL, D>::lookup(const vector_td<REAL, D>& u) const
    {
        return lookup(co_to_r2(u, rfp_));
    }

    template<template<class> class ARRAY, class REAL, unsigned int D>
    void SDC_kernel<ARRAY, REAL, D>::compute_lut()
    {
        // load based on radius squared
        REAL rfp2 = KernelProperties<REAL, D>::RFP * KernelProperties<REAL, D>::RFP;

        // length of LUT with oversampling
        size_t length = static_cast<size_t>(KernelProperties<REAL, D>::LUT_LENGTH * os_factor_);

        // allocate memory for table
        lut_.resize(length);

        // fill table
        for(size_t i = 0; i < length; i++)
        {
            lut_[i] = compute(sqrt(rfp2 * static_cast<REAL>(i) / static_cast<REAL>(length - 1)));
        }
    }

    template<template<class> class ARRAY, class REAL, unsigned int D>
    REAL SDC_kernel<ARRAY, REAL, D>::co_to_r2(const vector_td<REAL, D>& u, REAL rfp)
    {
        vector_td<REAL, D> v = u / rfp;
        REAL r2 = REAL(0);
        for (size_t d = 0; d < D; d++)
        {
            r2 += v[d] * v[d];
        }
        return r2;
    }

    template<template<class> class ARRAY, class REAL, unsigned int D>
    const vector_td<size_t, D>& SDC_kernel<ARRAY, REAL, D>::get_grid_size() const
    {
        return grid_size_;
    }

    template<template<class> class ARRAY, class REAL, unsigned int D>
    REAL SDC_kernel<ARRAY, REAL, D>::get_rfp() const
    {
        return rfp_;
    }
    
}   // namespace Gadgetron