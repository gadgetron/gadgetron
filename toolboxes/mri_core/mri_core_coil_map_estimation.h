
/** \file   mri_core_coil_map_estimation.h
    \brief  Implementation MRI coil sensitivity map estimation functions.

    ISMRMRD_SOUHEIL coil map estimation is based on:

        Inati SJ, Hansen MS, Kellman P.
        A solution to the phase problem in adaptive coil combination.
        In: ISMRM proceeding; April; Salt Lake City, Utah, USA; 2013. 2672.

        Kellman P, McVeigh ER.
        Image reconstruction in SNR units: A general method for SNR measurement.
        Magnetic Resonance in Medicine 2005;54(6):1439-1447.

    ISMRMRD_SOUHEIL_ITER coil map estimation is based on:

        Inati SJ, Hansen MS, Kellman P. 
        A Fast Optimal Method for Coil Sensitivity Estimation and Adaptive Coil Combination for Complex Images.
        In: ISMRM proceeding; May; Milan, Italy; 2014. 4407.

    \author Hui Xue
*/

#pragma once

#include "mri_core_export.h"
#include "hoNDArray.h"

namespace Gadgetron
{
    // --------------------------------------------------------------------------
    /// define the coil sensitivity map estimation algorithms
    // --------------------------------------------------------------------------
    enum ismrmrdCOILMAPALGO
    {
        ISMRMRD_Inati,
        ISMRMRD_Inati_iterative
    };

    EXPORTMRICORE std::string get_ismrmrd_coil_map_algo_name(ismrmrdCOILMAPALGO v);
    EXPORTMRICORE ismrmrdCOILMAPALGO get_ismrmrd_coil_map_algo(const std::string& name);

    /// define the coil map maker
    template <typename T> 
    class EXPORTMRICORE coilMapMaker
    {
    public:

        coilMapMaker();
        virtual ~coilMapMaker();

        /// complexIm: [RO E1 E2 CHA ...]
        /// if E2 == 1, the 2D coil map estimation will be assumed
        virtual void make_coil_map(const hoNDArray<T>& complexIm, hoNDArray<T>& coilMap) = 0;

        virtual void dump(std::ostream& os) const;
    };

    template <typename T> 
    class EXPORTMRICORE coilMapMakerInati : public coilMapMaker<T>
    {
    public:

        typedef coilMapMaker<T> BaseClass;

        coilMapMakerInati();
        virtual ~coilMapMakerInati();

        virtual void make_coil_map(const hoNDArray<T>& complexIm, hoNDArray<T>& coilMap);

        virtual void dump(std::ostream& os) const;

        /// parameters

        /// kernel size in the unit of pixel
        size_t ks_;
        /// number of iterations for power method
        size_t power_;
    };

    template <typename T> 
    class EXPORTMRICORE coilMapMakerInatiIter : public coilMapMaker<T>
    {
    public:

        typedef coilMapMaker<T> BaseClass;
        typedef typename realType<T>::Type value_type;

        coilMapMakerInatiIter();
        virtual ~coilMapMakerInatiIter();

        virtual void make_coil_map(const hoNDArray<T>& complexIm, hoNDArray<T>& coilMap);

        virtual void dump(std::ostream& os) const;

        /// parameters

        /// kernel size in the unit of pixel
        size_t ks_;
        /// kernel size in the unit of pixel for E2 dimension
        size_t kz_;
        /// number of iterations
        size_t iter_num_;
        /// threshold for iteration
        value_type thres_;
    };
}
