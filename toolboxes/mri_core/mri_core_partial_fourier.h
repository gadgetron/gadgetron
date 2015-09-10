
/** \file   mri_core_partial_fourier.h
    \brief  Implementation partial fourier handling functionalities for 2D and 3D MRI
    \author Hui Xue
*/

#pragma once

#include "mri_core_export.h"
#include "hoNDArray.h"

namespace Gadgetron
{
    // --------------------------------------------------------------------------
    /// define the partial fourier handling of ISMRMRD
    // --------------------------------------------------------------------------
    enum ismrmrdPARTIALFOURIERHANDLINGALGO
    {
        ISMRMRD_PF_None,
        ISMRMRD_PF_Zerofilling_filter,
        ISMRMRD_PF_Pocs,
        ISMRMRD_PF_FengHuang
    };

    EXPORTMRICORE std::string get_ismrmrd_partial_fourier_handling_algo_name(ismrmrdPARTIALFOURIERHANDLINGALGO v);
    EXPORTMRICORE ismrmrdPARTIALFOURIERHANDLINGALGO get_ismrmrd_partial_fourier_handling_algo(const std::string& name);

    
    // --------------------------------------------------------------------------
    /// define the partial fourier handler for cartesian sampling
    // --------------------------------------------------------------------------
    template <typename T>
    class EXPORTMRICORE partialFourierCartesianHandler
    {
    public:

        partialFourierCartesianHandler();
        virtual ~partialFourierCartesianHandler();

        /// complexIm: [RO E1 E2 ...]
        /// if E2==1, 2D imaging will be assumed
        virtual void partial_fourier(const hoNDArray<T>& kspace, hoNDArray<T>& res) = 0;

        virtual void dump(std::ostream& os) const;

        /// sampled region in kspace
        size_t start_RO_;
        size_t end_RO_;

        size_t start_E1_;
        size_t end_E1_;

        size_t start_E2_;
        size_t end_E2_;

        /// transition band width in pixel
        /// a transition band can help to remove the ringing caused by imperfect partial foureir computation
        size_t transit_band_RO_;
        size_t transit_band_E1_;
        size_t transit_band_E2_;

        /// whether to print out more information
        bool verbose_;
    };

    /// ------------------------------------------------------------------------
    /// filter handler
    /// perform the asymmetric partial fourier filter
    /// ------------------------------------------------------------------------
    template <typename T>
    class EXPORTMRICORE partialFourierCartesianFilterHandler : public partialFourierCartesianHandler<T>
    {
    public:

        typedef partialFourierCartesianHandler<T> BaseClass;
        typedef typename realType<T>::Type value_type;

        partialFourierCartesianFilterHandler();
        virtual ~partialFourierCartesianFilterHandler();

        virtual void partial_fourier(const hoNDArray<T>& kspace, hoNDArray<T>& res);

        virtual void dump(std::ostream& os) const;

        using BaseClass::start_RO_;
        using BaseClass::end_RO_;
        using BaseClass::start_E1_;
        using BaseClass::end_E1_;
        using BaseClass::start_E2_;
        using BaseClass::end_E2_;
        using BaseClass::verbose_;

        /// PF filter is tapered hanning
        /// whether to perform density compensation
        bool filter_pf_density_comp_;

        /// PF filter tapered width [0 1]
        double filter_pf_width_RO_;
        double filter_pf_width_E1_;
        double filter_pf_width_E2_;

        /// parfial fourier filter
        hoNDArray<T> filter_RO_pf_;
        hoNDArray<T> filter_E1_pf_;
        hoNDArray<T> filter_E2_pf_;
    };

    /// ------------------------------------------------------------------------
    /// POCS
    /// perform the iterative POCS reconstruction
    /// Magnetic Resonance Imaging : Physical Principles and Sequence Design.Page 296 - 297.
    /// E.Mark Haacke, Robert W.Brown, Michael R.Thompson, Ramesh Venkatesan.
    /// Wiley - Liss, ISBN - 10 : 0471351288.
    /// ------------------------------------------------------------------------
    template <typename T>
    class EXPORTMRICORE partialFourierCartesianPOCSHandler : public partialFourierCartesianHandler<T>
    {
    public:

        typedef partialFourierCartesianHandler<T> BaseClass;
        typedef typename realType<T>::Type value_type;

        partialFourierCartesianPOCSHandler();
        virtual ~partialFourierCartesianPOCSHandler();

        virtual void partial_fourier(const hoNDArray<T>& kspace, hoNDArray<T>& res);

        virtual void dump(std::ostream& os) const;

        using BaseClass::start_RO_;
        using BaseClass::end_RO_;
        using BaseClass::start_E1_;
        using BaseClass::end_E1_;
        using BaseClass::start_E2_;
        using BaseClass::end_E2_;
        using BaseClass::transit_band_RO_;
        using BaseClass::transit_band_E1_;
        using BaseClass::transit_band_E2_;
        using BaseClass::verbose_;

        /// number of iterations for POCS
        size_t iter_;
        /// threshold for iteration
        double thres_;
    };

    /// ------------------------------------------------------------------------
    /// Feng Huang method
    /// Feng Huang, Wei Lin, and Yu Li.
    /// Partial Fourier Reconstruction Through Data Fitting and Convolution in k-Space.
    /// Magnetic Resonance in Medicine, Vol 62, page 1261-1269, 2009.
    /// ------------------------------------------------------------------------
    template <typename T>
    class EXPORTMRICORE partialFourierCartesianFengHuangHandler : public partialFourierCartesianHandler<T>
    {
    public:

        typedef partialFourierCartesianHandler<T> BaseClass;
        typedef typename realType<T>::Type value_type;

        partialFourierCartesianFengHuangHandler();
        virtual ~partialFourierCartesianFengHuangHandler();

        virtual void partial_fourier(const hoNDArray<T>& kspace, hoNDArray<T>& res);

        virtual void dump(std::ostream& os) const;

        using BaseClass::start_RO_;
        using BaseClass::end_RO_;
        using BaseClass::start_E1_;
        using BaseClass::end_E1_;
        using BaseClass::start_E2_;
        using BaseClass::end_E2_;
        using BaseClass::transit_band_RO_;
        using BaseClass::transit_band_E1_;
        using BaseClass::transit_band_E2_;
        using BaseClass::verbose_;

        /// kernel size
        size_t kRO_;
        size_t kE1_;
        size_t kE2_;
        /// threshold for calibration regularization
        double thres_;
    };
}
