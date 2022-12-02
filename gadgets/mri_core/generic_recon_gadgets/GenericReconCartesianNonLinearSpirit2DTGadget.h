/** \file   GenericReconCartesianNonLinearSpirit2DTGadget.h
    \brief  This is the class gadget for 2DT cartesian non-linear Spirit reconstruction, working on the IsmrmrdReconData.
            The redundant dimension is dimension N.

            Ref: 
            Equation 3 and 4. [1] Hui Xue, Souheil Inati, Thomas Sangild Sorensen, Peter Kellman, Michael S. Hansen. 
            Distributed MRI Reconstruction using Gadgetron based Cloud Computing. Magenetic Resonance in Medicine, 73(3):1015-25, 2015.

    \author Hui Xue
*/

#pragma once

#include "GenericReconCartesianSpiritGadget.h"

namespace Gadgetron {

    class EXPORTGADGETSMRICORE GenericReconCartesianNonLinearSpirit2DTGadget : public GenericReconCartesianSpiritGadget
    {
    public:
        GADGET_DECLARE(GenericReconCartesianNonLinearSpirit2DTGadget);

        typedef GenericReconCartesianSpiritGadget BaseClass;
        typedef Gadgetron::GenericReconCartesianSpiritObj< std::complex<float> > ReconObjType;

        GenericReconCartesianNonLinearSpirit2DTGadget();
        ~GenericReconCartesianNonLinearSpirit2DTGadget();

        /// parameters for workflow
        // GADGET_PROPERTY(spirit_perform_linear                , bool,    "Spirit whether to perform linear reconstruction as the initialization", true);
        /// parameters for regularization
        GADGET_PROPERTY(spirit_parallel_imaging_lamda        , double,  "Spirit regularization strength for parallel imaging term", 1.0);
        GADGET_PROPERTY(spirit_image_reg_lamda               , double,  "Spirit regularization strength for imaging term", 0);
        GADGET_PROPERTY(spirit_data_fidelity_lamda           , double,  "Spirit regularization strength for data fidelity term", 1.0);
        /// parameters for non-linear iteration
        GADGET_PROPERTY(spirit_nl_iter_max                   , int,     "Spirit maximal number of iterations for nonlinear optimization", 0);
        GADGET_PROPERTY(spirit_nl_iter_thres                 , double,  "Spirit threshold to stop iteration for nonlinear optimization", 0);
        /// parameters for image domain regularization, wavelet type regularizer is used here
        GADGET_PROPERTY_LIMITS(spirit_reg_name               , std::string, "Spirit image domain regularizer", "db1", GadgetPropertyLimitsEnumeration, "db1", "db2", "db3", "db4", "db5");
        GADGET_PROPERTY(spirit_reg_level                     , int,     "Spirit image domain regularizer, number of transformation levels", 1);
        GADGET_PROPERTY(spirit_reg_keep_approx_coeff         , bool,    "Spirit whether to keep the approximation coefficients from being regularized", true);
        GADGET_PROPERTY(spirit_reg_keep_redundant_dimension_coeff, bool,    "Spirit whether to keep the boundary coefficients of N dimension from being regularized", false);
        GADGET_PROPERTY(spirit_reg_proximity_across_cha      , bool,    "Spirit whether to perform proximity operation across channels", false);
        GADGET_PROPERTY(spirit_reg_use_coil_sen_map          , bool,    "Spirit whether to use coil map in the imaging term", false);
        GADGET_PROPERTY(spirit_reg_estimate_noise_floor      , bool,    "Spirit whether to estimate noise floor for the imaging term", false);
        GADGET_PROPERTY(spirit_reg_minimal_num_images_for_noise_floor, int,    "Spirit minimal number of images for noise floor estimation", 16);
        /// W matrix of equation 3 and 4 in ref [1]
        GADGET_PROPERTY(spirit_reg_RO_weighting_ratio        , double,  "Spirit regularization weigthing ratio for RO", 1.0);
        GADGET_PROPERTY(spirit_reg_E1_weighting_ratio        , double,  "Spirit regularization weigthing ratio for E1", 1.0);
        GADGET_PROPERTY(spirit_reg_N_weighting_ratio         , double,  "Spirit regularization weigthing ratio for N", 0);

    protected:

        // --------------------------------------------------
        // variable for recon
        // --------------------------------------------------

        // --------------------------------------------------
        // gadget functions
        // --------------------------------------------------
        // default interface function
        virtual int process_config(ACE_Message_Block* mb);

        // --------------------------------------------------
        // recon step functions
        // --------------------------------------------------

        // unwrapping
        virtual void perform_unwrapping(IsmrmrdReconBit& recon_bit, ReconObjType& recon_obj, size_t encoding);

        // perform non-linear spirit unwrapping
        // kspace, kerIm, full_kspace: [RO E1 CHA N S SLC]
        void perform_nonlinear_spirit_unwrapping(hoNDArray< std::complex<float> >& kspace, hoNDArray< std::complex<float> >& kerIm, hoNDArray< std::complex<float> >& ref2DT, hoNDArray< std::complex<float> >& coilMap2DT, hoNDArray< std::complex<float> >& full_kspace, size_t e);
    };
}
