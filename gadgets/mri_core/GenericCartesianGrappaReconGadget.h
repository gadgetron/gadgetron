/** \file   GenericCartesianGrappaReconGadget.h
    \brief  This is the class gadget for both 2DT and 3DT cartesian grappa and grappaone reconstruction, working on the IsmrmrdReconData.

            This class implements a general cartesian recon chain: 

            Kspace -------------------------------------------------------> (padding for squared pixel recon)
            Ref                                                                           |
            |-> (padding for squared pixel recon)                                         |
                  |-> ref preparation for embedded/separate/interleaved modes             |
                        |-> coil compression                                              |
                              |-> calibration ---------------------------------------------> coil combination -> PF handling (GenericReconPartialFourierHandlingGadget)
                                                                                                                  |-> kspace filter (GenericReconKSpaceFilteringGadget)
                                                                                                                        |-> FOV adjustment and image resizing (GenericReconFieldOfViewAdjustmentGadget)
                                                                                                                             |-> Image array scaling (ImageArrayScalingGadget)

            The input is IsmrmrdReconData and output is single 2D or 3D ISMRMRD images

            If required, the gfactor map can be sent out

            If the  number of required destination channel is 1, the GrappaONE recon will be performed

            The image number computation logic is implemented in compute_image_number function, which can be overloaded

    \author Hui Xue
*/

#pragma once

#include <complex>
#include "gadgetron_mricore_export.h"
#include "Gadget.h"
#include "hoNDArray.h"
#include "ismrmrd/ismrmrd.h"
#include "ismrmrd/xml.h"
#include "ismrmrd/meta.h"
#include "GadgetronTimer.h"

#include "hoNDArray_utils.h"

#include "gtPlusIOAnalyze.h"

#include "GadgetStreamController.h"

#include "hoNDArray_utils.h"
#include "hoNDArray_elemwise.h"
#include "hoNDFFT.h"

#include "mri_core_data.h"
#include "mri_core_utility.h"
#include "mri_core_coil_map_estimation.h"

namespace Gadgetron {

    /// define the recon status of GenericCartesianGrappaRecon
    template <typename T>
    class EXPORTGADGETSMRICORE GenericCartesianGrappaReconObj
    {
    public:

        GenericCartesianGrappaReconObj();
        GenericCartesianGrappaReconObj(const GenericCartesianGrappaReconObj<T>& v);
        virtual ~GenericCartesianGrappaReconObj();

        GenericCartesianGrappaReconObj<T>& operator=(const GenericCartesianGrappaReconObj<T>& v);

        // ------------------------------------
        /// recon outputs
        // ------------------------------------
        /// reconstructed images, headers and meta attributes
        IsmrmrdImageArray recon_res_;

        /// gfactor, [RO E1 E2 uncombinedCHA+1 N S SLC]
        hoNDArray<typename realType<T>::Type> gfactor_;

        // ------------------------------------
        /// buffers used in the recon
        // ------------------------------------
        /// [RO E1 E2 CHA Nor1 Sor1 SLC]
        hoNDArray<T> ref_calib_;

        /// reference data ready for coil map estimation
        /// [RO E1 E2 CHA Nor1 Sor1 SLC]
        hoNDArray<T> ref_coil_map_;

        /// ref calibration array in dst channels
        hoNDArray<T> ref_calib_dst_;

        /// kspace data in dst channel
        hoNDArray<T> data_dst_;

        /// ref coil map array in dst channel
        hoNDArray<T> ref_coil_map_dst_;

        /// for combined imgae channel
        /// convolution kernel, [RO E1 E2 srcCHA - uncombinedCHA dstCHA - uncombinedCHA Nor1 Sor1 SLC]
        hoNDArray<T> kernel_;
        /// image domain kernel, [RO E1 E2 srcCHA - uncombinedCHA dstCHA - uncombinedCHA Nor1 Sor1 SLC]
        hoNDArray<T> kernelIm_;
        /// image domain unmixing coefficients, [RO E1 E2 srcCHA - uncombinedCHA Nor1 Sor1 SLC]
        hoNDArray<T> unmixing_coeff_;

        /// coil sensitivity map, [RO E1 E2 dstCHA - uncombinedCHA Nor1 Sor1 SLC]
        hoNDArray<T> coil_map_;
    };
}

namespace Gadgetron {

    class EXPORTGADGETSMRICORE GenericCartesianGrappaReconGadget : public Gadget1<IsmrmrdReconData>
    {
    public:
        GADGET_DECLARE(GenericCartesianGrappaReconGadget);

        typedef Gadget1<IsmrmrdReconData> BaseClass;
        typedef Gadgetron::GenericCartesianGrappaReconObj< std::complex<float> > ReconObjType;

        GenericCartesianGrappaReconGadget();
        ~GenericCartesianGrappaReconGadget();

        /// ------------------------------------------------------------------------------------
        /// parameters to control the reconstruction
        /// ------------------------------------------------------------------------------------

        /// image series
        GADGET_PROPERTY(image_series, int, "Image series number", 0);

        /// ------------------------------------------------------------------------------------
        /// debug and timing
        GADGET_PROPERTY(verbose, bool, "Whether to print more information", false);
        GADGET_PROPERTY(debug_folder, std::string, "If set, the debug output will be written out", "");
        GADGET_PROPERTY(perform_timing, bool, "Whether to perform timing on some computational steps", false);

        /// ------------------------------------------------------------------------------------
        /// image sending
        GADGET_PROPERTY(send_out_gfactor, bool, "Whether to send out gfactor map", false);
        GADGET_PROPERTY(send_out_snr_map, bool, "Whether to send out SNR map", false);

        /// ------------------------------------------------------------------------------------
        /// ref preparation
        /// whether to average all N for ref generation
        GADGET_PROPERTY(average_all_ref_N, bool, "Whether to average all N for ref generation", true);
        /// whether to average all S for ref generation
        GADGET_PROPERTY(average_all_ref_S, bool, "Whether to average all S for ref generation", false);

        /// ------------------------------------------------------------------------------------
        /// Grappa parameters
        GADGET_PROPERTY(grappa_kSize_RO, int, "Grappa kernel size RO", 5);
        GADGET_PROPERTY(grappa_kSize_E1, int, "Grappa kernel size E1", 4);
        GADGET_PROPERTY(grappa_kSize_E2, int, "Grappa kernel size E2", 4);
        GADGET_PROPERTY(grappa_reg_lamda, double, "Grappa regularization threshold", 0.0005);
        GADGET_PROPERTY(grappa_calib_over_determine_ratio, double, "Grappa calibration overdermination ratio", 45);

    protected:

        // --------------------------------------------------
        // variables for protocol
        // --------------------------------------------------

        // number of encoding spaces in the protocol
        size_t num_encoding_spaces_;

        // for every encoding space

        // acceleration factor for E1 and E2
        std::vector<double> acceFactorE1_;
        std::vector<double> acceFactorE2_;

        // calibration mode
        std::vector<Gadgetron::ismrmrdCALIBMODE> calib_mode_;

        // encoding space limits
        std::vector<ISMRMRD::EncodingCounters> meas_max_idx_;

        // --------------------------------------------------
        // variable for recon
        // --------------------------------------------------
        // record the recon kernel, coil maps etc. for every encoding space
        std::vector< ReconObjType > recon_obj_;

        // number of times the process function is called
        size_t process_called_times_;

        /// buffers used during recon
        hoNDArray< std::complex<float> > complex_im_recon_buf_;
        hoNDArray< std::complex<float> > data_recon_buf_;

        // filter used for ref coil map 
        hoNDArray< std::complex<float> > filter_RO_ref_coi_map_;
        hoNDArray< std::complex<float> > filter_E1_ref_coi_map_;
        hoNDArray< std::complex<float> > filter_E2_ref_coi_map_;

        // --------------------------------------------------
        // variables for debug and timing
        // --------------------------------------------------

        // debug folder
        std::string debug_folder_full_path_;

        // clock for timing
        Gadgetron::GadgetronTimer gt_timer1_;
        Gadgetron::GadgetronTimer gt_timer2_;
        Gadgetron::GadgetronTimer gt_timer3_;

        // exporter
        Gadgetron::gtPlus::gtPlusIOAnalyze gt_exporter_;

        // in verbose mode, more info is printed out
        bool verbose_;

        // --------------------------------------------------
        // gadget functions
        // --------------------------------------------------
        // default interface function
        virtual int process_config(ACE_Message_Block* mb);
        virtual int process(Gadgetron::GadgetContainerMessage< IsmrmrdReconData >* m1);

        // close call
        int close(unsigned long flags);

        // --------------------------------------------------
        // recon step functions
        // --------------------------------------------------

        // recon chain
        // prepare recon
        // if the squared pixel recon if prescribed, the new recon image size to ensure squared pixel will be computed
        // if needed, the incoming data matrix will be padded and ref_preparer_ recon size will be updated
        // then all recon components will be set with updated parameters
        // finally, the consistency of this recon will be checked; if necessary, changes will be made to ensure the parameter consistency
        virtual void prepare_recon(IsmrmrdReconBit& recon_bit, ReconObjType& recon_obj, size_t encoding);

        // prepare the ref data
        // this may be the most complicated part of the whole recon
        // for every calibration mode, the ref is prepared accordingly
        virtual void prepare_ref(IsmrmrdReconBit& recon_bit, ReconObjType& recon_obj, size_t encoding);

        // generate the destination channel
        virtual void generate_downstream_dst_channel(IsmrmrdReconBit& recon_bit, ReconObjType& recon_obj, size_t encoding);

        // estimate coil map
        virtual void perform_coil_map_estimation(IsmrmrdReconBit& recon_bit, ReconObjType& recon_obj, size_t encoding);

        /// complexIm: [RO E1 E2 CHA ...]
        /// if E2 == 1, the 2D coil map estimation will be assumed
        virtual void compute_coil_map(const hoNDArray< std::complex<float> >& complexIm, hoNDArray< std::complex<float> >& coilMap);

        // calibration, if only one dst channel is prescribed, the GrappaOne is used
        virtual void perform_calib(IsmrmrdReconBit& recon_bit, ReconObjType& recon_obj, size_t encoding);

        // unwrapping or coil combination
        virtual void perform_unwrapping(IsmrmrdReconBit& recon_bit, ReconObjType& recon_obj, size_t encoding);

        // compute image header
        virtual void compute_image_header(IsmrmrdReconBit& recon_bit, ReconObjType& recon_obj, size_t encoding);

        // send out the recon results
        virtual int send_out_image_array(IsmrmrdReconBit& recon_bit, IsmrmrdImageArray& res, size_t encoding, int series_num, const std::string& data_role);

        // --------------------------------------------------
        // implementation functions
        // --------------------------------------------------

        // perform one calibration
        virtual void perform_calib_impl(size_t n, size_t s, size_t slc, size_t e, const hoNDArray< std::complex<float> >& src, const hoNDArray< std::complex<float> >& dst, ReconObjType& recon_obj);

        // perform coil combination with uncombined channels
        // complex_im: [RO E1 E2 CHA N S SLC]
        // coil_map: [RO E1 E2 CHA-unCHA Nor1 Sor1 SLC]
        // res: [RO E1 E2 unCHA+1 N S SLC]
        void coil_combination(const hoNDArray< std::complex<float> >& complex_im, const hoNDArray< std::complex<float> >& coil_map, hoNDArray< std::complex<float> >& res);

        // --------------------------------------------------
        // utility functions
        // --------------------------------------------------
        void get_current_moment(std::string& procTime);
        void get_debug_folder_path(const std::string& debugFolder, std::string& debugFolderPath);

        // compute image number
        virtual size_t compute_image_number(ISMRMRD::ImageHeader& imheader, size_t encoding = 0, size_t CHA = 1, size_t cha = 0, size_t E2 = 1);
    };
}
