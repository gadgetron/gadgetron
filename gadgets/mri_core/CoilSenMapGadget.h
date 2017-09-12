/** \file   CoilSenMapGadget.h
    \brief  This gadget is a part of the depedency handling chain.

    This gadget will accumulate all received coil sensitivity data and buffer them into the gadgetron folder.
    If the incoming data has two SETs, the second SET is considered as the body coil readouts

    \author     Hui Xue
*/

#pragma once

#include "Gadget.h"
#include "hoNDArray.h"
#include "gadgetron_mricore_export.h"
#include "GadgetronTimer.h"

#include <ismrmrd/ismrmrd.h>
#include <ismrmrd/xml.h>
#include <complex>

#include "mri_core_def.h"
#include "mri_core_data.h"
#include "mri_core_utility.h"
#include "mri_core_acquisition_bucket.h"

#include "ImageIOAnalyze.h"

namespace Gadgetron {

    class EXPORTGADGETSMRICORE CoilSenMapGadget : public Gadget1<IsmrmrdAcquisitionBucket>
    {
    public:
        GADGET_DECLARE( CoilSenMapGadget );

        typedef Gadget1< IsmrmrdAcquisitionBucket > BaseClass;

        CoilSenMapGadget();
        virtual ~CoilSenMapGadget();

    protected:

        GADGET_PROPERTY(verbose, bool, "Whether to print more information", false);
        GADGET_PROPERTY(debug_folder, std::string, "If set, the debug output will be written out", "");
        GADGET_PROPERTY(perform_timing, bool, "Whether to perform timing on some computational steps", false);

        GADGET_PROPERTY( coil_sen_dependency_prefix, std::string, "Prefix of noise depencency file", "GadgetronCoilSenMap" );
        GADGET_PROPERTY( pass_nonconformant_data, bool, "Whether to pass data that does not conform", false );

        virtual int process_config( ACE_Message_Block* mb );
        virtual int process( GadgetContainerMessage<IsmrmrdAcquisitionBucket>* m1 );
        virtual int close(unsigned long flags);

        std::string generate_coil_sen_dependency_filename( const std::string& measurement_id );
        bool save_coil_sen_dependency();

        std::string coil_sen_dependency_folder_;
        std::string measurement_id_;
        std::string measurement_id_of_coil_sen_dependency_;
        std::string full_name_stored_coil_sen_dependency_;
        ISMRMRD::IsmrmrdHeader current_ismrmrd_header_;

        Gadgetron::ImageIOAnalyze gt_exporter_;
        std::string debug_folder_full_path_;

        // array for surface coil and body coil
        hoNDArray< std::complex<float> > scc_array_;
        std::vector<ISMRMRD::AcquisitionHeader> scc_header_;

        hoNDArray< std::complex<float> > body_array_;
        std::vector<ISMRMRD::AcquisitionHeader> body_header_;
    };
}
