/** \file   GenericReconBase.h
    \brief  This serves an optional base class gadget for the generic chain.
            Some common functionalities are implemented here and can be reused in specific recon gadgets.
            This gadget is instantiated for IsmrmrdReconData and IsmrmrdImageArray
    \author Hui Xue
*/

#pragma once

#include <boost/range/adaptor/strided.hpp>
#include <range/v3/action.hpp>
#include <range/v3/numeric.hpp>
#include <range/v3/view.hpp>

#include <complex>
#include "gadgetron_mricore_export.h"
#include "Gadget.h"
#include "GadgetronTimer.h"

#include "ismrmrd/ismrmrd.h"
#include "ismrmrd/xml.h"
#include "ismrmrd/meta.h"

#include "mri_core_def.h"
#include "mri_core_data.h"
#include "mri_core_utility.h"

#include "ImageIOAnalyze.h"

#include "gadgetron_sha1.h"

namespace Gadgetron {

    template <typename T> 
    class EXPORTGADGETSMRICORE GenericReconBase : public Gadget1<T>
    {
    public:
        GADGET_DECLARE(GenericReconBase);

        typedef Gadget1<T> BaseClass;

        GenericReconBase();
        ~GenericReconBase();

        /// ------------------------------------------------------------------------------------
        /// debug and timing
        GADGET_PROPERTY(verbose, bool, "Whether to print more information", false);
        GADGET_PROPERTY(debug_folder, std::string, "If set, the debug output will be written out", "");
        GADGET_PROPERTY(perform_timing, bool, "Whether to perform timing on some computational steps", false);

        /// ms for every time tick
        GADGET_PROPERTY(time_tick, float, "Time tick in ms", 2.5);

    protected:

        // number of encoding spaces in the protocol
        size_t num_encoding_spaces_;

        // number of times the process function is called
        size_t process_called_times_;

        // --------------------------------------------------
        // variables for debug and timing
        // --------------------------------------------------

        // clock for timing
        Gadgetron::GadgetronTimer gt_timer_local_;
        Gadgetron::GadgetronTimer gt_timer_;

        // debug folder
        std::string debug_folder_full_path_;

        // exporter
        Gadgetron::ImageIOAnalyze gt_exporter_;

        // --------------------------------------------------
        // gadget functions
        // --------------------------------------------------
        virtual int process_config(ACE_Message_Block* mb);
        virtual int process(GadgetContainerMessage<T>* m1);
    };

    class EXPORTGADGETSMRICORE GenericReconKSpaceReadoutBase :public GenericReconBase < ISMRMRD::AcquisitionHeader >
    {
    public:
        GADGET_DECLARE(GenericReconKSpaceReadoutBase);

        typedef GenericReconBase < ISMRMRD::AcquisitionHeader > BaseClass;

        GenericReconKSpaceReadoutBase();
        virtual ~GenericReconKSpaceReadoutBase();
    };

    class EXPORTGADGETSMRICORE GenericReconDataBase :public GenericReconBase < IsmrmrdReconData >
    {
    public:
        GADGET_DECLARE(GenericReconDataBase);

        typedef GenericReconBase < IsmrmrdReconData > BaseClass;

        GenericReconDataBase();
        virtual ~GenericReconDataBase();
    };

    class EXPORTGADGETSMRICORE GenericReconImageBase :public GenericReconBase < IsmrmrdImageArray >
    {
    public:
        GADGET_DECLARE(GenericReconImageBase);

        typedef GenericReconBase < IsmrmrdImageArray > BaseClass;

        GenericReconImageBase();
        virtual ~GenericReconImageBase();
    };

    class EXPORTGADGETSMRICORE GenericReconImageHeaderBase :public GenericReconBase < ISMRMRD::ImageHeader >
    {
    public:
        GADGET_DECLARE(GenericReconImageHeaderBase);

        typedef GenericReconBase < ISMRMRD::ImageHeader > BaseClass;

        GenericReconImageHeaderBase();
        virtual ~GenericReconImageHeaderBase();
    };
}
