/** \file   GenericReconBase.h
    \brief  This serves an optional base class gadget for the generic chain.
            Some common functionalities are implemented here and can be reused in specific recon gadgets.
            This gadget is instantiated for ReconData and ImageArray
    \author Hui Xue
*/

#pragma once

#include <complex>
#include "Gadget.h"
#include "GadgetronTimer.h"

#include "mri_core_def.h"
#include "mri_core_utility.h"
#include "mri_core_stream.h"

#include "ImageIOAnalyze.h"

#include "pingvin_sha1.h"

#include "GenericReconStreamDef.h"

namespace Gadgetron {

    template <typename T>
    class GenericReconBase : public Gadget1<T>
    {
    public:
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
        // data stream
        // --------------------------------------------------
        GenericReconMrdStreamer gt_streamer_;

        // --------------------------------------------------
        // gadget functions
        // --------------------------------------------------
        virtual int process_config(const mrd::Header& header);
        virtual int process(GadgetContainerMessage<T>* m1);
        virtual int close(unsigned long flags);
    };

    class GenericReconKSpaceReadoutBase :public GenericReconBase < mrd::AcquisitionHeader >
    {
    public:
        typedef GenericReconBase < mrd::AcquisitionHeader > BaseClass;

        GenericReconKSpaceReadoutBase();
        virtual ~GenericReconKSpaceReadoutBase();
        virtual int close(unsigned long flags) { return BaseClass::close(flags); }
    };

    class GenericReconDataBase :public GenericReconBase < mrd::ReconData >
    {
    public:
        typedef GenericReconBase < mrd::ReconData > BaseClass;

        GenericReconDataBase();
        virtual ~GenericReconDataBase();
        virtual int close(unsigned long flags) { return BaseClass::close(flags); }
    };

    class GenericReconImageBase :public GenericReconBase < mrd::ImageArray >
    {
    public:
        typedef GenericReconBase < mrd::ImageArray > BaseClass;

        GenericReconImageBase();
        virtual ~GenericReconImageBase();
        virtual int close(unsigned long flags) { return BaseClass::close(flags); }
    };

    class GenericReconImageHeaderBase :public GenericReconBase < mrd::ImageHeader >
    {
    public:
        typedef GenericReconBase < mrd::ImageHeader > BaseClass;

        GenericReconImageHeaderBase();
        virtual ~GenericReconImageHeaderBase();
        virtual int close(unsigned long flags) { return BaseClass::close(flags); }
    };
}
