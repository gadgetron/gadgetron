#ifndef PhysioInterpolationGadget_H
#define PhysioInterpolationGadget_H

#include "Gadget.h"
#include "GadgetMRIHeaders.h"
#include "hoNDArray.h"
#include "gadgetron_mricore_export.h"

#include <ismrmrd/ismrmrd.h>
#include <complex>

namespace Gadgetron{  

    class EXPORTGADGETSMRICORE PhysioInterpolationGadget : public Gadget2<ISMRMRD::ImageHeader, hoNDArray< std::complex<float> > >
    {
    public:
        GADGET_DECLARE(PhysioInterpolationGadget);

        PhysioInterpolationGadget();
        virtual ~PhysioInterpolationGadget();

        inline unsigned short get_number_of_phases() { return phases_to_reconstruct_; }

    protected:
        GADGET_PROPERTY(physiology_time_index, int, "Physiology time index", 0);
        GADGET_PROPERTY_LIMITS(mode, int, "Mode, 0=seperate series for each RR, 1=First complete RR only", 0, GadgetPropertyLimitsEnumeration, 0, 1);
        GADGET_PROPERTY(phases, unsigned short, "Number of cardiac phases", 30);
        GADGET_PROPERTY(first_beat_on_trigger, bool, "Indicates that acquisition was started on trigger", false);
        GADGET_PROPERTY_LIMITS(interp_method, std::string, "Interpolation method", "Spline", GadgetPropertyLimitsEnumeration, "Spline", "BSpline", "");
        GADGET_PROPERTY(time_stamp_resolution_, double, "Time stamp resolution in ms", 2.5);

        virtual int process_config(ACE_Message_Block* mb);

        virtual int process(GadgetContainerMessage< ISMRMRD::ImageHeader >* m1, GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2);

        virtual int close(unsigned long flags); //All the work is done here in this Gadget

        unsigned short phys_time_index_;
        unsigned short phases_to_reconstruct_;
        unsigned short mode_; //0=seperate series for each complete RR,
                              //1=First complete RR interval only

        // true, if the first beat is on trigger
        /// false, the first beat will be ignored
        bool first_beat_on_trigger_;

        // interpolation method, "Spline" or "BSpline"
        std::string interp_method_;

    private:

        std::vector< boost::shared_ptr< ACE_Message_Queue<ACE_MT_SYNCH> > > buffer_;
        std::vector< std::vector<float> > time_stamps_;

        size_t slc_limit_;

        bool image_with_attrib_;
    };
}

#endif //PhysioInterpolationGadget_H
