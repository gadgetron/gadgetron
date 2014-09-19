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
