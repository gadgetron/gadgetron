#ifndef GRAPPAUNMIXINGGADGET_H_
#define GRAPPAUNMIXINGGADGET_H_

#include "gadgetron_grappa_export.h"
#include "Gadget.h"
#include "hoNDArray.h"
#include "ismrmrd/ismrmrd.h"
#include "GrappaWeights.h"

#include <complex>

namespace Gadgetron {

    struct EXPORTGADGETSGRAPPA GrappaUnmixingJob {
        boost::shared_ptr<GrappaWeights<float> > weights_;
    };

    class EXPORTGADGETSGRAPPA GrappaUnmixingGadget
            : public Gadget3<GrappaUnmixingJob, ISMRMRD::ImageHeader, hoNDArray<std::complex<float> > > {
    public:
        GADGET_DECLARE(GrappaUnmixingGadget);

        GrappaUnmixingGadget();

        virtual ~GrappaUnmixingGadget();

    protected:
        virtual int process(GadgetContainerMessage<GrappaUnmixingJob> *unmixing_job_message,
                            GadgetContainerMessage<ISMRMRD::ImageHeader> *image_header_message,
                            GadgetContainerMessage<hoNDArray<std::complex<float> > > *image_data_message);
    };
}

#endif /* GRAPPAUNMIXINGGADGET_H_ */
