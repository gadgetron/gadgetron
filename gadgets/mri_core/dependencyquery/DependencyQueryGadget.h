#pragma once

#include "Gadget.h"
#include "hoNDArray.h"
#include "gadgetron_mricore_export.h"

#include <ismrmrd/ismrmrd.h>
#include <complex>
#include <ctime>
#include "GadgetMRIHeaders.h"
#include "ismrmrd/meta.h"

namespace Gadgetron
{
    class EXPORTGADGETSMRICORE DependencyQueryGadget : public Gadget2<ISMRMRD::AcquisitionHeader,hoNDArray< std::complex<float> > >
    {
    public:
        GADGET_DECLARE(DependencyQueryGadget);

        typedef std::complex<float> ValueType;

        typedef Gadget2<ISMRMRD::AcquisitionHeader,hoNDArray< ValueType > > BaseClass;

        DependencyQueryGadget();
        virtual ~DependencyQueryGadget();

        virtual int close(unsigned long flags);

    protected:
	GADGET_PROPERTY(noise_dependency_prefix, std::string, "Prefix on noise dependency file", "");
	GADGET_PROPERTY(noise_dependency_attrib_name, std::string, "Noise dependency attribute name", "");
	GADGET_PROPERTY(clean_storage_while_query, bool, "Clean storage while querying", false);
	GADGET_PROPERTY(time_limit_in_storage, float, "Time limit for storing noise dependency", 0);

        // if true, the old stored file will be deleted while querying
        bool clean_storage_while_query_;

        // in the unit of hours, how long a file is allowed to be in the storage
        double time_limit_in_storage_;

        // current time, year/month/day/hour/min/second
        std::time_t curr_time_UTC_;

        bool processed_in_close_;

        std::string noise_dependency_folder_;
        std::string noise_dependency_prefix_;

        std::string noise_dependency_attrib_name_;

        virtual int process_config(ACE_Message_Block* mb);

        virtual int process(GadgetContainerMessage<ISMRMRD::AcquisitionHeader>* m1,
            GadgetContainerMessage< hoNDArray< ValueType > >* m2);
    };
}
