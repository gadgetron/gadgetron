/** \file   DependencyQueryWriter.h
    \brief  MRI image writer with meta attributes.
    \author Hui Xue
*/

#ifndef DependencyQueryWriter_H
#define DependencyQueryWriter_H

#include "GadgetMessageInterface.h"
#include "GadgetMRIHeaders.h"
#include "hoNDMetaAttributes.h"
#include "gadgetron_mricore_export.h"

#include <ismrmrd.h>
#include <complex>

namespace Gadgetron
{

    class EXPORTGADGETSMRICORE DependencyQueryWriter : public GadgetMessageWriter
    {
    public:
        GADGETRON_WRITER_DECLARE(DependencyQueryWriter)
        virtual int write(ACE_SOCK_Stream* sock, ACE_Message_Block* mb);
    };

}
#endif
