/** \file   MRIImageAttribWriter.h
    \brief  MRI image writer with meta attributes.
    \author Hui Xue
*/

#ifndef MRIImageAttribWriter_H
#define MRIImageAttribWriter_H

#include "GadgetMessageInterface.h"
#include "GadgetMRIHeaders.h"
#include "ismrmrd/meta.h"
#include "gadgetron_mricore_export.h"

#include <ismrmrd/ismrmrd.h>
#include <complex>

namespace Gadgetron{

    template<typename T> class MRIImageAttribWriter : public GadgetMessageWriter
    {
    public:
        virtual int write(ACE_SOCK_Stream* sock, ACE_Message_Block* mb);
    };

    class EXPORTGADGETSMRICORE MRIImageAttribWriterUSHORT : public MRIImageAttribWriter<ACE_UINT16>
    {
    public:
        GADGETRON_WRITER_DECLARE(MRIImageAttribWriterUSHORT);
    };

    class EXPORTGADGETSMRICORE MRIImageAttribWriterSHORT : public MRIImageAttribWriter<ACE_INT16>
    {
    public:
        GADGETRON_WRITER_DECLARE(MRIImageAttribWriterSHORT);
    };

    class EXPORTGADGETSMRICORE MRIImageAttribWriterFLOAT : public MRIImageAttribWriter<float>
    {
    public:
        GADGETRON_WRITER_DECLARE(MRIImageAttribWriterFLOAT);
    };

    class EXPORTGADGETSMRICORE MRIImageAttribWriterCPLX : public MRIImageAttribWriter< std::complex<float> >
    {
    public:
        GADGETRON_WRITER_DECLARE(MRIImageAttribWriterCPLX);
    };
}
#endif
