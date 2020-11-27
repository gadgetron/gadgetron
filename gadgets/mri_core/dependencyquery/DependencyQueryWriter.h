/** \file   DependencyQueryWriter.h
    \brief  MRI image writer with meta attributes.
    \author Hui Xue
*/

#ifndef DependencyQueryWriter_H
#define DependencyQueryWriter_H

#include "GadgetMRIHeaders.h"
#include "ismrmrd/meta.h"
#include "gadgetron_mricore_export.h"

#include <ismrmrd/ismrmrd.h>
#include <complex>
#include "Dependency.h"

namespace Gadgetron
{

    class EXPORTGADGETSMRICORE DependencyQueryWriter : public Core::TypedWriter<DependencyQuery::Dependency>
    {
    protected:
        void serialize(std::ostream &stream, const DependencyQuery::Dependency&) override;

    public:
//        virtual int write(ACE_SOCK_Stream* sock, ACE_Message_Block* mb);
    };

}
#endif
