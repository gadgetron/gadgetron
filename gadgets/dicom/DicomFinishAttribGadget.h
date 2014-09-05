/** \file       DicomFinishAttribGadget.h
    \brief      Assemble the dicom images and send out 

                The dicom image is sent out with message id -> dicom image -> dicom image name -> meta attributes
    \author     Hui Xue
*/

#ifndef DICOMFINISHATTRIBGADGET_H
#define DICOMFINISHATTRIBGADGET_H

#include "gadgetron_dicom_export.h"

#include "Gadget.h"
#include "hoNDArray.h"
#include "ismrmrd_meta.h"
#include "GadgetMRIHeaders.h"
#include "ismrmrd.h"
#include "GadgetStreamController.h"

#include "dcmtk/config/osconfig.h"
#include "dcmtk/ofstd/ofstdinc.h"
#define INCLUDE_CSTDLIB
#define INCLUDE_CSTDIO
#define INCLUDE_CSTRING
#include "dcmtk/dcmdata/dctk.h"
#include "dcmtk/dcmdata/dcostrmb.h"

#include <string>
#include <map>
#include <complex>

namespace Gadgetron
{

template <typename T>
class EXPORTGADGETSDICOM DicomFinishAttribGadget : public Gadget3<ISMRMRD::ImageHeader, hoNDArray< T >, ISMRMRD::MetaContainer >
{
    public:

        typedef Gadget3<ISMRMRD::ImageHeader, hoNDArray< T >, ISMRMRD::MetaContainer > BaseClass;

        DicomFinishAttribGadget<T>()
            : BaseClass()
            , dcmFile()
            , seriesIUIDRoot()
        { }

    protected:

        virtual int process_config(ACE_Message_Block * mb);
        virtual int process(GadgetContainerMessage<ISMRMRD::ImageHeader>* m1, GadgetContainerMessage< hoNDArray< T > >* m2, GadgetContainerMessage<ISMRMRD::MetaContainer>* m3);

    private:
        DcmFileFormat dcmFile;
        std::string seriesIUIDRoot;
        long initialSeriesNumber;
        std::map <unsigned int, std::string> seriesIUIDs;
};

class EXPORTGADGETSDICOM DicomFinishAttribGadgetUSHORT :
    public DicomFinishAttribGadget<ACE_UINT16>
{
    public:
        GADGET_DECLARE(DicomFinishAttribGadgetUSHORT);
};

class EXPORTGADGETSDICOM DicomFinishAttribGadgetFLOAT :
    public DicomFinishAttribGadget<float>
{
    public:
        GADGET_DECLARE(DicomFinishAttribGadgetFLOAT);
};

} /* namespace Gadgetron */

#endif // DICOMFINISHATTRIBGADGET_H
