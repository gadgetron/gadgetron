/** \file   DicomFinishGadget.h
\brief      Assemble the dicom images and send out

The dicom image is sent out with message id -> dicom image -> dicom image name -> meta attributes
\author     Hui Xue
*/

#ifndef DICOMFINISHGADGET_H
#define DICOMFINISHGADGET_H

#include "gadgetron_dicom_export.h"

#include "Gadget.h"
#include "hoNDArray.h"
#include "ismrmrd/meta.h"
#include "GadgetMRIHeaders.h"
#include "ismrmrd/ismrmrd.h"

#include "dcmtk/config/osconfig.h"
#include "dcmtk/ofstd/ofstdinc.h"
#include "dcmtk/dcmdata/dctk.h"
#include "dcmtk/dcmdata/dcostrmb.h"

#include "mri_core_def.h"

#include "dicom_ismrmrd_utility.h"

#include <string>
#include <map>
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <cstring>

namespace Gadgetron
{
    class EXPORTGADGETSDICOM DicomFinishGadget : public Gadget1< ISMRMRD::ImageHeader >
    {
    public:

        typedef Gadget1<ISMRMRD::ImageHeader> BaseClass;

        GADGET_DECLARE(DicomFinishGadget);

        DicomFinishGadget()
            : BaseClass()
            , dcmFile()
            , initialSeriesNumber(0)
            , seriesIUIDRoot()
        { }

    protected:

        virtual int process_config(ACE_Message_Block * mb);
        virtual int process(GadgetContainerMessage<ISMRMRD::ImageHeader>* m1);

        template <typename T>
        int write_data_attrib(GadgetContainerMessage<ISMRMRD::ImageHeader>* m1, GadgetContainerMessage< hoNDArray< T > >* m2)
        {

            GadgetContainerMessage< ISMRMRD::MetaContainer >* m3 = AsContainerMessage< ISMRMRD::MetaContainer >(m2->cont());

            std::string filename;

            if (m3)
            {
                ISMRMRD::MetaContainer* img_attrib = m3->getObjectPtr();

                size_t n;

                size_t num = img_attrib->length(GADGETRON_DATA_ROLE);

                std::vector<std::string> dataRole;
                if (num == 0)
                {
                    dataRole.push_back("Image");
                }
                else
                {
                    dataRole.resize(num);
                    for (n = 0; n < num; n++)
                    {
                        dataRole[n] = std::string(img_attrib->as_str(GADGETRON_DATA_ROLE, n));
                    }
                }

                long imageNumber = img_attrib->as_long(GADGETRON_IMAGENUMBER, 0);

                long slc, con, phs, rep, set, ave;
                slc = m1->getObjectPtr()->slice;
                con = m1->getObjectPtr()->contrast;
                phs = m1->getObjectPtr()->phase;
                rep = m1->getObjectPtr()->repetition;
                set = m1->getObjectPtr()->set;
                ave = m1->getObjectPtr()->average;

                std::ostringstream ostr;

                for (n = 0; n < dataRole.size(); n++)
                {
                    ostr << dataRole[n] << "_";
                }

                ostr << "SLC" << slc << "_"
                    << "CON" << con << "_"
                    << "PHS" << phs << "_"
                    << "REP" << rep << "_"
                    << "SET" << set << "_"
                    << "AVE" << ave << "_"
                    << imageNumber;

                filename = ostr.str();
            }
            else
            {
                std::ostringstream ostr;
                ostr << "Image_" << m1->getObjectPtr()->image_index << "_" << m1->getObjectPtr()->image_series_index;
                filename = ostr.str();
            }

            GadgetContainerMessage<std::string>* mfilename = new GadgetContainerMessage<std::string>();
            *(mfilename->getObjectPtr()) = filename;

            // --------------------------------------------------

            unsigned short series_number = m1->getObjectPtr()->image_series_index + 1;

            // Try to find an already-generated Series Instance UID in our map
            std::map<unsigned int, std::string>::iterator it = seriesIUIDs.find(series_number);

            if (it == seriesIUIDs.end()) {
                // Didn't find a Series Instance UID for this series number
                char prefix[32];
                char newuid[96];
                if (seriesIUIDRoot.length() > 20) {
                    memcpy(prefix, seriesIUIDRoot.c_str(), 20);
                    prefix[20] = '\0';
                    dcmGenerateUniqueIdentifier(newuid, prefix);
                }
                else {
                    dcmGenerateUniqueIdentifier(newuid);
                }
                seriesIUIDs[series_number] = std::string(newuid);
            }

            // --------------------------------------------------

            if(m3)
            {
                Gadgetron::write_ismrmd_image_into_dicom(*m1->getObjectPtr(), *m2->getObjectPtr(), xml, *m3->getObjectPtr(), seriesIUIDs[series_number], dcmFile);
            }
            else
            {
                ISMRMRD::MetaContainer attrib;
                Gadgetron::write_ismrmd_image_into_dicom(*m1->getObjectPtr(), *m2->getObjectPtr(), xml, attrib, seriesIUIDs[series_number], dcmFile);
            }

            // --------------------------------------------------
            /* release the old data array */
            m2->cont(NULL); // still need m3
            m1->release();

            GadgetContainerMessage<DcmFileFormat>* mdcm = new GadgetContainerMessage<DcmFileFormat>();

            *mdcm->getObjectPtr() = dcmFile;
            mdcm->cont(mfilename);

            if (m3)
            {
                mfilename->cont(m3);
            }

            this->next()->putq(mdcm);

            return GADGET_OK;
        }

    private:
        ISMRMRD::IsmrmrdHeader xml;
        DcmFileFormat dcmFile;
        std::string seriesIUIDRoot;
        long initialSeriesNumber;
        std::map <unsigned int, std::string> seriesIUIDs;
    };

} /* namespace Gadgetron */

#endif // DICOMFINISHGADGET_H
