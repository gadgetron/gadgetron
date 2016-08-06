/** \file   DicomIsmrmrdImageDumpGadget.h
    \brief  Dump the incoming ismrmrd images as dicom images
    \author Hui Xue
*/

#ifndef DicomIsmrmrdImageDumpGadget_H
#define DicomIsmrmrdImageDumpGadget_H

#include "gadgetron_dicom_export.h"

#include "Gadget.h"
#include "hoNDArray.h"
#include "ismrmrd/meta.h"
#include "GadgetMRIHeaders.h"
#include "ismrmrd/ismrmrd.h"
#include "GadgetStreamController.h"

#include "dcmtk/config/osconfig.h"
#include "dcmtk/ofstd/ofstdinc.h"
#define INCLUDE_CSTDLIB
#define INCLUDE_CSTDIO
#define INCLUDE_CSTRING
#include "dcmtk/dcmdata/dctk.h"
#include "dcmtk/dcmdata/dcostrmb.h"

#include "mri_core_def.h"

#include "dicom_ismrmrd_utility.h"

#include <string>
#include <map>
#include <complex>

namespace Gadgetron
{
    class EXPORTGADGETSDICOM DicomIsmrmrdImageDumpGadget : public Gadget1< ISMRMRD::ImageHeader >
    {
    public:

        typedef Gadget1<ISMRMRD::ImageHeader> BaseClass;

        GADGET_DECLARE(DicomIsmrmrdImageDumpGadget);

        DicomIsmrmrdImageDumpGadget();
        virtual ~DicomIsmrmrdImageDumpGadget();

    protected:

        GADGET_PROPERTY(folder,           std::string, "Base folder for dump file", ".");
        GADGET_PROPERTY(file_prefix,      std::string, "Prefix for dump file", "ISMRMRD_DICOM_DUMP");
        GADGET_PROPERTY(append_id,        bool,        "ISMRMRD measurement ID to file name prefix (if available)", true);
        GADGET_PROPERTY(append_timestamp, bool,        "Append timestamp to file name prefix", true);

        virtual int process_config(ACE_Message_Block * mb);
        virtual int process(GadgetContainerMessage<ISMRMRD::ImageHeader>* m1);

        template <typename T>
        int dump_image(GadgetContainerMessage<ISMRMRD::ImageHeader>* m1, GadgetContainerMessage< hoNDArray< T > >* m2)
        {
            try{
                GadgetContainerMessage< ISMRMRD::MetaContainer >* m3 = AsContainerMessage< ISMRMRD::MetaContainer >(m2->cont());

                std::string filename;

                ISMRMRD::MetaContainer* img_attrib = NULL;
                if (m3)
                {
                    img_attrib = m3->getObjectPtr();

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

                // --------------------------------------------------

                unsigned short series_number = m1->getObjectPtr()->image_series_index + 1;

                // Try to find an already-generated Series Instance UID in our map
                std::map<unsigned int, std::string>::iterator it = seriesIUIDs.find(series_number);

                if (it == seriesIUIDs.end()) {
                    // Didn't find a Series Instance UID for this series number
                    char newuid[96];
                    dcmGenerateUniqueIdentifier(newuid);
                    seriesIUIDs[series_number] = std::string(newuid);
                }

                // --------------------------------------------------

                ISMRMRD::ImageHeader mh = *m1->getObjectPtr();

                if(m3)
                {
                    Gadgetron::write_ismrmd_image_into_dicom(mh, *m2->getObjectPtr(), xml, *m3->getObjectPtr(), seriesIUIDs[series_number], dcmFile);
                }
                else
                {
                    ISMRMRD::MetaContainer attrib;
                    Gadgetron::write_ismrmd_image_into_dicom(mh, *m2->getObjectPtr(), xml, attrib, seriesIUIDs[series_number], dcmFile);
                }

                DcmTagKey key;
                DcmDataset *dataset = dcmFile.getDataset();

                if(!device_.empty())
                {
                    key.set(0x0018, 0x1000);
                    write_dcm_string(dataset, key, device_.c_str());
                }

                if(!patient_.empty())
                {
                    // patient id
                    key.set(0x0010, 0x0020);
                    write_dcm_string(dataset, key, patient_.c_str());

                    //// patient name
                    //key.set(0x0010, 0x0010);
                    //write_dcm_string(dataset, key, patient_.c_str());
                }

                if(!study_.empty())
                {
                    key.set(0x0020, 0x000D);
                    write_dcm_string(dataset, key, study_.c_str());

                    std::string studyID = study_ + "19000101";
                    key.set(0x0020, 0x0010);
                    write_dcm_string(dataset, key, studyID.c_str());
                }

                if(!patient_string_.empty())
                {
                    // patient name
                    key.set(0x0010, 0x0010);
                    write_dcm_string(dataset, key, patient_string_.c_str());
                }

                if(!protocol_name_.empty())
                {
                    if(img_attrib!=NULL)
                    {
                        size_t N = img_attrib->length(GADGETRON_SEQUENCEDESCRIPTION);
                        if(N>0)
                        {
                            std::string str(std::string(img_attrib->as_str(GADGETRON_SEQUENCEDESCRIPTION, 0)));
                            for (size_t n=1; n<N; n++)
                            {
                                str.append("_");
                                str.append( std::string(img_attrib->as_str(GADGETRON_SEQUENCEDESCRIPTION, n)) );
                            }

                            std::string str_seq_description = protocol_name_ + str;

                            key.set(0x0008, 0x103E);
                            write_dcm_string(dataset, key, str_seq_description.c_str());
                        }
                    }
                }

                if(!xml.studyInformation.is_present() || xml.studyInformation.get().studyDate.is_present())
                {
                    // set current date as study date
                    key.set(0x0008, 0x0020);
                    std::string d = this->get_date_string();
                    write_dcm_string(dataset, key, d.c_str());
                }

                // --------------------------------------------------
                // dump dicom image

                DcmFileFormat mdcm = dcmFile;

                mdcm.transferInit();

                long buffer_length = mdcm.calcElementLength(EXS_LittleEndianExplicit, EET_ExplicitLength) * 2;
                std::vector<char> bufferChar(buffer_length);
                char* buffer = &bufferChar[0];

                DcmOutputBufferStream out_stream(buffer, buffer_length);

                OFCondition status;

                status = mdcm.write(out_stream, EXS_LittleEndianExplicit, EET_ExplicitLength, NULL);
                if (!status.good())
                {
                    GERROR_STREAM("Dump dicom image failed : " << filename);
                    GADGET_THROW("Failed to write DcmFileFormat to DcmOutputStream ... ");
                }

                void *serialized = NULL;
                offile_off_t serialized_length = 0;
                out_stream.flushBuffer(serialized, serialized_length);

                mdcm.transferEnd();

                std::string file_full_name = dicom_folder_ + "/" + filename + ".dcm";

                std::ofstream outfile;
                outfile.open(file_full_name.c_str(), std::ios::out | std::ios::binary);

                outfile.write( (char*)serialized, serialized_length);
                outfile.close();
            }
            catch(...)
            {
                GERROR_STREAM("Exceptions happened in DicomIsmrmrdImageDumpGadget::dump_image(...) ... ");
                return GADGET_FAIL;
            }

            return GADGET_OK;
        }

    private:
        ISMRMRD::IsmrmrdHeader xml;
        DcmFileFormat dcmFile;
        std::map <unsigned int, std::string> seriesIUIDs;
        std::string dicom_folder_;
        std::string device_;
        std::string patient_;
        std::string patient_string_;
        std::string patient_position_;
        std::string study_;
        std::string measurement_;
        std::string protocol_name_;

        std::string get_date_string();
        std::string get_time_string();
    };

} /* namespace Gadgetron */

#endif // DicomIsmrmrdImageDumpGadget_H
