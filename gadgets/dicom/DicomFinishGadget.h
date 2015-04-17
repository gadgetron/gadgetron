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
#include "GadgetStreamController.h"

#include "dcmtk/config/osconfig.h"
#include "dcmtk/ofstd/ofstdinc.h"
#define INCLUDE_CSTDLIB
#define INCLUDE_CSTDIO
#define INCLUDE_CSTRING
#include "dcmtk/dcmdata/dctk.h"
#include "dcmtk/dcmdata/dcostrmb.h"

#include "mri_core_def.h"

#include <string>
#include <map>
#include <complex>

namespace Gadgetron
{

    // Used for windowing using short ints
#define PIX_RANGE_MAX    (+32767)
#define PIX_RANGE_MIN    (-32768)


    // Writes a DICOM string value at the given location in the header
    // Saves keystrokes
#define WRITE_DCM_STRING(k, s)    \
    do {                                                                    \
        status = dataset->putAndInsertString(k, s);            \
        if (!status.good()) {                                               \
            GDEBUG("Failed to insert DICOM field (0x%04X,0x%04X) at "\
                "line %u\n", k.getGroup(), k.getElement(), __LINE__);       \
            return GADGET_FAIL;                                             \
                        }                                                                   \
            } while (0)

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

            GadgetContainerMessage<hoNDArray< ACE_INT16 > > *pixels = new GadgetContainerMessage<hoNDArray< ACE_INT16 > >();
            boost::shared_ptr< std::vector<size_t> > dims = m2->getObjectPtr()->get_dimensions();

            try {
                pixels->getObjectPtr()->create(dims.get());
            }
            catch (bad_alloc&) {
                GDEBUG("Unable to create short storage in DicomFinishGadget");
                return GADGET_FAIL;
            }

            /* create ImageHeader and hoNDArray pointers for better readability */
            hoNDArray<ACE_INT16>* data = pixels->getObjectPtr();

            /* grab pointers to both the original and new data arrays
            * The original is of type T
            * The new is of type ACE_INT16 */
            T *src = m2->getObjectPtr()->get_data_ptr();
            ACE_INT16 *dst = data->get_data_ptr();

            /* Convert/cast each element in the data array
            * and simultaneously find the min/max pixel value, which
            * will be used later for some crude windowing */
            T min_pix_val, max_pix_val, sum_pix_val = 0;
            if (pixels->getObjectPtr()->get_number_of_elements() > 0)
            {
                min_pix_val = src[0];
                max_pix_val = src[0];
            }

            for (unsigned long i = 0; i < pixels->getObjectPtr()->get_number_of_elements(); i++)
            {
                T pix_val = src[i];
                // search for minimum and maximum pixel values
                if (pix_val < min_pix_val) min_pix_val = pix_val;
                if (pix_val > max_pix_val) max_pix_val = pix_val;
                sum_pix_val += pix_val / 4; // scale by 25% to avoid overflow

                // copy/cast the pixel value to a short int
                dst[i] = static_cast<ACE_INT16>(pix_val);
            }
            T mean_pix_val = (T)((sum_pix_val * 4) / (T)pixels->getObjectPtr()->get_number_of_elements());

            /* replace the old 'message2' with the new data */
            m1->cont(pixels);
            /* release the old data array */
            m2->cont(NULL);
            m2->release();
            /* update the image data_type.
            * There is currently no SIGNED SHORT type so this will have to suffice */
            m1->getObjectPtr()->data_type = ISMRMRD::ISMRMRD_USHORT;

            unsigned int BUFSIZE = 1024;
            char *buf = new char[BUFSIZE];
            OFCondition status;
            DcmTagKey key;
            DcmDataset *dataset = dcmFile.getDataset();

            // Echo Number
            // TODO: it is often the case the img->contrast is not properly set
            // likely due to the allocated ISMRMRD::ImageHeader being uninitialized
            key.set(0x0018, 0x0086);
            ACE_OS::snprintf(buf, BUFSIZE, "%d", m1->getObjectPtr()->contrast);
            WRITE_DCM_STRING(key, buf);

            // Acquisition Matrix ... Image Dimensions
            // Defined as: [frequency rows, frequency columns, phase rows, phase columns]
            // But at this point in the gadget I don't know the frequency encode direction
            key.set(0x0018, 0x1310);
            ACE_UINT16 im_dim[4] = { 0, 0, 0, 0 };
            /* if (frequence_encode_dir == "ROW")) {
            // frequency encoding direction is ROW
            im_dim[1] = img->matrix_size[0];
            im_dim[2] = img->matrix_size[1];
            } */
            // frequency encoding direction is COLUMN
            /*im_dim[0] = img->matrix_size[0];
            im_dim[3] = img->matrix_size[1];*/

            im_dim[1] = m1->getObjectPtr()->matrix_size[0];
            im_dim[2] = m1->getObjectPtr()->matrix_size[1];

            status = dataset->putAndInsertUint16Array(key, im_dim, 4);
            if (!status.good())
            {
                GDEBUG("Failed to stuff image dimensions\n");
                return GADGET_FAIL;
            }

            // Series Number
            // Only write a number if the image_series_index is positive and non-zero
            key.set(0x0020, 0x0011);
            ACE_OS::snprintf(buf, BUFSIZE, "%ld", this->initialSeriesNumber * 100 + m1->getObjectPtr()->image_series_index);
            WRITE_DCM_STRING(key, buf);

            // Image Number
            key.set(0x0020, 0x0013);
            ACE_OS::snprintf(buf, BUFSIZE, "%d", m1->getObjectPtr()->image_index + 1);
            WRITE_DCM_STRING(key, buf);

            // Image Position (Patient)
            float corner[3];

            corner[0] = m1->getObjectPtr()->position[0] -
                (m1->getObjectPtr()->field_of_view[0] / 2.0f) * m1->getObjectPtr()->read_dir[0] -
                (m1->getObjectPtr()->field_of_view[1] / 2.0f) * m1->getObjectPtr()->phase_dir[0];
            corner[1] = m1->getObjectPtr()->position[1] -
                (m1->getObjectPtr()->field_of_view[0] / 2.0f) * m1->getObjectPtr()->read_dir[1] -
                (m1->getObjectPtr()->field_of_view[1] / 2.0f) * m1->getObjectPtr()->phase_dir[1];
            corner[2] = m1->getObjectPtr()->position[2] -
                (m1->getObjectPtr()->field_of_view[0] / 2.0f) * m1->getObjectPtr()->read_dir[2] -
                (m1->getObjectPtr()->field_of_view[1] / 2.0f) * m1->getObjectPtr()->phase_dir[2];

            key.set(0x0020, 0x0032);
            ACE_OS::snprintf(buf, BUFSIZE, "%.4f\\%.4f\\%.4f", corner[0], corner[1], corner[2]);
            WRITE_DCM_STRING(key, buf);

            // Image Orientation
            // read_dir, phase_dir, and slice_dir were calculated in
            // a DICOM/patient coordinate system, so just plug them in
            key.set(0x0020, 0x0037);
            ACE_OS::snprintf(buf, BUFSIZE, "%.4f\\%.4f\\%.4f\\%.4f\\%.4f\\%.4f",
                m1->getObjectPtr()->read_dir[0], m1->getObjectPtr()->read_dir[1], m1->getObjectPtr()->read_dir[2],
                m1->getObjectPtr()->phase_dir[0], m1->getObjectPtr()->phase_dir[1], m1->getObjectPtr()->phase_dir[2]);
            WRITE_DCM_STRING(key, buf);

            // Slice Location
            key.set(0x0020, 0x1041);
            ACE_OS::snprintf(buf, BUFSIZE, "%f", m1->getObjectPtr()->position[2]);
            WRITE_DCM_STRING(key, buf);

            // Columns
            key.set(0x0028, 0x0010);
            ACE_OS::snprintf(buf, BUFSIZE, "%d", m1->getObjectPtr()->matrix_size[1]);
            WRITE_DCM_STRING(key, buf);

            // Rows
            key.set(0x0028, 0x0011);
            ACE_OS::snprintf(buf, BUFSIZE, "%d", m1->getObjectPtr()->matrix_size[0]);
            WRITE_DCM_STRING(key, buf);

            //Number of frames
            if (m1->getObjectPtr()->matrix_size[2] > 1){ //Only write if we have more than 1 frame
            	key.set(0x0028,0x0008);
            	ACE_OS::snprintf(buf,BUFSIZE,"%d",m1->getObjectPtr()->matrix_size[2]);
            	WRITE_DCM_STRING(key,buf);
            }

            // Simple windowing using pixel values calculated earlier...
            int mid_pix_val = (int)(max_pix_val + min_pix_val) / 2;
            int window_center = (int)(mid_pix_val + mean_pix_val) / 2;
            int window_width_left = (int)(window_center - min_pix_val);
            int window_width_right = (int)(max_pix_val - window_center);
            int window_width = (window_width_right > window_width_left) ?
            window_width_right : window_width_left;

            // Window Center
            key.set(0x0028, 0x1050);
            ACE_OS::snprintf(buf, BUFSIZE, "%d", window_center);
            WRITE_DCM_STRING(key, buf);

            // Window Width
            key.set(0x0028, 0x1051);
            ACE_OS::snprintf(buf, BUFSIZE, "%d", window_width);
            WRITE_DCM_STRING(key, buf);

            // ACR_NEMA_2C_VariablePixelDataGroupLength
            key.set(0x7fe0, 0x0000);
            status = dataset->insertEmptyElement(key);
            if (!status.good()) {
                GDEBUG("Failed to write 0x7fe0 Group Length\n");
                return GADGET_FAIL;
            }

            // Pixel Data
            if ((unsigned long)m1->getObjectPtr()->matrix_size[0] * (unsigned long)m1->getObjectPtr()->matrix_size[1]*(unsigned long)m1->getObjectPtr()->matrix_size[2] !=
                data->get_number_of_elements()) {
                GDEBUG("Mismatch in image dimensions and available data\n");
                return GADGET_FAIL;
            }
            key.set(0x7fe0, 0x0010);
            status = dataset->putAndInsertUint16Array(key, (unsigned short *)data->get_data_ptr(), (unsigned long)data->get_number_of_elements());
            if (!status.good()) {
                GDEBUG("Failed to stuff Pixel Data\n");
                return GADGET_FAIL;
            }

            // Series Instance UID = generated here
            key.set(0x0020, 0x000E);
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
            WRITE_DCM_STRING(key, seriesIUIDs[series_number].c_str());

            // At a minimum, to put the DICOM image back into the database,
            // you must change the SOPInstanceUID.
            key.set(0x0008, 0x0018);        // SOPInstanceUID
            const char *root;
            if (seriesIUIDRoot.length() > 0) {
                root = std::string(seriesIUIDRoot, 0, 20).c_str();
            }
            else {
                root = "1.2.840.113619.2.156";
            }
            char newuid[65];
            dcmGenerateUniqueIdentifier(newuid, root);
            WRITE_DCM_STRING(key, newuid);

            //// set the private fields to store meta attributes
            //key.set(0x0051, 0x0000);
            //status = dataset->insertEmptyElement(key);
            //if (!status.good())
            //{
            //    GDEBUG("Failed to write 0x0051 Group Length\n");
            //    return GADGET_FAIL;
            //}

            //key.set(0x0051, 0x0019);
            //WRITE_DCM_STRING(key, buf+sizeof(size_t_type));

            //delete [] meta_buf;

            /* clean up the char[] we created for ACE_OS::snprintf */
            delete[] buf;

            GadgetContainerMessage<DcmFileFormat>* mdcm = new GadgetContainerMessage<DcmFileFormat>();

            *mdcm->getObjectPtr() = dcmFile;

            GadgetContainerMessage<GadgetMessageIdentifier>* mb =
                new GadgetContainerMessage<GadgetMessageIdentifier>();

            mb->getObjectPtr()->id = GADGET_MESSAGE_DICOM_WITHNAME;

            mb->cont(mdcm);
            mdcm->cont(mfilename);

            if (m3)
            {
                mfilename->cont(m3);
            }

            int ret = this->controller_->output_ready(mb);

            if ((ret < 0))
            {
                GDEBUG("Failed to return message to controller\n");
                return GADGET_FAIL;
            }

            return GADGET_OK;
        }

    private:
        DcmFileFormat dcmFile;
        std::string seriesIUIDRoot;
        long initialSeriesNumber;
        std::map <unsigned int, std::string> seriesIUIDs;
    };

} /* namespace Gadgetron */

#endif // DICOMFINISHGADGET_H
