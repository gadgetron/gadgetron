/** \file   dicom_ismrmrd_utility.cpp
    \brief  Implement some utility functions to convert ismrmrd image into dicom image
    \author Hui Xue
*/

#include "dicom_ismrmrd_utility.h"
#include <stdio.h>
#include "boost/date_time/gregorian/gregorian.hpp"

namespace Gadgetron
{
    void write_dcm_string(DcmDataset *dataset, DcmTagKey& key, const char* s)
    {
        OFCondition status = dataset->putAndInsertString(key, s);
        if (!status.good())
        {
            GDEBUG("Failed to insert DICOM field (0x%04X,0x%04X) at line %u\n", key.getGroup(), key.getElement(), __LINE__);
            GADGET_THROW("write_dcm_string failed ... ");
        }
    }

    void fill_dicom_image_from_ismrmrd_header(const ISMRMRD::IsmrmrdHeader& h, DcmFileFormat& dcmFile)
    {
        try
        {
            OFCondition status;
            DcmTagKey key;
            long BUFSIZE = 1024;
            std::vector<char> bufVec(BUFSIZE);
            char *buf = &bufVec[0]; 

            // Ensure DICOM dictionary is loaded
            if (!dcmDataDict.isDictionaryLoaded())
            {
                GADGET_THROW("Dictionary not loaded!  Set DCMDICTPATH");
            }

            // -------------------------------------------------

            ISMRMRD::ExperimentalConditions exp_cond = h.experimentalConditions;

            // -------------------------------------------------

            ISMRMRD::SubjectInformation patient_info = *h.subjectInformation;
            if(!h.subjectInformation.is_present() || !h.subjectInformation.get().patientName.is_present())
                patient_info.patientName.set("XXXXXXXX");

            if(!h.subjectInformation.is_present() || !h.subjectInformation.get().patientWeight_kg.is_present())
                patient_info.patientWeight_kg.set(1000);

            if(!h.subjectInformation.is_present() || !h.subjectInformation.get().patientID.is_present())
                patient_info.patientID.set("XXXXXXXX");

            if(!h.subjectInformation.is_present() || !h.subjectInformation.get().patientBirthdate.is_present())
                patient_info.patientBirthdate.set("1900-01-01");

            if(!h.subjectInformation.is_present() || !h.subjectInformation.get().patientGender.is_present())
                patient_info.patientGender.set("o");

            // -------------------------------------------------

            ISMRMRD::StudyInformation study_info = *h.studyInformation;

            if(!h.studyInformation.is_present() || !h.studyInformation.get().studyDate.is_present())
                study_info.studyDate.set("1900-01-01");

            if(!h.studyInformation.is_present() || !h.studyInformation.get().studyTime.is_present())
                study_info.studyTime.set("12:12:12");

            if(!h.studyInformation.is_present() || !h.studyInformation.get().studyID.is_present())
                study_info.studyID.set("XXXXXXXX");

            if(!h.studyInformation.is_present() || !h.studyInformation.get().accessionNumber.is_present())
                study_info.accessionNumber.set(0);

            if(!h.studyInformation.is_present() || !h.studyInformation.get().referringPhysicianName.is_present())
                study_info.referringPhysicianName.set("XXXXXXXX");

            if(!h.studyInformation.is_present() || !h.studyInformation.get().studyDescription.is_present())
                study_info.studyDescription.set("XXXXXXXX");

            if(!h.studyInformation.is_present() || !h.studyInformation.get().studyInstanceUID.is_present())
                study_info.studyInstanceUID.set("XXXXXXXX");

            // -------------------------------------------------

            if (!h.measurementInformation) {
                GWARN("Header missing MeasurementInformation parameters");
            }

            auto meas_info = h.measurementInformation;

            // -------------------------------------------------

            if (!h.acquisitionSystemInformation) {
                GWARN("Header missing AcquisitionSystemInformation parameters");
            }

            auto sys_info = h.acquisitionSystemInformation;

            if (!h.sequenceParameters) {
                GWARN("Header missing SequenceTiming parameters");
            }

            // -------------------------------------------------

            auto seq_info = h.sequenceParameters;

            if (h.encoding.size() == 0) {
                GDEBUG_STREAM("Number of encoding spaces: " << h.encoding.size());
                GADGET_THROW("Encoding description is required");
            }

            // -------------------------------------------------

            ISMRMRD::EncodingSpace e_space = h.encoding[0].encodedSpace;
            ISMRMRD::EncodingSpace r_space = h.encoding[0].reconSpace;
            ISMRMRD::EncodingLimits e_limits = h.encoding[0].encodingLimits;

            // -------------------------------------------------

            DcmDataset *dataset = dcmFile.getDataset();
            DcmMetaInfo *metainfo = dcmFile.getMetaInfo();

            // Set the Application Entity Title in the DICOM Meta Info section
            // The rest of the Meta Info will be automatically populated by DCMTK
            if (sys_info && sys_info->stationName) {
                status = metainfo->putAndInsertString(DcmTagKey(0x0002, 0x0016),
                    sys_info.get().stationName->c_str());
                if (!status.good()) {
                    GADGET_THROW("Failed to set AET in MetaInfo");
                }
            }
            else {
                status = metainfo->putAndInsertString(DcmTagKey(0x0002, 0x0016), "none");
                if (!status.good()) {
                    GADGET_THROW("Failed to set AET in MetaInfo");
                }
            }

            // Group Length
            key.set(0x0008, 0x0000);
            status = dataset->insertEmptyElement(key);
            if (status.bad()) {
                GADGET_THROW("Failed to write 0x0008 Group Length");
            }

            // Specific Character Set
            key.set(0x0008, 0x0005);
            write_dcm_string(dataset, key, "ISO_IR 192");

            // Image Type
            // ORIGINAL or DERIVED describes origin of pixel data
            // PRIMARY or SECONDARY describes image creation time (during or after exam)
            // OTHER, etc. are implementation-specific
            key.set(0x0008, 0x0008);
            write_dcm_string(dataset, key, "ORIGINAL\\PRIMARY\\OTHER");

            // SOPClassUID
            key.set(0x0008, 0x0016);
            write_dcm_string(dataset, key, UID_MRImageStorage);

            // Study Date
            if (study_info.studyDate) {
                key.set(0x0008, 0x0020);
                std::string d(study_info.studyDate.get());
                d.erase(std::remove(d.begin(), d.end(), '-'), d.end());     // erase all occurrences of '-'
                write_dcm_string(dataset, key, d.c_str());
            }

            // Series, Acquisition, Content Date
            if (meas_info && meas_info->seriesDate) {
                key.set(0x0008, 0x0021);
                std::string d(meas_info.get().seriesDate.get());
                d.erase(std::remove(d.begin(), d.end(), '-'), d.end());
                write_dcm_string(dataset, key, d.c_str());

                key.set(0x0008, 0x0022);
                write_dcm_string(dataset, key, d.c_str());

                key.set(0x0008, 0x0023);
                write_dcm_string(dataset, key, d.c_str());
            }

            // Study Time
            if (study_info.studyTime) {
                key.set(0x0008, 0x0030);
                std::string t(study_info.studyTime.get());
                t.erase(std::remove(t.begin(), t.end(), ':'), t.end());
                write_dcm_string(dataset, key, t.c_str());
            }

            // Series, Acquisition, Content Time
            if (meas_info->seriesTime) {
                key.set(0x0008, 0x0031);
                std::string t(meas_info.get().seriesTime.get());
                t.erase(std::remove(t.begin(), t.end(), ':'), t.end());
                write_dcm_string(dataset, key, t.c_str());

                key.set(0x0008, 0x0032);
                write_dcm_string(dataset, key, t.c_str());

                key.set(0x0008, 0x0033);
                write_dcm_string(dataset, key, t.c_str());
            }

            // Accession Number
            key.set(0x0008, 0x0050);
            if (study_info.accessionNumber) {
                snprintf(buf, BUFSIZE, "%ld", *study_info.accessionNumber);
                write_dcm_string(dataset, key, buf);
            }
            else {
                const char* p = NULL;
                write_dcm_string(dataset, key, p);
            }

            // Modality
            // TODO: this is hardcoded!!
            key.set(0x0008, 0x0060);
            write_dcm_string(dataset, key, "MR");

            // Manufacturer
            key.set(0x0008, 0x0070);
            if (sys_info && sys_info->systemVendor) {
                write_dcm_string(dataset, key, sys_info.get().systemVendor->c_str());
            }
            else {
                write_dcm_string(dataset, key, "UNKNOWN");
            }

            // Institution Name
            key.set(0x0008, 0x0080);
            if (sys_info && sys_info->institutionName) {
                write_dcm_string(dataset, key, sys_info.get().institutionName->c_str());
            }
            else {
                write_dcm_string(dataset, key, "UNKNOWN");
            }

            // Referring Physician's Name
            key.set(0x0008, 0x0090);
            if (study_info.referringPhysicianName) {
                write_dcm_string(dataset, key, study_info.referringPhysicianName->c_str());
            }
            else {
                write_dcm_string(dataset, key, "");
            }

            // Station Name
            key.set(0x0008, 0x1010);
            if (sys_info && sys_info->stationName) {
                write_dcm_string(dataset, key, sys_info.get().stationName->c_str());
            }
            else {
                write_dcm_string(dataset, key, "");
            }

            // Study Description
            key.set(0x0008, 0x1030);
            if (study_info.studyDescription) {
                write_dcm_string(dataset, key, study_info.studyDescription->c_str());
            }
            else {
                write_dcm_string(dataset, key, "");
            }

            // Series Description
            key.set(0x0008, 0x103E);
            if (meas_info->seriesDescription) {
                write_dcm_string(dataset, key, meas_info.get().seriesDescription->c_str());
            }
            else {
                write_dcm_string(dataset, key, "");
            }

            // Manufacturer's Model Name
            key.set(0x0008, 0x1090);
            if (sys_info && sys_info->systemModel) {
                write_dcm_string(dataset, key, sys_info.get().systemModel->c_str());
            }
            else {
                write_dcm_string(dataset, key, "");
            }

            // Referenced SOP Instance UIDs
            std::vector<ISMRMRD::ReferencedImageSequence> refs(meas_info->referencedImageSequence);
            if (refs.size() > 0) {
                DcmItem *ref_sequence;
                std::vector<ISMRMRD::ReferencedImageSequence>::iterator it;
                for (it = refs.begin(); it != refs.end(); ++it) {
                    std::string ref_uid(it->referencedSOPInstanceUID);
                    if (ref_uid.length() > 0) {   // Only write non-empty strings
                        if (dataset->findOrCreateSequenceItem(key, ref_sequence, -2).good()) {
                            // Write the Referenced SOPClassUID (MRImageStorage)
                            key.set(0x0008, 0x1150);
                            ((DcmDataset *)ref_sequence)->putAndInsertString(key, UID_MRImageStorage);
                            // Write the Referenced SOPInstanceUID
                            key.set(0x0008, 0x1155);
                            ((DcmDataset *)ref_sequence)->putAndInsertString(key, ref_uid.c_str());
                        }
                    }
                }
            }

            // Group Length
            key.set(0x0010, 0x0000);
            status = dataset->insertEmptyElement(key);
            if (!status.good()) {
                GADGET_THROW("Failed to write 0x0010 Group Length");
            }

            // Patient Name
            key.set(0x0010, 0x0010);
            if (patient_info.patientName) {
                write_dcm_string(dataset, key, patient_info.patientName->c_str());
            }
            else {
                write_dcm_string(dataset, key, "None");
            }

            // Patient ID
            key.set(0x0010, 0x0020);
            if (patient_info.patientID) {
                write_dcm_string(dataset, key, patient_info.patientID->c_str());
            }
            else {
                write_dcm_string(dataset, key, "0");
            }

            // Patient Birthdate
            key.set(0x0010, 0x0030);
            if (patient_info.patientBirthdate) {
                std::string d(patient_info.patientBirthdate.get());
                d.erase(std::remove(d.begin(), d.end(), '-'), d.end());
                write_dcm_string(dataset, key, d.c_str());
            }
            else {
                status = dataset->insertEmptyElement(key);
            }

            // Patient Sex
            key.set(0x0010, 0x0040);
            if (patient_info.patientGender) {
                if (*patient_info.patientGender == "O") {
                    status = dataset->insertEmptyElement(key);
                }
                else {
                    write_dcm_string(dataset, key, patient_info.patientGender->c_str());
                }
            }
            else {
                write_dcm_string(dataset, key, "");
            }

            // Patient Age
            key.set(0x0010, 0x1010);
            if (patient_info.patientBirthdate && meas_info->seriesDate) {
                boost::gregorian::date bday(boost::gregorian::from_simple_string(patient_info.patientBirthdate.get()));
                boost::gregorian::date seriesDate(boost::gregorian::from_simple_string(meas_info.get().seriesDate.get()));

                boost::gregorian::days age = seriesDate - bday;

                long age_in_years = age.days() / 365;

                snprintf(buf, BUFSIZE, "%03ldY", age_in_years);
                write_dcm_string(dataset, key, buf);
            }
            else {
                write_dcm_string(dataset, key, "000Y");
            }

            // Patient Weight
            key.set(0x0010, 0x1030);
            if (patient_info.patientWeight_kg) {
                snprintf(buf, BUFSIZE, "%f", *patient_info.patientWeight_kg);
                write_dcm_string(dataset, key, buf);
            }
            else {
                write_dcm_string(dataset, key, "0.0");
            }

            // Group Length
            key.set(0x0018, 0x0000);
            status = dataset->insertEmptyElement(key);
            if (!status.good()) {
                GADGET_THROW("Failed to write 0x0018 Group Length");
            }

            // Scanning Sequence, Sequence Variant, Scan Options, Acquisition Type
            std::string scanningSequence("RM");
            std::string sequenceVariant("NONE");
            std::string scanOptions("NONE");
            std::string mrAcquisitionType("2D");
            if (h.userParameters) {
                ISMRMRD::UserParameters user_params = h.userParameters.get();
                std::vector<ISMRMRD::UserParameterString> strings = user_params.userParameterString;
                std::vector<ISMRMRD::UserParameterString>::iterator it;

                for (it = strings.begin(); it != strings.end(); ++it) {
                    if (it->name == "scanningSequence") {
                        scanningSequence = it->value;
                    }
                    else if (it->name == "sequenceVariant") {
                        sequenceVariant = it->value;
                    }
                    else if (it->name == "scanOptions") {
                        scanOptions = it->value;
                    }
                    else if (it->name == "mrAcquisitionType") {
                        mrAcquisitionType = it->value;
                    }
                }
            }
            key.set(0x0018, 0x0020);
            write_dcm_string(dataset, key, scanningSequence.c_str());
            key.set(0x0018, 0x0021);
            write_dcm_string(dataset, key, sequenceVariant.c_str());
            key.set(0x0018, 0x0022);
            write_dcm_string(dataset, key, scanOptions.c_str());
            key.set(0x0018, 0x0023);
            write_dcm_string(dataset, key, mrAcquisitionType.c_str());

            // Angio Flag
            // TODO: hardcoded
            key.set(0x0018, 0x0025);
            write_dcm_string(dataset, key, "N");

            // Slice Thickness
            // This will need updated if the "reconSpace.fieldOfView_mm.z" field
            // is changed in the ISMRMRD populating code (client)
            if (r_space.matrixSize.z == 0) r_space.matrixSize.z = 1;
            key.set(0x0018, 0x0050);
            snprintf(buf, BUFSIZE, "%f", r_space.fieldOfView_mm.z / r_space.matrixSize.z);
            write_dcm_string(dataset, key, buf);

            // Spacing Between Slices
            key.set(0x0018, 0x0088);
            snprintf(buf, BUFSIZE, "%f", r_space.fieldOfView_mm.z);
            write_dcm_string(dataset, key, buf);

            // Repetition Time
            if (seq_info && seq_info.get().TR.is_present() && seq_info.get().TR.get().size() > 0)
            {
                key.set(0x0018, 0x0080);
                snprintf(buf, BUFSIZE, "%f", seq_info.get().TR.get().front());
                write_dcm_string(dataset, key, buf);
            }

            // Echo Time
            if (seq_info && seq_info.get().TE.is_present() && seq_info.get().TE.get().size() > 0)
            {
                key.set(0x0018, 0x0081);
                snprintf(buf, BUFSIZE, "%f", seq_info.get().TE.get().front());
                write_dcm_string(dataset, key, buf);
            }

            // Inversion Time
            if (seq_info && seq_info.get().TI.is_present() && seq_info.get().TI.get().size()>0)
            {
                key.set(0x0018, 0x0082);
                snprintf(buf, BUFSIZE, "%f", seq_info.get().TI.get().front());
                write_dcm_string(dataset, key, buf);
            }

            // Flip Angle
            if (seq_info && seq_info.get().flipAngle_deg.is_present() && seq_info.get().flipAngle_deg.get().size()>0)
            {
                key.set(0x0018, 0x1314);
                snprintf(buf, BUFSIZE, "%ld", (long)seq_info.get().flipAngle_deg.get().front());
                write_dcm_string(dataset, key, buf);
            }

            // Imaging Frequency in tenths of MHz ???
            key.set(0x0018, 0x0084);
            snprintf(buf, BUFSIZE, "%f", (float)exp_cond.H1resonanceFrequency_Hz / 10000000.);
            write_dcm_string(dataset, key, buf);

            // Magnetic Field Strength (T)
            key.set(0x0018, 0x0087);
            if (sys_info && sys_info->systemFieldStrength_T) {
                snprintf(buf, BUFSIZE, "%f", *sys_info->systemFieldStrength_T);
                write_dcm_string(dataset, key, buf);
            }
            else {
                write_dcm_string(dataset, key, "3.0");
            }

            key.set(0x0018, 0x0091);
            // Echo Train Length
            if (h.encoding[0].echoTrainLength) {
                snprintf(buf, BUFSIZE, "%ld", (long)*h.encoding[0].echoTrainLength);
                write_dcm_string(dataset, key, buf);
            }
            else {
                write_dcm_string(dataset, key, "1");
            }

            // Percent Sampling
            // TODO: hardcoded
            key.set(0x0018, 0x0093);
            write_dcm_string(dataset, key, "100");

            // Percent Phase FOV
            // TODO: hardcoded
            key.set(0x0018, 0x0094);
            write_dcm_string(dataset, key, "100");

            // Protocol Name
            if (meas_info && meas_info->protocolName) {
                key.set(0x0018, 0x1030);
                write_dcm_string(dataset, key, meas_info.get().protocolName.get().c_str());
            }
            else {
                write_dcm_string(dataset, key, "");
            }

            // Trigger Time - TODO: use Image Meta Data
            key.set(0x0018, 0x1060);
            write_dcm_string(dataset, key, "0.0");

            // Reconstruction Diameter (FOV) - TODO: ?
            key.set(0x0018, 0x1100);

            // Frequency Encoding Direction - TODO: use Image Meta Data
            key.set(0x0018, 0x1312);
            write_dcm_string(dataset, key, "ROW");

            if (meas_info) {
                // Patient Position
                key.set(0x0018, 0x5100);
                write_dcm_string(dataset, key, meas_info.get().patientPosition.c_str());
            }
            /****************************************/
            // Group Length
            key.set(0x0020, 0x0000);
            status = dataset->insertEmptyElement(key);
            if (!status.good()) {
                GADGET_THROW("Failed to write 0x0020 Group Length");
            }

            // Study Instance UID
            key.set(0x0020, 0x000D);
            if (study_info.studyInstanceUID) {
                write_dcm_string(dataset, key, study_info.studyInstanceUID->c_str());
            }

            // Study ID
            if (study_info.studyID) {
                key.set(0x0020, 0x0010);
                write_dcm_string(dataset, key, study_info.studyID->c_str());
            }
            else {
                write_dcm_string(dataset, key, "0");
            }

            // Frame of Reference UID
            if (meas_info && meas_info->frameOfReferenceUID) {
                key.set(0x0020, 0x0052);
                write_dcm_string(dataset, key, meas_info.get().frameOfReferenceUID->c_str());
            }

            /****************************************/
            // Group Length
            key.set(0x0028, 0x0000);
            status = dataset->insertEmptyElement(key);
            if (!status.good()) {
                GADGET_THROW("Failed to write 0x0028 Group Length");
            }

            // Samples Per Pixel
            key.set(0x0028, 0x0002);
            // TODO: hardcoded
            write_dcm_string(dataset, key, "1");

            // Photometric Interpretation
            key.set(0x0028, 0x0004);
            // TODO: hardcoded
            write_dcm_string(dataset, key, "MONOCHROME2");

            // Pixel Spacing (Array of len 2)
            key.set(0x0028, 0x0030);
            float pixel_spacing_X = r_space.fieldOfView_mm.x / r_space.matrixSize.x;
            float pixel_spacing_Y = r_space.fieldOfView_mm.y / r_space.matrixSize.y;
            snprintf(buf, BUFSIZE, "%.3f\\%.3f", pixel_spacing_X, pixel_spacing_Y);
            write_dcm_string(dataset, key, buf);

            // Bits Allocated
            key.set(0x0028, 0x0100);
            write_dcm_string(dataset, key, "16");
            // Bits Stored
            key.set(0x0028, 0x0101);
            write_dcm_string(dataset, key, "16");
            // High Bit
            key.set(0x0028, 0x0102);
            write_dcm_string(dataset, key, "15");
            // Pixel Representation
            key.set(0x0028, 0x0103);
            write_dcm_string(dataset, key, "1");
        }
        catch(...)
        {
            GADGET_THROW("Exceptions happened in fill_dicom_image_from_ismrmrd_header(...) ... ");
        }
    }

    template<typename T> 
    void write_ismrmd_image_into_dicom(const ISMRMRD::ImageHeader& m1, const hoNDArray<T>& m2, std::string& seriesIUID, DcmFileFormat& dcmFile)
    {
        try
        {
            hoNDArray< int16_t > data(m2.dimensions());

            const T *src = m2.get_data_ptr();
            auto dst = data.get_data_ptr();

            T min_pix_val, max_pix_val, sum_pix_val = 0;
            if (m2.get_number_of_elements() > 0)
            {
                min_pix_val = src[0];
                max_pix_val = src[0];
            }

            for (unsigned long i = 0; i < m2.get_number_of_elements(); i++)
            {
                T pix_val = src[i];
                // search for minimum and maximum pixel values
                if (pix_val < min_pix_val) min_pix_val = pix_val;
                if (pix_val > max_pix_val) max_pix_val = pix_val;
                sum_pix_val += pix_val / 4; // scale by 25% to avoid overflow

                dst[i] = static_cast<int16_t>(pix_val);
            }
            T mean_pix_val = (T)((sum_pix_val * 4) / (T)data.get_number_of_elements());

            unsigned int BUFSIZE = 1024;
            std::vector<char> bufVec(BUFSIZE);
            char *buf = &bufVec[0];

            OFCondition status;
            DcmTagKey key;
            DcmDataset *dataset = dcmFile.getDataset();

            // Echo Number
            // TODO: it is often the case the img->contrast is not properly set
            // likely due to the allocated ISMRMRD::ImageHeader being uninitialized
            key.set(0x0018, 0x0086);
            snprintf(buf, BUFSIZE, "%d", m1.contrast);
            write_dcm_string(dataset, key, buf);

            // Acquisition Matrix ... Image Dimensions
            // Defined as: [frequency rows, frequency columns, phase rows, phase columns]
            // But at this point in the gadget I don't know the frequency encode direction
            key.set(0x0018, 0x1310);
            uint16_t im_dim[4] = { 0, 0, 0, 0 };
            im_dim[1] = m1.matrix_size[0];
            im_dim[2] = m1.matrix_size[1];

            status = dataset->putAndInsertUint16Array(key, im_dim, 4);
            if (!status.good())
            {
                GADGET_THROW("Failed to stuff image dimensions");
            }

            // Series Number
            // Only write a number if the image_series_index is positive and non-zero
            key.set(0x0020, 0x0011);
            snprintf(buf, BUFSIZE, "%ld", (long int)m1.image_series_index);
            write_dcm_string(dataset, key, buf);

            // Image Number
            key.set(0x0020, 0x0013);
            snprintf(buf, BUFSIZE, "%d", m1.image_index + 1);
            write_dcm_string(dataset, key, buf);

            // Image Position (Patient)
            float corner[3];

            corner[0] = m1.position[0] -
                (m1.field_of_view[0] / 2.0f) * m1.read_dir[0] -
                (m1.field_of_view[1] / 2.0f) * m1.phase_dir[0];
            corner[1] = m1.position[1] -
                (m1.field_of_view[0] / 2.0f) * m1.read_dir[1] -
                (m1.field_of_view[1] / 2.0f) * m1.phase_dir[1];
            corner[2] = m1.position[2] -
                (m1.field_of_view[0] / 2.0f) * m1.read_dir[2] -
                (m1.field_of_view[1] / 2.0f) * m1.phase_dir[2];

            key.set(0x0020, 0x0032);
            snprintf(buf, BUFSIZE, "%.4f\\%.4f\\%.4f", corner[0], corner[1], corner[2]);
            write_dcm_string(dataset, key, buf);

            // Image Orientation
            // read_dir, phase_dir, and slice_dir were calculated in
            // a DICOM/patient coordinate system, so just plug them in
            key.set(0x0020, 0x0037);
            snprintf(buf, BUFSIZE, "%.4f\\%.4f\\%.4f\\%.4f\\%.4f\\%.4f",
                m1.read_dir[0], m1.read_dir[1], m1.read_dir[2],
                m1.phase_dir[0], m1.phase_dir[1], m1.phase_dir[2]);
            write_dcm_string(dataset, key, buf);

            // Slice Location
            key.set(0x0020, 0x1041);
            snprintf(buf, BUFSIZE, "%f", m1.position[2]);
            write_dcm_string(dataset, key, buf);

            // Columns
            key.set(0x0028, 0x0010);
            snprintf(buf, BUFSIZE, "%d", m1.matrix_size[1]);
            write_dcm_string(dataset, key, buf);

            // Rows
            key.set(0x0028, 0x0011);
            snprintf(buf, BUFSIZE, "%d", m1.matrix_size[0]);
            write_dcm_string(dataset, key, buf);

            //Number of frames
            if (m1.matrix_size[2] > 1){ //Only write if we have more than 1 frame
                key.set(0x0028,0x0008);
                snprintf(buf,BUFSIZE,"%d",m1.matrix_size[2]);
                write_dcm_string(dataset, key,buf);
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
            snprintf(buf, BUFSIZE, "%d", window_center);
            write_dcm_string(dataset, key, buf);

            // Window Width
            key.set(0x0028, 0x1051);
            snprintf(buf, BUFSIZE, "%d", window_width);
            write_dcm_string(dataset, key, buf);

            // ACR_NEMA_2C_VariablePixelDataGroupLength
            key.set(0x7fe0, 0x0000);
            status = dataset->insertEmptyElement(key);
            if (!status.good()) {
                GADGET_THROW("Failed to write 0x7fe0 Group Length");
            }

            // Pixel Data
            if ((unsigned long)m1.matrix_size[0] * (unsigned long)m1.matrix_size[1]*(unsigned long)m1.matrix_size[2] !=
                data.get_number_of_elements()) {
                GADGET_THROW("Mismatch in image dimensions and available data");
            }
            key.set(0x7fe0, 0x0010);
            status = dataset->putAndInsertUint16Array(key, (unsigned short *)data.get_data_ptr(), (unsigned long)data.get_number_of_elements());
            if (!status.good()) {
                GADGET_THROW("Failed to stuff Pixel Data");
            }

            // Series Instance UID = generated here
            key.set(0x0020, 0x000E);
            write_dcm_string(dataset, key, seriesIUID.c_str());

            // At a minimum, to put the DICOM image back into the database,
            // you must change the SOPInstanceUID.
            key.set(0x0008, 0x0018);        // SOPInstanceUID
            const char *root = "1.2.840.113619.2.156";
            char newuid[65];
            dcmGenerateUniqueIdentifier(newuid, root);
            write_dcm_string(dataset, key, newuid);
        }
        catch(...)
        {
            GADGET_THROW("Exceptions happened in write_ismrmd_image_into_dicom(...) ... ");
        }
    }

    template EXPORTGADGETSDICOM void write_ismrmd_image_into_dicom(const ISMRMRD::ImageHeader& m1, const hoNDArray<short>& m2, std::string& seriesIUID, DcmFileFormat& dcmFile);
    template EXPORTGADGETSDICOM void write_ismrmd_image_into_dicom(const ISMRMRD::ImageHeader& m1, const hoNDArray<unsigned short>& m2, std::string& seriesIUID, DcmFileFormat& dcmFile);
    template EXPORTGADGETSDICOM void write_ismrmd_image_into_dicom(const ISMRMRD::ImageHeader& m1, const hoNDArray<int>& m2, std::string& seriesIUID, DcmFileFormat& dcmFile);
    template EXPORTGADGETSDICOM void write_ismrmd_image_into_dicom(const ISMRMRD::ImageHeader& m1, const hoNDArray<unsigned int>& m2, std::string& seriesIUID, DcmFileFormat& dcmFile);
    template EXPORTGADGETSDICOM void write_ismrmd_image_into_dicom(const ISMRMRD::ImageHeader& m1, const hoNDArray<float>& m2, std::string& seriesIUID, DcmFileFormat& dcmFile);
    template EXPORTGADGETSDICOM void write_ismrmd_image_into_dicom(const ISMRMRD::ImageHeader& m1, const hoNDArray<double>& m2, std::string& seriesIUID, DcmFileFormat& dcmFile);

    template<typename T> 
    void write_ismrmd_image_into_dicom(const ISMRMRD::ImageHeader& m1, const hoNDArray<T>& m2, ISMRMRD::IsmrmrdHeader& h, ISMRMRD::MetaContainer& attrib, std::string& seriesIUID, DcmFileFormat& dcmFile)
    {
        try
        {
            Gadgetron::write_ismrmd_image_into_dicom(m1, m2, seriesIUID, dcmFile);

            DcmDataset *dataset = dcmFile.getDataset();
            DcmTagKey key;

            // ---------------------------------------------------------------------
            // check the attrib and filling corresponding dicom fields
            // ---------------------------------------------------------------------

            long BUFSIZE = 1024;
            std::vector<char> bufVec(BUFSIZE);
            char *buf = &bufVec[0]; 

            size_t n;

            // ----------------------------------
            // image comment
            // ----------------------------------
            size_t N = attrib.length(GADGETRON_IMAGECOMMENT);
            if(N>0)
            {
                std::string str(std::string(attrib.as_str(GADGETRON_IMAGECOMMENT, 0)));
                for (n=1; n<N; n++)
                {
                    str.append("_");
                    str.append( std::string(attrib.as_str(GADGETRON_IMAGECOMMENT, n)) );
                }

                key.set(0x0020, 0x4000);
                write_dcm_string(dataset, key, str.c_str());
            }

            // ----------------------------------
            // sequence description
            // ----------------------------------
            N = attrib.length(GADGETRON_SEQUENCEDESCRIPTION);
            if(N>0)
            {
                std::string str(std::string(attrib.as_str(GADGETRON_SEQUENCEDESCRIPTION, 0)));
                for (n=1; n<N; n++)
                {
                    str.append("_");
                    str.append( std::string(attrib.as_str(GADGETRON_SEQUENCEDESCRIPTION, n)) );
                }

                key.set(0x0008, 0x103E);
                write_dcm_string(dataset, key, str.c_str());
            }

            // ----------------------------------
            // TR
            // ----------------------------------
            if(h.sequenceParameters.is_present())
            {
                if(h.sequenceParameters.get().TR.is_present())
                {
                    float v = h.sequenceParameters.get().TR.get()[0];

                    key.set(0x0018, 0x0080);
                    snprintf(buf, BUFSIZE, "%f", v);
                    write_dcm_string(dataset, key, buf);
                }
            }

            // ----------------------------------
            // TE
            // ----------------------------------
            if ( attrib.length(GADGETRON_IMAGE_ECHOTIME) > 0 )
            {
                double v = attrib.as_double(GADGETRON_IMAGE_ECHOTIME, 0);

                key.set(0x0018, 0x0081);
                snprintf(buf, BUFSIZE, "%f", v);
                write_dcm_string(dataset, key, buf);
            }
            else
            {
                if(h.sequenceParameters.is_present())
                {
                    if(h.sequenceParameters.get().TE.is_present())
                    {
                        float v = h.sequenceParameters.get().TE.get()[0];

                        key.set(0x0018, 0x0081);
                        snprintf(buf, BUFSIZE, "%f", v);
                        write_dcm_string(dataset, key, buf);
                    }
                }
            }

            // ----------------------------------
            // Trigger Time
            // ----------------------------------
            uint32_t ticks = m1.physiology_time_stamp[0];
            double millisec = ticks * 2.5;
            if ( millisec > 0 && millisec<60000 )
            {
                key.set(0x0018, 0x1060);
                snprintf(buf, BUFSIZE, "%f", millisec);
                write_dcm_string(dataset, key, buf);
            }

            // ----------------------------------
            // TI
            // ----------------------------------
            if ( attrib.length(GADGETRON_IMAGE_INVERSIONTIME) > 0 )
            {
                double v = attrib.as_double(GADGETRON_IMAGE_INVERSIONTIME, 0);

                key.set(0x0018, 0x0082);
                snprintf(buf, BUFSIZE, "%f", v);
                write_dcm_string(dataset, key, buf);
            }
            else
            {
                if(h.sequenceParameters.is_present())
                {
                    if(h.sequenceParameters.get().TI.is_present())
                    {
                        float v = h.sequenceParameters.get().TI.get()[0];

                        key.set(0x0018, 0x0082);
                        snprintf(buf, BUFSIZE, "%f", v);
                        write_dcm_string(dataset, key, buf);
                    }
                }
            }

            // ----------------------------------
            // WindowCenter, WindowWidth
            // ----------------------------------
            if ( attrib.length(GADGETRON_IMAGE_WINDOWCENTER) > 0 )
            {
                double v = attrib.as_double(GADGETRON_IMAGE_WINDOWCENTER, 0);

                key.set(0x0028, 0x1050);
                snprintf(buf, BUFSIZE, "%d", (int)v);
                write_dcm_string(dataset, key, buf);
            }

            if ( attrib.length(GADGETRON_IMAGE_WINDOWWIDTH) > 0 )
            {
                double v = attrib.as_double(GADGETRON_IMAGE_WINDOWWIDTH, 0);

                key.set(0x0028, 0x1051);
                snprintf(buf, BUFSIZE, "%d", (int)v);
                write_dcm_string(dataset, key, buf);
            }
        }
        catch(...)
        {
            GADGET_THROW("Exceptions happened in write_ismrmd_image_into_dicom(attrib) ... ");
        }
    }

    template EXPORTGADGETSDICOM void write_ismrmd_image_into_dicom(const ISMRMRD::ImageHeader& m1, const hoNDArray<short>& m2, ISMRMRD::IsmrmrdHeader& h, ISMRMRD::MetaContainer& attrib, std::string& seriesIUID, DcmFileFormat& dcmFile);
    template EXPORTGADGETSDICOM void write_ismrmd_image_into_dicom(const ISMRMRD::ImageHeader& m1, const hoNDArray<unsigned short>& m2, ISMRMRD::IsmrmrdHeader& h, ISMRMRD::MetaContainer& attrib, std::string& seriesIUID, DcmFileFormat& dcmFile);
    template EXPORTGADGETSDICOM void write_ismrmd_image_into_dicom(const ISMRMRD::ImageHeader& m1, const hoNDArray<int>& m2, ISMRMRD::IsmrmrdHeader& h, ISMRMRD::MetaContainer& attrib, std::string& seriesIUID, DcmFileFormat& dcmFile);
    template EXPORTGADGETSDICOM void write_ismrmd_image_into_dicom(const ISMRMRD::ImageHeader& m1, const hoNDArray<unsigned int>& m2, ISMRMRD::IsmrmrdHeader& h, ISMRMRD::MetaContainer& attrib, std::string& seriesIUID, DcmFileFormat& dcmFile);
    template EXPORTGADGETSDICOM void write_ismrmd_image_into_dicom(const ISMRMRD::ImageHeader& m1, const hoNDArray<float>& m2, ISMRMRD::IsmrmrdHeader& h, ISMRMRD::MetaContainer& attrib, std::string& seriesIUID, DcmFileFormat& dcmFile);
    template EXPORTGADGETSDICOM void write_ismrmd_image_into_dicom(const ISMRMRD::ImageHeader& m1, const hoNDArray<double>& m2, ISMRMRD::IsmrmrdHeader& h, ISMRMRD::MetaContainer& attrib, std::string& seriesIUID, DcmFileFormat& dcmFile);
}
