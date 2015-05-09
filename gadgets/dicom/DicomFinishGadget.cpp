#include <vector>
#include "boost/date_time/gregorian/gregorian.hpp"

#include "DicomFinishGadget.h"
#include "ismrmrd/xml.h"

namespace Gadgetron {

    int DicomFinishGadget::process_config(ACE_Message_Block* mb)
    {
        OFCondition status;
        DcmTagKey key;
        long BUFSIZE = 1024;
        char *buf = new char[BUFSIZE];  // used for writing numbers as strings in DCMTK

        ISMRMRD::IsmrmrdHeader h;
        deserialize(mb->rd_ptr(), h);

        // Ensure DICOM dictionary is loaded
        if (!dcmDataDict.isDictionaryLoaded()) {
            GDEBUG("Dictionary not loaded!  Set DCMDICTPATH\n");
            return GADGET_FAIL;
        }

        ISMRMRD::ExperimentalConditions exp_cond = h.experimentalConditions;

        ISMRMRD::SubjectInformation patient_info = *h.subjectInformation;
        if (!h.subjectInformation) {
            GWARN("Header missing SubjectInformation parameters\n");

            patient_info.patientName.set("XXXXXXXX");
            patient_info.patientWeight_kg.set(1000);
            patient_info.patientID.set("XXXXXXXX");
            patient_info.patientBirthdate.set("19000101");
            patient_info.patientGender.set("o");
        }

        ISMRMRD::StudyInformation study_info = *h.studyInformation;
        if (!h.studyInformation) {
            GWARN("Header missing StudyInformation parameters\n");

            study_info.studyDate.set("19000101");
            study_info.studyTime.set("121212");
            study_info.studyID.set("XXXXXXXX");
            study_info.accessionNumber.set(0);
            study_info.referringPhysicianName.set("XXXXXXXX");
            study_info.studyDescription.set("XXXXXXXX");
            study_info.studyInstanceUID.set("XXXXXXXX");
        }

        if (!h.measurementInformation) {
            GDEBUG("Header missing MeasurementInformation parameters\n");
            return GADGET_FAIL;
        }

        ISMRMRD::MeasurementInformation meas_info = *h.measurementInformation;

        if (!h.acquisitionSystemInformation) {
            GDEBUG("Header missing AcquisitionSystemInformation parameters\n");
            return GADGET_FAIL;
        }

        ISMRMRD::AcquisitionSystemInformation sys_info = *h.acquisitionSystemInformation;

        if (!h.sequenceParameters) {
            GDEBUG("Header missing SequenceTiming parameters\n");
            return GADGET_FAIL;
        }

        ISMRMRD::SequenceParameters seq_info = *h.sequenceParameters;

        if (h.encoding.size() == 0) {
            GDEBUG("Number of encoding spaces: %d\n", h.encoding.size());
            GDEBUG("This Gadget needs an encoding description\n");
            return GADGET_FAIL;
        }


        ISMRMRD::EncodingSpace e_space = h.encoding[0].encodedSpace;
        ISMRMRD::EncodingSpace r_space = h.encoding[0].reconSpace;
        ISMRMRD::EncodingLimits e_limits = h.encoding[0].encodingLimits;

        DcmDataset *dataset = dcmFile.getDataset();
        DcmMetaInfo *metainfo = dcmFile.getMetaInfo();

        // Store initial Series Number for later
        if (meas_info.initialSeriesNumber) {
            this->initialSeriesNumber = (long)*meas_info.initialSeriesNumber;
        }
        else {
            this->initialSeriesNumber = 0;
        }


        // Set the Application Entity Title in the DICOM Meta Info section
        // The rest of the Meta Info will be automatically populated by DCMTK
        if (sys_info.stationName) {
            status = metainfo->putAndInsertString(DcmTagKey(0x0002, 0x0016),
                sys_info.stationName->c_str());
            if (!status.good()) {
                GDEBUG("Failed to set AET in MetaInfo\n");
                return GADGET_FAIL;
            }
        }
        else {
            status = metainfo->putAndInsertString(DcmTagKey(0x0002, 0x0016), "none");
            if (!status.good()) {
                GDEBUG("Failed to set AET in MetaInfo\n");
                return GADGET_FAIL;
            }
        }

        // Group Length
        key.set(0x0008, 0x0000);
        status = dataset->insertEmptyElement(key);
        if (status.bad()) {
            GDEBUG("Failed to write 0x0008 Group Length\n");
            return GADGET_FAIL;
        }

        // Specific Character Set
        key.set(0x0008, 0x0005);
        WRITE_DCM_STRING(key, "ISO_IR 100");

        // Image Type
        // ORIGINAL or DERIVED describes origin of pixel data
        // PRIMARY or SECONDARY describes image creation time (during or after exam)
        // OTHER, etc. are implementation-specific
        key.set(0x0008, 0x0008);
        WRITE_DCM_STRING(key, "ORIGINAL\\PRIMARY\\OTHER");

        // SOPClassUID
        key.set(0x0008, 0x0016);
        WRITE_DCM_STRING(key, UID_MRImageStorage);

        // Study Date
        if (study_info.studyDate) {
            key.set(0x0008, 0x0020);
            std::string d(study_info.studyDate.get());
            d.erase(std::remove(d.begin(), d.end(), '-'), d.end());     // erase all occurrences of '-'
            WRITE_DCM_STRING(key, d.c_str());
        }

        // Series, Acquisition, Content Date
        if (meas_info.seriesDate) {
            key.set(0x0008, 0x0021);
            std::string d(meas_info.seriesDate.get());
            d.erase(std::remove(d.begin(), d.end(), '-'), d.end());
            WRITE_DCM_STRING(key, d.c_str());

            key.set(0x0008, 0x0022);
            WRITE_DCM_STRING(key, d.c_str());

            key.set(0x0008, 0x0023);
            WRITE_DCM_STRING(key, d.c_str());
        }

        // Study Time
        if (study_info.studyTime) {
            key.set(0x0008, 0x0030);
            std::string t(study_info.studyTime.get());
            t.erase(std::remove(t.begin(), t.end(), ':'), t.end());
            WRITE_DCM_STRING(key, t.c_str());
        }

        // Series, Acquisition, Content Time
        if (meas_info.seriesTime) {
            key.set(0x0008, 0x0031);
            std::string t(meas_info.seriesTime.get());
            t.erase(std::remove(t.begin(), t.end(), ':'), t.end());
            WRITE_DCM_STRING(key, t.c_str());

            key.set(0x0008, 0x0032);
            WRITE_DCM_STRING(key, t.c_str());

            key.set(0x0008, 0x0033);
            WRITE_DCM_STRING(key, t.c_str());
        }

        // Accession Number
        key.set(0x0008, 0x0050);
        if (study_info.accessionNumber) {
            ACE_OS::snprintf(buf, BUFSIZE, "%ld", *study_info.accessionNumber);
            WRITE_DCM_STRING(key, buf);
        }
        else {
            WRITE_DCM_STRING(key, 0);
        }

        // Modality
        // TODO: this is hardcoded!!
        key.set(0x0008, 0x0060);
        WRITE_DCM_STRING(key, "MR");

        // Manufacturer
        key.set(0x0008, 0x0070);
        if (sys_info.systemVendor) {
            WRITE_DCM_STRING(key, sys_info.systemVendor->c_str());
        }
        else {
            WRITE_DCM_STRING(key, "UNKNOWN");
        }

        // Institution Name
        key.set(0x0008, 0x0080);
        if (sys_info.institutionName) {
            WRITE_DCM_STRING(key, sys_info.institutionName->c_str());
        }
        else {
            WRITE_DCM_STRING(key, "UNKNOWN");
        }

        // Referring Physician's Name
        key.set(0x0008, 0x0090);
        if (study_info.referringPhysicianName) {
            WRITE_DCM_STRING(key, study_info.referringPhysicianName->c_str());
        }
        else {
            WRITE_DCM_STRING(key, "");
        }

        // Station Name
        key.set(0x0008, 0x1010);
        if (sys_info.stationName) {
            WRITE_DCM_STRING(key, sys_info.stationName->c_str());
        }
        else {
            WRITE_DCM_STRING(key, "");
        }

        // Study Description
        key.set(0x0008, 0x1030);
        if (study_info.studyDescription) {
            WRITE_DCM_STRING(key, study_info.studyDescription->c_str());
        }
        else {
            WRITE_DCM_STRING(key, "");
        }

        // Series Description
        key.set(0x0008, 0x103E);
        if (meas_info.seriesDescription) {
            WRITE_DCM_STRING(key, meas_info.seriesDescription->c_str());
        }
        else {
            WRITE_DCM_STRING(key, "");
        }

        // Manufacturer's Model Name
        key.set(0x0008, 0x1090);
        if (sys_info.systemModel) {
            WRITE_DCM_STRING(key, sys_info.systemModel->c_str());
        }
        else {
            WRITE_DCM_STRING(key, "");
        }

        // Referenced SOP Instance UIDs
        std::vector<ISMRMRD::ReferencedImageSequence> refs(meas_info.referencedImageSequence);
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
            GDEBUG("Failed to write 0x0010 Group Length\n");
            return GADGET_FAIL;
        }

        // Patient Name
        key.set(0x0010, 0x0010);
        if (patient_info.patientName) {
            WRITE_DCM_STRING(key, patient_info.patientName->c_str());
        }
        else {
            WRITE_DCM_STRING(key, "None");
        }

        // Patient ID
        key.set(0x0010, 0x0020);
        if (patient_info.patientID) {
            WRITE_DCM_STRING(key, patient_info.patientID->c_str());
        }
        else {
            WRITE_DCM_STRING(key, "0");
        }

        // Patient Birthdate
        key.set(0x0010, 0x0030);
        if (patient_info.patientBirthdate) {
            std::string d(patient_info.patientBirthdate.get());
            d.erase(std::remove(d.begin(), d.end(), '-'), d.end());
            WRITE_DCM_STRING(key, d.c_str());
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
                WRITE_DCM_STRING(key, patient_info.patientGender->c_str());
            }
        }
        else {
            WRITE_DCM_STRING(key, "");
        }

        // Patient Age
        key.set(0x0010, 0x1010);
        if (patient_info.patientBirthdate && meas_info.seriesDate) {
            boost::gregorian::date bday(boost::gregorian::from_simple_string(patient_info.patientBirthdate.get()));
            boost::gregorian::date seriesDate(boost::gregorian::from_simple_string(meas_info.seriesDate.get()));

            boost::gregorian::days age = seriesDate - bday;

            long age_in_years = age.days() / 365;

            ACE_OS::snprintf(buf, BUFSIZE, "%03ldY", age_in_years);
            WRITE_DCM_STRING(key, buf);
        }
        else {
            WRITE_DCM_STRING(key, "000Y");
        }

        // Patient Weight
        key.set(0x0010, 0x1030);
        if (patient_info.patientWeight_kg) {
            ACE_OS::snprintf(buf, BUFSIZE, "%f", *patient_info.patientWeight_kg);
            WRITE_DCM_STRING(key, buf);
        }
        else {
            WRITE_DCM_STRING(key, "0.0");
        }

        // Group Length
        key.set(0x0018, 0x0000);
        status = dataset->insertEmptyElement(key);
        if (!status.good()) {
            GDEBUG("Failed to write 0x0018 Group Length\n");
            return GADGET_FAIL;
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
        WRITE_DCM_STRING(key, scanningSequence.c_str());
        key.set(0x0018, 0x0021);
        WRITE_DCM_STRING(key, sequenceVariant.c_str());
        key.set(0x0018, 0x0022);
        WRITE_DCM_STRING(key, scanOptions.c_str());
        key.set(0x0018, 0x0023);
        WRITE_DCM_STRING(key, mrAcquisitionType.c_str());

        // Angio Flag
        // TODO: hardcoded
        key.set(0x0018, 0x0025);
        WRITE_DCM_STRING(key, "N");

        // Slice Thickness
        // This will need updated if the "reconSpace.fieldOfView_mm.z" field
        // is changed in the ISMRMRD populating code (client)
        key.set(0x0018, 0x0050);
        ACE_OS::snprintf(buf, BUFSIZE, "%f", r_space.fieldOfView_mm.z / std::max(r_space.matrixSize.z, (unsigned short)1));
        WRITE_DCM_STRING(key, buf);

        // Spacing Between Slices
        key.set(0x0018, 0x0088);
        ACE_OS::snprintf(buf, BUFSIZE, "%f", r_space.fieldOfView_mm.z);
        WRITE_DCM_STRING(key, buf);

        // Repetition Time
        if (seq_info.TR.is_present() && seq_info.TR.get().size() > 0)
        {
            key.set(0x0018, 0x0080);
            ACE_OS::snprintf(buf, BUFSIZE, "%f", seq_info.TR.get().front());
            WRITE_DCM_STRING(key, buf);
        }

        // Echo Time
        if (seq_info.TE.is_present() && seq_info.TE.get().size() > 0)
        {
            key.set(0x0018, 0x0081);
            ACE_OS::snprintf(buf, BUFSIZE, "%f", seq_info.TE.get().front());
            WRITE_DCM_STRING(key, buf);
        }

        // Inversion Time
        if (seq_info.TI.is_present() && seq_info.TI.get().size()>0)
        {
            key.set(0x0018, 0x0082);
            ACE_OS::snprintf(buf, BUFSIZE, "%f", seq_info.TI.get().front());
            WRITE_DCM_STRING(key, buf);
        }

        // Flip Angle
        if (seq_info.flipAngle_deg.is_present() && seq_info.flipAngle_deg.get().size()>0)
        {
            key.set(0x0018, 0x1314);
            ACE_OS::snprintf(buf, BUFSIZE, "%ld", (long)seq_info.flipAngle_deg.get().front());
            WRITE_DCM_STRING(key, buf);
        }

        // Imaging Frequency in tenths of MHz ???
        key.set(0x0018, 0x0084);
        ACE_OS::snprintf(buf, BUFSIZE, "%f", (float)exp_cond.H1resonanceFrequency_Hz / 10000000.);
        WRITE_DCM_STRING(key, buf);

        // Magnetic Field Strength (T)
        key.set(0x0018, 0x0087);
        if (sys_info.systemFieldStrength_T) {
            ACE_OS::snprintf(buf, BUFSIZE, "%f", *sys_info.systemFieldStrength_T);
            WRITE_DCM_STRING(key, buf);
        }
        else {
            WRITE_DCM_STRING(key, "3.0");
        }


        // Echo Train Length
        if (h.encoding[0].echoTrainLength) {
            key.set(0x0018, 0x0091);
            ACE_OS::snprintf(buf, BUFSIZE, "%ld", (long)*h.encoding[0].echoTrainLength);
            WRITE_DCM_STRING(key, buf);
        }
        else {
            WRITE_DCM_STRING(key, "1");
        }

        // Percent Sampling
        // TODO: hardcoded
        key.set(0x0018, 0x0093);
        WRITE_DCM_STRING(key, "100");

        // Percent Phase FOV
        // TODO: hardcoded
        key.set(0x0018, 0x0094);
        WRITE_DCM_STRING(key, "100");

        // Protocol Name
        if (meas_info.protocolName) {
            key.set(0x0018, 0x1030);
            WRITE_DCM_STRING(key, meas_info.protocolName.get().c_str());
        }
        else {
            WRITE_DCM_STRING(key, "");
        }

        // Trigger Time - TODO: use Image Meta Data
        key.set(0x0018, 0x1060);
        WRITE_DCM_STRING(key, "0.0");

        // Reconstruction Diameter (FOV) - TODO: ?
        key.set(0x0018, 0x1100);

        // Frequency Encoding Direction - TODO: use Image Meta Data
        key.set(0x0018, 0x1312);
        WRITE_DCM_STRING(key, "ROW");

        // Patient Position
        key.set(0x0018, 0x5100);
        WRITE_DCM_STRING(key, meas_info.patientPosition.c_str());

        /****************************************/
        // Group Length
        key.set(0x0020, 0x0000);
        status = dataset->insertEmptyElement(key);
        if (!status.good()) {
            GDEBUG("Failed to write 0x0020 Group Length\n");
            return GADGET_FAIL;
        }

        // Study Instance UID
        key.set(0x0020, 0x000D);
        if (study_info.studyInstanceUID) {
            WRITE_DCM_STRING(key, study_info.studyInstanceUID->c_str());
        }

        // Study ID
        if (study_info.studyID) {
            key.set(0x0020, 0x0010);
            WRITE_DCM_STRING(key, study_info.studyID->c_str());
        }
        else {
            WRITE_DCM_STRING(key, "0");
        }

        // Store Series Instance UID for later
        if (meas_info.seriesInstanceUIDRoot) {
            seriesIUIDRoot = *meas_info.seriesInstanceUIDRoot;
        }

        // Frame of Reference UID
        if (meas_info.frameOfReferenceUID) {
            key.set(0x0020, 0x0052);
            WRITE_DCM_STRING(key, meas_info.frameOfReferenceUID->c_str());
        }

        /****************************************/
        // Group Length
        key.set(0x0028, 0x0000);
        status = dataset->insertEmptyElement(key);
        if (!status.good()) {
            GDEBUG("Failed to write 0x0028 Group Length\n");
            return GADGET_FAIL;
        }

        // Samples Per Pixel
        key.set(0x0028, 0x0002);
        // TODO: hardcoded
        WRITE_DCM_STRING(key, "1");

        // Photometric Interpretation
        key.set(0x0028, 0x0004);
        // TODO: hardcoded
        WRITE_DCM_STRING(key, "MONOCHROME2");

        // Pixel Spacing (Array of len 2)
        key.set(0x0028, 0x0030);
        float pixel_spacing_X = r_space.fieldOfView_mm.x / r_space.matrixSize.x;
        float pixel_spacing_Y = r_space.fieldOfView_mm.y / r_space.matrixSize.y;
        ACE_OS::snprintf(buf, BUFSIZE, "%.3f\\%.3f", pixel_spacing_X, pixel_spacing_Y);
        WRITE_DCM_STRING(key, buf);

        // Bits Allocated
        key.set(0x0028, 0x0100);
        WRITE_DCM_STRING(key, "16");
        // Bits Stored
        key.set(0x0028, 0x0101);
        WRITE_DCM_STRING(key, "16");
        // High Bit
        key.set(0x0028, 0x0102);
        WRITE_DCM_STRING(key, "15");
        // Pixel Representation
        key.set(0x0028, 0x0103);
        WRITE_DCM_STRING(key, "1");

        //GDEBUG("Finished populating DICOM fields\n");

        /* clean up the buffer we created for ACE_OS::snprintf */
        delete[] buf;

        return GADGET_OK;
    }

    int DicomFinishGadget::process(GadgetContainerMessage<ISMRMRD::ImageHeader>* m1)
    {
        if (!this->controller_) {
            GERROR("Cannot return result to controller, no controller set");
            return GADGET_FAIL;
        }

        // --------------------------------------------------
        ISMRMRD::ImageHeader *img = m1->getObjectPtr();

        uint16_t data_type = img->data_type;

        if (data_type == ISMRMRD::ISMRMRD_USHORT)
        {
            GadgetContainerMessage< hoNDArray< unsigned short > >* datamb = AsContainerMessage< hoNDArray< unsigned short > >(m1->cont());
            if (!datamb)
            {
                GERROR("DicomFinishGadget::process, invalid image message objects\n");
                return GADGET_FAIL;
            }

            if (this->write_data_attrib(m1, datamb) != GADGET_OK)
            {
                GERROR("DicomFinishGadget::write_data_attrib failed for unsigned short ... \n");
                return GADGET_FAIL;
            }
        }
        else if (data_type == ISMRMRD::ISMRMRD_SHORT)
        {
            GadgetContainerMessage< hoNDArray< short > >* datamb = AsContainerMessage< hoNDArray< short > >(m1->cont());
            if (!datamb)
            {
                GERROR("DicomFinishGadget::process, invalid image message objects\n");
                return GADGET_FAIL;
            }

            if (this->write_data_attrib(m1, datamb) != GADGET_OK)
            {
                GERROR("DicomFinishGadget::write_data_attrib failed for short ... \n");
                return GADGET_FAIL;
            }
        }
        else if (data_type == ISMRMRD::ISMRMRD_UINT)
        {
            GadgetContainerMessage< hoNDArray< unsigned int > >* datamb = AsContainerMessage< hoNDArray< unsigned int > >(m1->cont());
            if (!datamb)
            {
                GERROR("DicomFinishGadget::process, invalid image message objects\n");
                return GADGET_FAIL;
            }

            if (this->write_data_attrib(m1, datamb) != GADGET_OK)
            {
                GERROR("DicomFinishGadget::write_data_attrib failed for unsigned int ... \n");
                return GADGET_FAIL;
            }
        }
        else if (data_type == ISMRMRD::ISMRMRD_INT)
        {
            GadgetContainerMessage< hoNDArray< int > >* datamb = AsContainerMessage< hoNDArray< int > >(m1->cont());
            if (!datamb)
            {
                GERROR("DicomFinishGadget::process, invalid image message objects\n");
                return GADGET_FAIL;
            }

            if (this->write_data_attrib(m1, datamb) != GADGET_OK)
            {
                GERROR("DicomFinishGadget::write_data_attrib failed for int ... \n");
                return GADGET_FAIL;
            }
        }
        else if (data_type == ISMRMRD::ISMRMRD_FLOAT)
        {
            GadgetContainerMessage< hoNDArray< float > >* datamb = AsContainerMessage< hoNDArray< float > >(m1->cont());
            if (!datamb)
            {
                GERROR("DicomFinishGadget::process, invalid image message objects\n");
                return GADGET_FAIL;
            }

            if (this->write_data_attrib(m1, datamb) != GADGET_OK)
            {
                GERROR("DicomFinishGadget::write_data_attrib failed for float ... \n");
                return GADGET_FAIL;
            }
        }
        else if (data_type == ISMRMRD::ISMRMRD_DOUBLE)
        {
            GadgetContainerMessage< hoNDArray< double > >* datamb = AsContainerMessage< hoNDArray< double > >(m1->cont());
            if (!datamb)
            {
                GERROR("DicomFinishGadget::process, invalid image message objects\n");
                return GADGET_FAIL;
            }

            if (this->write_data_attrib(m1, datamb) != GADGET_OK)
            {
                GERROR("DicomFinishGadget::write_data_attrib failed for double ... \n");
                return GADGET_FAIL;
            }
        }
        else if (data_type == ISMRMRD::ISMRMRD_CXFLOAT)
        {
            GERROR("DicomFinishGadget::process, does not supprot ISMRMRD_CXFLOAT data type\n");
            return GADGET_FAIL;
        }
        else if (data_type == ISMRMRD::ISMRMRD_CXDOUBLE)
        {
            GERROR("DicomFinishGadget::process, does not supprot ISMRMRD_CXDOUBLE data type\n");
            return GADGET_FAIL;
        }

        return GADGET_OK;
    }

    GADGET_FACTORY_DECLARE(DicomFinishGadget)

} /* namespace Gadgetron */
