
#include <vector>
#include "boost/date_time/gregorian/gregorian.hpp"            

#include "DicomFinishGadget.h"
#include "ismrmrd_xml.h"

// Used for windowing using short ints
#define PIX_RANGE_MAX    (+32767)
#define PIX_RANGE_MIN    (-32768)


// Writes a DICOM string value at the given location in the header
// Saves keystrokes
#define WRITE_DCM_STRING(k, s)    \
    do {                                                                    \
        status = dataset->putAndInsertString(k, s);            \
        if (!status.good()) {                                               \
            GADGET_DEBUG2("Failed to insert DICOM field (0x%04X,0x%04X) at "\
                "line %u\n", k.getGroup(), k.getElement(), __LINE__);       \
            return GADGET_FAIL;                                             \
        }                                                                   \
    } while (0)

namespace Gadgetron {

template <typename T>
int DicomFinishGadget<T>::process_config(ACE_Message_Block* mb)
{
    OFCondition status;
    DcmTagKey key;
    long BUFSIZE = 1024;
    char *buf = new char[BUFSIZE];  // used for writing numbers as strings in DCMTK

    ISMRMRD::IsmrmrdHeader h;
    deserialize(mb->rd_ptr(), h);

    // Ensure DICOM dictionary is loaded
    if (!dcmDataDict.isDictionaryLoaded()) {
        GADGET_DEBUG1("Dictionary not loaded!  Set DCMDICTPATH\n");
        return GADGET_FAIL;
    }

    ISMRMRD::ExperimentalConditions exp_cond = h.experimentalConditions;

    if (!h.subjectInformation) {
        GADGET_DEBUG1("Header missing SubjectInformation parameters\n");
        return GADGET_FAIL;
    }

    ISMRMRD::SubjectInformation patient_info = *h.subjectInformation;

    if (!h.studyInformation) {
      GADGET_DEBUG1("Header missing StudyInformation parameters\n");
      return GADGET_FAIL;
    }

    ISMRMRD::StudyInformation study_info = *h.studyInformation;

    if (!h.measurementInformation) {
        GADGET_DEBUG1("Header missing MeasurementInformation parameters\n");
        return GADGET_FAIL;
    }

    ISMRMRD::MeasurementInformation meas_info = *h.measurementInformation;

    if (!h.acquisitionSystemInformation) {
        GADGET_DEBUG1("Header missing AcquisitionSystemInformation parameters\n");
        return GADGET_FAIL;
    }

    ISMRMRD::AcquisitionSystemInformation sys_info = *h.acquisitionSystemInformation;

    if (!h.sequenceParameters) {
        GADGET_DEBUG1("Header missing SequenceTiming parameters\n");
        return GADGET_FAIL;
    }

    ISMRMRD::SequenceParameters seq_info = *h.sequenceParameters;

    if (h.encoding.size() == 0) {
      GADGET_DEBUG2("Number of encoding spaces: %d\n", h.encoding.size());
      GADGET_DEBUG1("This Gadget needs an encoding description\n");
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
    } else {
      this->initialSeriesNumber = 0;
    }


    // Set the Application Entity Title in the DICOM Meta Info section
    // The rest of the Meta Info will be automatically populated by DCMTK
    if (sys_info.stationName) {
        status = metainfo->putAndInsertString(DcmTagKey(0x0002,0x0016),
                sys_info.stationName->c_str());
        if (!status.good()) {
            GADGET_DEBUG1("Failed to set AET in MetaInfo\n");
            return GADGET_FAIL;
        }
    } else {
        status = metainfo->putAndInsertString(DcmTagKey(0x0002,0x0016), "none");
        if (!status.good()) {
            GADGET_DEBUG1("Failed to set AET in MetaInfo\n");
            return GADGET_FAIL;
        }
    }

    // Group Length
    key.set(0x0008, 0x0000);
    status = dataset->insertEmptyElement(key);
    if (status.bad()) {
        GADGET_DEBUG1("Failed to write 0x0008 Group Length\n");
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
    } else {
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
    } else {
        WRITE_DCM_STRING(key, "UNKNOWN");
    }

    // Institution Name
    key.set(0x0008, 0x0080);
    if (sys_info.institutionName) {
        WRITE_DCM_STRING(key, sys_info.institutionName->c_str());
    } else {
        WRITE_DCM_STRING(key, "UNKNOWN");
    }

    // Referring Physician's Name
    key.set(0x0008, 0x0090);
    if (study_info.referringPhysicianName) {
        WRITE_DCM_STRING(key, study_info.referringPhysicianName->c_str());
    } else {
        WRITE_DCM_STRING(key, "");
    }

    // Station Name
    key.set(0x0008, 0x1010);
    if (sys_info.stationName) {
        WRITE_DCM_STRING(key, sys_info.stationName->c_str());
    } else {
        WRITE_DCM_STRING(key, "");
    }

    // Study Description
    key.set(0x0008, 0x1030);
    if (study_info.studyDescription) {
      WRITE_DCM_STRING(key, study_info.studyDescription->c_str());
    } else {
      WRITE_DCM_STRING(key, "");
    }

    // Series Description
    key.set(0x0008, 0x103E);
    if (meas_info.seriesDescription) {
        WRITE_DCM_STRING(key, meas_info.seriesDescription->c_str());
    } else {
        WRITE_DCM_STRING(key, "");
    }

    // Manufacturer's Model Name
    key.set(0x0008, 0x1090);
    if (sys_info.systemModel) {
        WRITE_DCM_STRING(key, sys_info.systemModel->c_str());
    } else {
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
        GADGET_DEBUG1("Failed to write 0x0010 Group Length\n");
        return GADGET_FAIL;
    }

    // Patient Name
    key.set(0x0010, 0x0010);
    if (patient_info.patientName) {
        WRITE_DCM_STRING(key, patient_info.patientName->c_str());
    } else {
        WRITE_DCM_STRING(key, "None");
    }

    // Patient ID
    key.set(0x0010, 0x0020);
    if (patient_info.patientID) {
        WRITE_DCM_STRING(key, patient_info.patientID->c_str());
    } else {
        WRITE_DCM_STRING(key, "0");
    }

    // Patient Birthdate
    key.set(0x0010, 0x0030);
    if (patient_info.patientBirthdate) {
        //ACE_OS::snprintf(buf, BUFSIZE, "%04d%02d%02d", patient_info.patientBirthdate().get().year(),
        //        patient_info.patientBirthdate().get().month(), patient_info.patientBirthdate().get().day());
        //WRITE_DCM_STRING(key, buf);
    } else {
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
    } else {
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
    } else {
        WRITE_DCM_STRING(key, "000Y");
    }

    // Patient Weight
    key.set(0x0010, 0x1030);
    if (patient_info.patientWeight_kg) {
        ACE_OS::snprintf(buf, BUFSIZE, "%f", *patient_info.patientWeight_kg);
        WRITE_DCM_STRING(key, buf);
    } else {
        WRITE_DCM_STRING(key, "0.0");
    }

    // Group Length
    key.set(0x0018, 0x0000);
    status = dataset->insertEmptyElement(key);
    if (!status.good()) {
        GADGET_DEBUG1("Failed to write 0x0018 Group Length\n");
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
            } else if (it->name == "sequenceVariant") {
                sequenceVariant = it->value;
            } else if (it->name == "scanOptions") {
                scanOptions = it->value;
            } else if (it->name == "mrAcquisitionType") {
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
    ACE_OS::snprintf(buf, BUFSIZE, "%f", r_space.fieldOfView_mm.z);
    WRITE_DCM_STRING(key, buf);

    // Repetition Time
    key.set(0x0018, 0x0080);
    ACE_OS::snprintf(buf, BUFSIZE, "%f", seq_info.TR.front());
    WRITE_DCM_STRING(key, buf);

    // Echo Time
    key.set(0x0018, 0x0081);
    ACE_OS::snprintf(buf, BUFSIZE, "%f", seq_info.TE.front());
    WRITE_DCM_STRING(key, buf);

    // Inversion Time
    key.set(0x0018, 0x0082);
    ACE_OS::snprintf(buf, BUFSIZE, "%f", seq_info.TI.front());
    WRITE_DCM_STRING(key, buf);

    // Flip Angle
    key.set(0x0018, 0x1314);
    ACE_OS::snprintf(buf, BUFSIZE, "%ld", (long)seq_info.flipAngle_deg.front());
    WRITE_DCM_STRING(key, buf);

    // Imaging Frequency in tenths of MHz ???
    key.set(0x0018, 0x0084);
    ACE_OS::snprintf(buf, BUFSIZE, "%f", (float)exp_cond.H1resonanceFrequency_Hz / 10000000.);
    WRITE_DCM_STRING(key, buf);

    // Magnetic Field Strength (T)
    key.set(0x0018, 0x0087);
    if (sys_info.systemFieldStrength_T) {
        ACE_OS::snprintf(buf, BUFSIZE, "%f", *sys_info.systemFieldStrength_T);
        WRITE_DCM_STRING(key, buf);
    } else {
        WRITE_DCM_STRING(key, "3.0");
    }

    // Spacing Between Slices
    key.set(0x0018, 0x0088);
    ACE_OS::snprintf(buf, BUFSIZE, "%f", r_space.fieldOfView_mm.z);
    WRITE_DCM_STRING(key, buf);

    // Echo Train Length
    if (h.encoding[0].echoTrainLength) {
        key.set(0x0018, 0x0091);
        ACE_OS::snprintf(buf, BUFSIZE, "%ld", (long)*h.encoding[0].echoTrainLength);
        WRITE_DCM_STRING(key, buf);
    } else {
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
    } else {
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
        GADGET_DEBUG1("Failed to write 0x0020 Group Length\n");
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
    } else {
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
        GADGET_DEBUG1("Failed to write 0x0028 Group Length\n");
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

    //GADGET_DEBUG1("Finished populating DICOM fields\n");

    /* clean up the buffer we created for ACE_OS::snprintf */
    delete[] buf;

    return GADGET_OK;
}

template <typename T>
int DicomFinishGadget<T>::process(GadgetContainerMessage<ISMRMRD::ImageHeader>* m1, GadgetContainerMessage< hoNDArray< T > >* m2)
{
    if (!this->controller_)
    {
        ACE_DEBUG( (LM_DEBUG,
                    ACE_TEXT("Cannot return result to controller, no controller set")) );
        return -1;
    }

    GadgetContainerMessage<hoNDArray< ACE_INT16 > > *pixels =
            new GadgetContainerMessage<hoNDArray< ACE_INT16 > >();
    boost::shared_ptr< std::vector<size_t> > dims = m2->getObjectPtr()->get_dimensions();

    try {
        pixels->getObjectPtr()->create(dims.get());
    } catch (bad_alloc& err) {
        GADGET_DEBUG1("Unable to create short storage in DicomFinishGadget");
        return GADGET_FAIL;
    }

    /* create ImageHeader and hoNDArray pointers for better readability */
    ISMRMRD::ImageHeader *img = m1->getObjectPtr();
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
    if (pixels->getObjectPtr()->get_number_of_elements() > 0) {
        min_pix_val = src[0];
        max_pix_val = src[0];
    }
    for (unsigned long i = 0; i < pixels->getObjectPtr()->get_number_of_elements(); i++) {
        T pix_val = src[i];
        // search for minimum and maximum pixel values
        if (pix_val < min_pix_val) min_pix_val = pix_val;
        if (pix_val > max_pix_val) max_pix_val = pix_val;
        sum_pix_val += pix_val / 4; // scale by 25% to avoid overflow

        // copy/cast the pixel value to a short int
        dst[i] = static_cast<ACE_INT16>(pix_val);
    }
    T mean_pix_val = (sum_pix_val * 4) / pixels->getObjectPtr()->get_number_of_elements();

    /* replace the old 'message2' with the new data */
    m1->cont(pixels);
    /* release the old data array */
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
    ACE_OS::snprintf(buf, BUFSIZE, "%d", img->contrast);
    WRITE_DCM_STRING(key, buf);

    // Acquisition Matrix ... Image Dimensions
    // Defined as: [frequency rows, frequency columns, phase rows, phase columns]
    key.set(0x0018, 0x1310);
    ACE_UINT16 im_dim[4] = {0,0,0,0};
    // frequency encoding direction is always the COLUMN
    im_dim[0] = img->matrix_size[0];
    im_dim[3] = img->matrix_size[1];
    status = dataset->putAndInsertUint16Array(key, im_dim, 4);
    if (!status.good()) {
        GADGET_DEBUG1("Failed to stuff image dimensions\n");
        return GADGET_FAIL;
    }

    // Series Number
    // Only write a number if the image_series_index is positive and non-zero
    key.set(0x0020, 0x0011);
    ACE_OS::snprintf(buf, BUFSIZE, "%ld", this->initialSeriesNumber * 100 + img->image_series_index);
    WRITE_DCM_STRING(key, buf);

    // Image Number
    key.set(0x0020, 0x0013);
    ACE_OS::snprintf(buf, BUFSIZE, "%d", img->image_index + 1);
    WRITE_DCM_STRING(key, buf);

    // Image Position (Patient)
    float corner[3];

    corner[0] = img->position[0] -
            (img->field_of_view[0] / 2.0) * img->read_dir[0] -
            (img->field_of_view[1] / 2.0) * img->phase_dir[0];
    corner[1] = img->position[1] -
            (img->field_of_view[0] / 2.0) * img->read_dir[1] -
            (img->field_of_view[1] / 2.0) * img->phase_dir[1];
    corner[2] = img->position[2] -
            (img->field_of_view[0] / 2.0) * img->read_dir[2] -
            (img->field_of_view[1] / 2.0) * img->phase_dir[2];

    key.set(0x0020, 0x0032);
    ACE_OS::snprintf(buf, BUFSIZE, "%.4f\\%.4f\\%.4f", corner[0], corner[1], corner[2]);
    WRITE_DCM_STRING(key, buf);

    // Image Orientation
    // read_dir, phase_dir, and slice_dir were calculated in
    // a DICOM/patient coordinate system, so just plug them in
    key.set(0x0020, 0x0037);
    ACE_OS::snprintf(buf, BUFSIZE, "%.4f\\%.4f\\%.4f\\%.4f\\%.4f\\%.4f",
            img->read_dir[0], img->read_dir[1], img->read_dir[2],
            img->phase_dir[0], img->phase_dir[1], img->phase_dir[2]);
    WRITE_DCM_STRING(key, buf);

    // Slice Location
    key.set(0x0020, 0x1041);
    ACE_OS::snprintf(buf, BUFSIZE, "%f", img->position[2]);
    WRITE_DCM_STRING(key, buf);

    // Columns
    key.set(0x0028, 0x0010);
    ACE_OS::snprintf(buf, BUFSIZE, "%d", img->matrix_size[0]);
    WRITE_DCM_STRING(key, buf);

    // Rows
    key.set(0x0028, 0x0011);
    ACE_OS::snprintf(buf, BUFSIZE, "%d", img->matrix_size[1]);
    WRITE_DCM_STRING(key, buf);

    // Simple windowing using pixel values calculated earlier...
    int mid_pix_val = (max_pix_val + min_pix_val) / 2;
    int window_center = (mid_pix_val + mean_pix_val) / 2;
    int window_width_left = window_center - min_pix_val;
    int window_width_right = max_pix_val - window_center;
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
        GADGET_DEBUG1("Failed to write 0x7fe0 Group Length\n");
        return GADGET_FAIL;
    }

    // Pixel Data
    if ((unsigned long)img->matrix_size[0] * (unsigned long)img->matrix_size[1] !=
                data->get_number_of_elements()) {
        GADGET_DEBUG1("Mismatch in image dimensions and available data\n");
        return GADGET_FAIL;
    }
    key.set(0x7fe0, 0x0010);
    status = dataset->putAndInsertUint16Array(key, (unsigned short *)data->get_data_ptr(),
            data->get_number_of_elements());
    if (!status.good()) {
        GADGET_DEBUG1("Failed to stuff Pixel Data\n");
        return GADGET_FAIL;
    }

    // Series Instance UID = generated here
    key.set(0x0020, 0x000E);
    unsigned short series_number = img->image_series_index + 1;

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
        } else {
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
    } else {
       root = "1.2.840.113619.2.156";
    }
    char newuid[65];
    dcmGenerateUniqueIdentifier(newuid, root);
    WRITE_DCM_STRING(key, newuid);

    /* clean up the char[] we created for ACE_OS::snprintf */
    delete[] buf;

    GadgetContainerMessage<DcmFileFormat>* mdcm = new GadgetContainerMessage<DcmFileFormat>();

    *mdcm->getObjectPtr() = dcmFile;

    GadgetContainerMessage<GadgetMessageIdentifier>* mb =
        new GadgetContainerMessage<GadgetMessageIdentifier>();

    mb->getObjectPtr()->id = GADGET_MESSAGE_DICOM;

    mb->cont(mdcm);

    int ret =  this->controller_->output_ready(mb);

    //GADGET_DEBUG1("Finished Finishing DICOM\n");

    if ( (ret < 0) ) {
        GADGET_DEBUG1("Failed to return message to controller\n");
        return GADGET_FAIL;
    }

    return GADGET_OK;
}

//Declare factories for the various template instances
GADGET_FACTORY_DECLARE(DicomFinishGadgetFLOAT)
GADGET_FACTORY_DECLARE(DicomFinishGadgetUSHORT)

} /* namespace Gadgetron */
