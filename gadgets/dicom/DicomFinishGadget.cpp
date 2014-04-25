
#include <vector>

#include "GadgetIsmrmrdReadWrite.h"
#include "DicomFinishGadget.h"

using namespace std;

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

    // Parse ISMRMRD XML header
    boost::shared_ptr<ISMRMRD::ismrmrdHeader> cfg = parseIsmrmrdXMLHeader(string(mb->rd_ptr()));

    //GADGET_DEBUG1("Processing XML config in DicomFinishAttribGadget\n");

    // Ensure DICOM dictionary is loaded
    if (!dcmDataDict.isDictionaryLoaded()) {
        GADGET_DEBUG1("Dictionary not loaded!  Set DCMDICTPATH\n");
        return GADGET_FAIL;
    }

    ISMRMRD::experimentalConditionsType exp_cond = cfg->experimentalConditions();

    if (!cfg->subjectInformation().present()) {
        GADGET_DEBUG1("Header missing SubjectInformation parameters\n");
        return GADGET_FAIL;
    }
    ISMRMRD::subjectInformationType patient_info = cfg->subjectInformation().get();

    if (!cfg->studyInformation().present()) {
        GADGET_DEBUG1("Header missing StudyInformation parameters\n");
        return GADGET_FAIL;
    }
    ISMRMRD::studyInformationType study_info = cfg->studyInformation().get();

    if (!cfg->measurementInformation().present()) {
        GADGET_DEBUG1("Header missing MeasurementInformation parameters\n");
        return GADGET_FAIL;
    }
    ISMRMRD::measurementInformationType meas_info = cfg->measurementInformation().get();

    if (!cfg->acquisitionSystemInformation().present()) {
        GADGET_DEBUG1("Header missing AcquisitionSystemInformation parameters\n");
        return GADGET_FAIL;
    }
    ISMRMRD::acquisitionSystemInformationType sys_info = cfg->acquisitionSystemInformation().get();

    if (!cfg->sequenceParameters().present()) {
        GADGET_DEBUG1("Header missing SequenceTiming parameters\n");
        return GADGET_FAIL;
    }
    ISMRMRD::sequenceParametersType seq_info = cfg->sequenceParameters().get();

    // Ensure that the XML header contains the DICOM parameters
    if (!cfg->dicomParameters().present()) {
        GADGET_DEBUG1("Header missing DICOM parameters\n");
        return GADGET_OK;
    }

    ISMRMRD::ismrmrdHeader::encoding_sequence e_seq = cfg->encoding();
    ISMRMRD::encodingSpaceType e_space = (*e_seq.begin()).encodedSpace();
    ISMRMRD::encodingSpaceType r_space = (*e_seq.begin()).reconSpace();
    ISMRMRD::encodingLimitsType e_limits = (*e_seq.begin()).encodingLimits();

    ISMRMRD::dicomParametersType dcm_params = cfg->dicomParameters().get();
    ISMRMRD::MRImageModule mr_image(dcm_params.MRImageModule().get());

    DcmDataset *dataset = dcmFile.getDataset();
    DcmMetaInfo *metainfo = dcmFile.getMetaInfo();


    // Store initial Series Number for later
    if (meas_info.initialSeriesNumber().present()) {
        this->initialSeriesNumber = meas_info.initialSeriesNumber().get();
    } else {
        this->initialSeriesNumber = 0;
    }


    // Set the Application Entity Title in the DICOM Meta Info section
    // The rest of the Meta Info will be automatically populated by DCMTK
    if (sys_info.stationName().present()) {
        status = metainfo->putAndInsertString(DcmTagKey(0x0002,0x0016),
                sys_info.stationName().get().c_str());
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
    key.set(0x0008, 0x0008);
    if (mr_image.imageType().present()) {
        WRITE_DCM_STRING(key, mr_image.imageType().get().c_str());
    } else {
        WRITE_DCM_STRING(key, "ORIGINAL\\PRIMARY\\OTHER");
    }

    // SOPClassUID
    key.set(0x0008, 0x0016);
    WRITE_DCM_STRING(key, UID_MRImageStorage);

    // Study Date
    key.set(0x0008, 0x0020);
    ACE_OS::snprintf(buf, BUFSIZE, "%04d%02d%02d", study_info.studyDate().get().year(), study_info.studyDate().get().month(), study_info.studyDate().get().day());
    WRITE_DCM_STRING(key, buf);

    // Series Date
    key.set(0x0008, 0x0021);
    ACE_OS::snprintf(buf, BUFSIZE, "%04d%02d%02d", meas_info.seriesDate().get().year(), meas_info.seriesDate().get().month(), meas_info.seriesDate().get().day());
    WRITE_DCM_STRING(key, buf);
    // Acquisition Date
    key.set(0x0008, 0x0022);
    WRITE_DCM_STRING(key, buf);
    // Content Date
    key.set(0x0008, 0x0023);
    WRITE_DCM_STRING(key, buf);

    // Study Time
    key.set(0x0008, 0x0030);
    ACE_OS::snprintf(buf, BUFSIZE, "%02d%02d%02d", study_info.studyTime().get().hours(), study_info.studyTime().get().minutes(), (int)study_info.studyTime().get().seconds());
    WRITE_DCM_STRING(key, buf);

    // Series Time
    key.set(0x0008, 0x0031);
    ACE_OS::snprintf(buf, BUFSIZE, "%02d%02d%02d", meas_info.seriesTime().get().hours(), meas_info.seriesTime().get().minutes(), (int)meas_info.seriesTime().get().seconds());
    WRITE_DCM_STRING(key, buf);

    // Acquisition Time
    key.set(0x0008, 0x0032);
    WRITE_DCM_STRING(key, buf);

    // Content Time
    key.set(0x0008, 0x0033);
    WRITE_DCM_STRING(key, buf);

    // Accession Number
    key.set(0x0008, 0x0050);
    if (study_info.accessionNumber().present()) {
        ACE_OS::snprintf(buf, BUFSIZE, "%d", (int)study_info.accessionNumber().get());
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
    if (sys_info.systemVendor().present()) {
        WRITE_DCM_STRING(key, sys_info.systemVendor().get().c_str());
    } else {
        WRITE_DCM_STRING(key, "UNKNOWN");
    }

    // Institution Name
    key.set(0x0008, 0x0080);
    if (sys_info.institutionName().present()) {
        WRITE_DCM_STRING(key, sys_info.institutionName().get().c_str());
    } else {
        WRITE_DCM_STRING(key, "UNKNOWN");
    }

    // Referring Physician's Name
    key.set(0x0008, 0x0090);
    if (study_info.referringPhysicianName().present()) {
        WRITE_DCM_STRING(key, study_info.referringPhysicianName().get().c_str());
    } else {
        WRITE_DCM_STRING(key, "");
    }

    // Station Name
    key.set(0x0008, 0x1010);
    if (sys_info.stationName().present()) {
        WRITE_DCM_STRING(key, sys_info.stationName().get().c_str());
    } else {
        WRITE_DCM_STRING(key, "");
    }

    // Study Description
    key.set(0x0008, 0x1030);
    if (study_info.studyDescription().present()) {
        WRITE_DCM_STRING(key, study_info.studyDescription().get().c_str());
    } else {
        WRITE_DCM_STRING(key, "");
    }

    // Series Description
    key.set(0x0008, 0x103E);
    if (meas_info.seriesDescription().present()) {
        WRITE_DCM_STRING(key, meas_info.seriesDescription().get().c_str());
    } else {
        WRITE_DCM_STRING(key, "");
    }

    // Manufacturer's Model Name
    key.set(0x0008, 0x1090);
    if (sys_info.systemModel().present()) {
        WRITE_DCM_STRING(key, sys_info.systemModel().get().c_str());
    } else {
        WRITE_DCM_STRING(key, "");
    }

    // Referenced SOP Instance UIDs
    if (dcm_params.referencedImageSequence().present()) {
        ISMRMRD::referencedImageSequence refs = dcm_params.referencedImageSequence().get();
        DcmItem *ref_sequence;
        string ref_uid;
        for (unsigned int i = 0; i < refs.referencedSOPInstanceUID().size(); i++) {
            ref_uid = refs.referencedSOPInstanceUID()[i];

            if (ref_uid.length() > 0) {   // Only write non-empty strings
                if (dataset->findOrCreateSequenceItem(key, ref_sequence, -2 /* append */).good()) {
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
    if (patient_info.patientName().present()) {
        WRITE_DCM_STRING(key, patient_info.patientName().get().c_str());
    } else {
        WRITE_DCM_STRING(key, "None");
    }

    // Patient ID
    key.set(0x0010, 0x0020);
    if (patient_info.patientID().present()) {
        WRITE_DCM_STRING(key, patient_info.patientID().get().c_str());
    } else {
        WRITE_DCM_STRING(key, "0");
    }

    // Patient Birthdate
    key.set(0x0010, 0x0030);
    if (patient_info.patientBirthdate().present()) {
        ACE_OS::snprintf(buf, BUFSIZE, "%04d%02d%02d", patient_info.patientBirthdate().get().year(),
                patient_info.patientBirthdate().get().month(), patient_info.patientBirthdate().get().day());
        WRITE_DCM_STRING(key, buf);
    } else {
        status = dataset->insertEmptyElement(key);
    }

    // Patient Sex
    key.set(0x0010, 0x0040);
    if (patient_info.patientGender().present()) {
        if (patient_info.patientGender().get() == "O") {
            status = dataset->insertEmptyElement(key);
        }
        else {
            WRITE_DCM_STRING(key, patient_info.patientGender().get().c_str());
        }
    } else {
        WRITE_DCM_STRING(key, "");
    }

    // Patient Age
    key.set(0x0010, 0x1010);
    if (patient_info.patientBirthdate().present()) {
        ACE_OS::snprintf(buf, BUFSIZE, "%03uY", meas_info.seriesDate().get().year() - patient_info.patientBirthdate().get().year());
        WRITE_DCM_STRING(key, buf);
    } else {
        WRITE_DCM_STRING(key, "000Y");
    }

    // Patient Weight
    key.set(0x0010, 0x1030);
    if (patient_info.patientWeight_kg().present()) {
        ACE_OS::snprintf(buf, BUFSIZE, "%f", patient_info.patientWeight_kg().get());
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

    // Scanning Sequence
    if (mr_image.scanningSequence().present()) {
        key.set(0x0018, 0x0020);
        WRITE_DCM_STRING(key, mr_image.scanningSequence().get().c_str());
    } else {
        WRITE_DCM_STRING(key, "RM");
    }

    // Sequence Variant
    if (mr_image.sequenceVariant().present()) {
        key.set(0x0018, 0x0021);
        WRITE_DCM_STRING(key, mr_image.sequenceVariant().get().c_str());
    } else {
        WRITE_DCM_STRING(key, "NONE");
    }

    // Scan Options
    if (mr_image.scanOptions().present()) {
        key.set(0x0018, 0x0022);
        WRITE_DCM_STRING(key, mr_image.scanOptions().get().c_str());
    } else {
        WRITE_DCM_STRING(key, "NONE");
    }

    // Acquisition Type
    if (mr_image.mrAcquisitionType().present()) {
        key.set(0x0018, 0x0023);
        WRITE_DCM_STRING(key, mr_image.mrAcquisitionType().get().c_str());
    } else {
        WRITE_DCM_STRING(key, "2D");
    }

    // Angio Flag
    // TODO: hardcoded
    key.set(0x0018, 0x0025);
    WRITE_DCM_STRING(key, "N");

    // Slice Thickness
    // This will need updated if the "reconSpace.fieldOfView_mm.z" field
    // is changed in the ISMRMRD populating code (client)
    key.set(0x0018, 0x0050);
    ACE_OS::snprintf(buf, BUFSIZE, "%f", cfg->encoding().front().reconSpace().fieldOfView_mm().z());
    WRITE_DCM_STRING(key, buf);

    // Repetition Time
    key.set(0x0018, 0x0080);
    ACE_OS::snprintf(buf, BUFSIZE, "%f", seq_info.TR().front());
    WRITE_DCM_STRING(key, buf);

    // Echo Time
    key.set(0x0018, 0x0081);
    ACE_OS::snprintf(buf, BUFSIZE, "%f", seq_info.TE().front());
    WRITE_DCM_STRING(key, buf);

    // Inversion Time
    key.set(0x0018, 0x0082);
    ACE_OS::snprintf(buf, BUFSIZE, "%f", seq_info.TI().front());
    WRITE_DCM_STRING(key, buf);

    // Imaging Frequency in tenths of MHz ???
    key.set(0x0018, 0x0084);
    ACE_OS::snprintf(buf, BUFSIZE, "%f", (float)exp_cond.H1resonanceFrequency_Hz() / 10000000.);
    WRITE_DCM_STRING(key, buf);

    // Magnetic Field Strength (T)
    key.set(0x0018, 0x0087);
    if (sys_info.systemFieldStrength_T().present()) {
        ACE_OS::snprintf(buf, BUFSIZE, "%f", sys_info.systemFieldStrength_T().get());
        WRITE_DCM_STRING(key, buf);
    } else {
        WRITE_DCM_STRING(key, "3.0");
    }

    // Spacing Between Slices
    key.set(0x0018, 0x0088);
    ACE_OS::snprintf(buf, BUFSIZE, "%f", cfg->encoding().front().reconSpace().fieldOfView_mm().z());
    WRITE_DCM_STRING(key, buf);

    // Echo Train Length
    if (mr_image.echoTrainLength().present()) {
        key.set(0x0018, 0x0091);
        ACE_OS::snprintf(buf, BUFSIZE, "%ld", (long)mr_image.echoTrainLength().get());
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
    if (meas_info.protocolName().present()) {
        key.set(0x0018, 0x1030);
        WRITE_DCM_STRING(key, meas_info.protocolName().get().c_str());
    } else {
        WRITE_DCM_STRING(key, "");
    }

    // Trigger Time
    if (mr_image.triggerTime().present()) {
        key.set(0x0018, 0x1060);
        ACE_OS::snprintf(buf, BUFSIZE, "%f", mr_image.triggerTime().get());
        WRITE_DCM_STRING(key, buf);
    } else {
        WRITE_DCM_STRING(key, "0.0");
    }

    // Reconstruction Diameter (FOV)
    // TODO: hmm
    key.set(0x0018, 0x1100);

    // Frequency Encoding Direction
    if (mr_image.freqEncodingDirection().present()) {
        key.set(0x0018, 0x1312);
        WRITE_DCM_STRING(key, mr_image.freqEncodingDirection().get().c_str());
    } else {
        WRITE_DCM_STRING(key, "COL");
    }

    // Flip Angle
    if (mr_image.flipAngle_deg().present()) {
        key.set(0x0018, 0x1314);
        ACE_OS::snprintf(buf, BUFSIZE, "%d", (int)mr_image.flipAngle_deg().get());
        WRITE_DCM_STRING(key, buf);
    } else {
        WRITE_DCM_STRING(key, "0");
    }

    // Patient Position
    key.set(0x0018, 0x5100);
    WRITE_DCM_STRING(key, meas_info.patientPosition().c_str());

    // Group Length
    key.set(0x0020, 0x0000);
    status = dataset->insertEmptyElement(key);
    if (!status.good()) {
        GADGET_DEBUG1("Failed to write 0x0020 Group Length\n");
        return GADGET_FAIL;
    }

    // Study Instance UID
    key.set(0x0020, 0x000D);
    WRITE_DCM_STRING(key, dcm_params.studyInstanceUID().c_str());

    // Study ID
    if (study_info.studyID().present()) {
        key.set(0x0020, 0x0010);
        WRITE_DCM_STRING(key, study_info.studyID().get().c_str());
    } else {
        WRITE_DCM_STRING(key, "0");
    }

    // Store Series Instance UID for later
    if (dcm_params.seriesInstanceUIDRoot().present()) {
        seriesIUIDRoot = dcm_params.seriesInstanceUIDRoot().get();
    }

    // Frame of Reference UID
    if (dcm_params.frameOfReferenceUID().present()) {
        key.set(0x0020, 0x0052);
        WRITE_DCM_STRING(key, dcm_params.frameOfReferenceUID().get().c_str());
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
    float pixel_spacing_X = r_space.fieldOfView_mm().x() / r_space.matrixSize().x();
    float pixel_spacing_Y = r_space.fieldOfView_mm().y() / r_space.matrixSize().y();
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
    /* update the image_data_type.
     * There is currently no SIGNED SHORT type so this will have to suffice */
    m1->getObjectPtr()->image_data_type = ISMRMRD::DATA_UNSIGNED_SHORT;

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
    std::map<unsigned int, string>::iterator it = seriesIUIDs.find(series_number);

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
        seriesIUIDs[series_number] = string(newuid);
    }
    WRITE_DCM_STRING(key, seriesIUIDs[series_number].c_str());

    // At a minimum, to put the DICOM image back into the database,
    // you must change the SOPInstanceUID.
    key.set(0x0008, 0x0018);        // SOPInstanceUID
    const char *root;
    if (seriesIUIDRoot.length() > 0) {
        root = string(seriesIUIDRoot, 0, 20).c_str();
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
