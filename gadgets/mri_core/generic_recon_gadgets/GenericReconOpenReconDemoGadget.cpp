
#include "GenericReconOpenReconDemoGadget.h"
#include <cmath>
#include <sstream>

namespace Gadgetron {

    GenericReconOpenReconDemoGadget::GenericReconOpenReconDemoGadget() : triggered_in_close_(false), BaseClass()
    {
    }

    GenericReconOpenReconDemoGadget::~GenericReconOpenReconDemoGadget()
    {
    }

    int GenericReconOpenReconDemoGadget::process_config(ACE_Message_Block *mb)
    {
        GADGET_CHECK_RETURN(BaseClass::process_config(mb) == GADGET_OK, GADGET_FAIL);

        // -------------------------------------------------

        ISMRMRD::IsmrmrdHeader h;
        try
        {
            deserialize(mb->rd_ptr(), h);
        }
        catch (...)
        {
            GDEBUG("Error parsing ISMRMRD Header");
        }

        this->num_encoding_spaces_ = h.encoding.size();
        GDEBUG_STREAM("Number of encoding spaces: " << this->num_encoding_spaces_);

        return GADGET_OK;
    }

    int GenericReconOpenReconDemoGadget::process(Gadgetron::GadgetContainerMessage<IsmrmrdReconData> *m1)
    {
        process_called_times_++;

        IsmrmrdReconData *recon_bit_ = m1->getObjectPtr();
        if (recon_bit_->rbit_.size() > num_encoding_spaces_)
        {
            GWARN_STREAM("Incoming recon_bit has more encoding spaces than the protocol : " << recon_bit_->rbit_.size()
                                                                                            << " instead of "
                                                                                            << num_encoding_spaces_);
        }

        // save an acq header
        this->acq_header_ = recon_bit_->rbit_[0].data_.headers_[0];

        m1->release();

        return GADGET_OK;
    }


void GenericReconOpenReconDemoGadget::fill_image_header(IsmrmrdImageArray& res, bool is_rgb)
{
        size_t RO  = res.data_.get_size(0);
        size_t E1  = res.data_.get_size(1);
        size_t E2  = res.data_.get_size(2);
        size_t CHA = res.data_.get_size(3);
        size_t N   = res.data_.get_size(4);
        size_t S   = res.data_.get_size(5);
        size_t SLC = res.data_.get_size(6);

        res.headers_.create(N, S, SLC);
        if (res.meta_.size() < N*S*SLC )
            res.meta_.resize(N * S * SLC);

        size_t n, s, slc, im_ind(1);

        ISMRMRD::AcquisitionHeader& acq_header = this->acq_header_;

        for (slc = 0; slc < SLC; slc++)
        {
            for (s = 0; s < S; s++)
            {
                for (n = 0; n < N; n++)
                {
                    ISMRMRD::ImageHeader& im_header = res.headers_(n, s, slc);
                    ISMRMRD::MetaContainer& meta    = res.meta_[n + s * N + slc * N * S];

                    im_header.version         = acq_header.version;
                    im_header.data_type       = ISMRMRD::ISMRMRD_CXFLOAT;
                    im_header.measurement_uid = acq_header.measurement_uid;

                    im_header.matrix_size[0] = (uint16_t)RO; // pixel size will be 1.0
                    im_header.matrix_size[1] = (uint16_t)E1;
                    im_header.matrix_size[2] = (uint16_t)E2;

                    im_header.field_of_view[0] = RO; // pixel size will be 1.0
                    im_header.field_of_view[1] = E1;
                    im_header.field_of_view[2] = E2;

                    im_header.channels = (uint16_t)CHA;
                    if (is_rgb)
                    {
                        GADGET_CHECK_THROW(CHA==3);
                    }

                    im_header.position[0] = acq_header.position[0];
                    im_header.position[1] = acq_header.position[1];
                    im_header.position[2] = acq_header.position[2];

                    im_header.read_dir[0] = acq_header.read_dir[0];
                    im_header.read_dir[1] = acq_header.read_dir[1];
                    im_header.read_dir[2] = acq_header.read_dir[2];

                    im_header.phase_dir[0] = acq_header.phase_dir[0];
                    im_header.phase_dir[1] = acq_header.phase_dir[1];
                    im_header.phase_dir[2] = acq_header.phase_dir[2];

                    im_header.slice_dir[0] = acq_header.slice_dir[0];
                    im_header.slice_dir[1] = acq_header.slice_dir[1];
                    im_header.slice_dir[2] = acq_header.slice_dir[2];

                    im_header.patient_table_position[0] = acq_header.patient_table_position[0];
                    im_header.patient_table_position[1] = acq_header.patient_table_position[1];
                    im_header.patient_table_position[2] = acq_header.patient_table_position[2];

                    im_header.average    = acq_header.idx.average;
                    im_header.slice      = slc;
                    im_header.contrast   = acq_header.idx.contrast;
                    im_header.phase      = acq_header.idx.phase;
                    im_header.repetition = s;
                    im_header.set        = n;

                    im_header.acquisition_time_stamp = acq_header.acquisition_time_stamp;

                    im_header.physiology_time_stamp[0] = acq_header.physiology_time_stamp[0];
                    im_header.physiology_time_stamp[1] = acq_header.physiology_time_stamp[1];
                    im_header.physiology_time_stamp[2] = acq_header.physiology_time_stamp[2];

                    if (is_rgb)
                    {
                        im_header.image_type         = 6; // ISMRMRD_IMTYPE_RGB
                    }
                    else
                    {
                        im_header.image_type         = ISMRMRD::ISMRMRD_IMTYPE_MAGNITUDE;
                    }
                    im_header.image_index        = (uint16_t)(im_ind);
                    im_ind++;
                    im_header.image_series_index = 0;

                    memcpy(im_header.user_int, acq_header.user_int, sizeof(int32_t) * ISMRMRD::ISMRMRD_USER_INTS);
                    memcpy(im_header.user_float, acq_header.user_float, sizeof(float) * ISMRMRD::ISMRMRD_USER_FLOATS);

                    im_header.attribute_string_len = 0;

                    meta.set("PatientPosition", (double)im_header.position[0]);
                    meta.append("PatientPosition", (double)im_header.position[1]);
                    meta.append("PatientPosition", (double)im_header.position[2]);

                    meta.set("read_dir", (double)im_header.read_dir[0]);
                    meta.append("read_dir", (double)im_header.read_dir[1]);
                    meta.append("read_dir", (double)im_header.read_dir[2]);

                    meta.set("phase_dir", (double)im_header.phase_dir[0]);
                    meta.append("phase_dir", (double)im_header.phase_dir[1]);
                    meta.append("phase_dir", (double)im_header.phase_dir[2]);

                    meta.set("slice_dir", (double)im_header.slice_dir[0]);
                    meta.append("slice_dir", (double)im_header.slice_dir[1]);
                    meta.append("slice_dir", (double)im_header.slice_dir[2]);

                    meta.set("patient_table_position", (double)im_header.patient_table_position[0]);
                    meta.append("patient_table_position", (double)im_header.patient_table_position[1]);
                    meta.append("patient_table_position", (double)im_header.patient_table_position[2]);

                    meta.set("acquisition_time_stamp", (long)im_header.acquisition_time_stamp);

                    meta.set("physiology_time_stamp", (long)im_header.physiology_time_stamp[0]);
                    meta.append("physiology_time_stamp", (long)im_header.physiology_time_stamp[1]);
                    meta.append("physiology_time_stamp", (long)im_header.physiology_time_stamp[2]);

                    size_t ui;
                    for (ui = 0; ui < ISMRMRD::ISMRMRD_USER_INTS; ui++)
                    {
                        std::ostringstream str;
                        str << "user_int_" << ui;
                        meta.append(str.str().c_str(), (long)res.headers_(n, s, slc).user_int[ui]);
                    }

                    for (ui = 0; ui < ISMRMRD::ISMRMRD_USER_FLOATS; ui++)
                    {
                        std::ostringstream str;
                        str << "user_float_" << ui;
                        meta.append(str.str().c_str(), (long)res.headers_(n, s, slc).user_float[ui]);
                    }

                    meta.set("measurementID", this->measurement_id_.c_str());
                    meta.set("protocolName", this->protocol_name_.c_str());
                    meta.set("patientID", this->patient_.c_str());
                    meta.set("studyID", this->study_.c_str());
                    meta.set("measurementNumber", this->measurement_.c_str());
                    meta.set("deviceID", this->device_.c_str());
                    meta.set("patient_position", this->patient_position_.c_str());
                }
            }
        }
    }

    int GenericReconOpenReconDemoGadget::close(unsigned long flags)
    {
        GDEBUG_STREAM("GenericReconOpenReconDemoGadget - close(flags) : " << flags);
        if (BaseClass::close(flags) != GADGET_OK) return GADGET_FAIL;

        if (!this->triggered_in_close_)
        {
            this->triggered_in_close_ = true;

            GDEBUG_STREAM("GenericReconOpenReconDemoGadget - assemble and send out images ... ");

            // RGB images
            IsmrmrdImageArray rgb_3D;
            this->create_rgb_3D_image(rgb_3D);
            this->fill_image_header(rgb_3D, true);

            if (!debug_folder_full_path_.empty())
            {
                hoNDArray<std::complex<float>> res;
                res = rgb_3D.data_;
                res.squeeze();
                gt_exporter_.export_array_complex(res, debug_folder_full_path_ + "rgb_3D");
            }

            IsmrmrdImageArray rgb_2D;
            this->create_rgb_2D_image(rgb_2D);
            this->fill_image_header(rgb_2D, true);

            if (!debug_folder_full_path_.empty())
            {
                hoNDArray<std::complex<float>> res;
                res = rgb_2D.data_;
                res.squeeze();
                gt_exporter_.export_array_complex(res, debug_folder_full_path_ + "rgb_2D");
            }

            // Images with contours and landmarks
            IsmrmrdImageArray im_contours_landmarks;
            this->create_image_with_contours(im_contours_landmarks);
            this->fill_image_header(im_contours_landmarks, false);
            if (!debug_folder_full_path_.empty())
            {
                im_contours_landmarks.data_.print(std::cout);
                gt_exporter_.export_array_complex(im_contours_landmarks.data_, debug_folder_full_path_ + "im_contours_landmarks");
            }

            // Report page
            IsmrmrdImageArray reports;
            this->create_reports(reports);
            this->fill_image_header(reports, false);
            if (!debug_folder_full_path_.empty())
            {
                gt_exporter_.export_array_complex(reports.data_, debug_folder_full_path_ + "reports");
            }

            this->send_out_image_array(rgb_3D, 0, 100, GADGETRON_IMAGE_REGULAR);
            this->send_out_image_array(rgb_2D, 0, 101, GADGETRON_IMAGE_REGULAR);
            this->send_out_image_array(im_contours_landmarks, 0, 102, GADGETRON_IMAGE_REGULAR);
            this->send_out_image_array(reports, 0, 103, GADGETRON_IMAGE_RECON_FIGURE);
        }

        return GADGET_OK;
    }

    void GenericReconOpenReconDemoGadget::create_rgb_3D_image(IsmrmrdImageArray& im)
    {
        size_t RO = 256;
        size_t E1 = 144;
        size_t E2 = 36;
        size_t CHA = 3; // for rgb
        size_t N = 1;
        size_t S = 1;
        size_t LOC = 1;

        im.data_.create(RO, E1, E2, CHA, N, S, LOC);
        Gadgetron::clear(im.data_);
        im.headers_.create(N, S, LOC);
        im.meta_.resize(N*S*LOC);

        size_t ro, e1, e2;
        for (e2=0; e2<E2; e2++)
        {
            bool first_half_e2 = (e2<E2/2) ? true : false;
            for (e1=0; e1<E1; e1++)
            {
                bool first_half_e1 = (e1<E1/2) ? true : false;
                for (ro=0; ro<RO; ro++)
                {
                    bool first_half_ro = (ro<RO/2) ? true : false;

                    if (first_half_e2 && first_half_e1 && first_half_ro) 
                    {
                        im.data_(ro, e1, e2, 0, 0, 0, 0) = 255.0; // red
                    }
                    else if (first_half_e2 && first_half_e1 && !first_half_ro)
                    {
                        im.data_(ro, e1, e2, 1, 0, 0, 0) = 255.0; // green
                    }
                    else if (first_half_e2 && !first_half_e1 && first_half_ro)
                    {
                        im.data_(ro, e1, e2, 2, 0, 0, 0) = 255.0; // blue
                    }
                    else if (first_half_e2 && !first_half_e1 && !first_half_ro)
                    {
                        im.data_(ro, e1, e2, 1, 0, 0, 0) = 255.0;
                        im.data_(ro, e1, e2, 2, 0, 0, 0) = 255.0;
                    }
                    else if (!first_half_e2 && first_half_e1 && first_half_ro) 
                    {
                        im.data_(ro, e1, e2, 0, 0, 0, 0) = 156.0;
                        im.data_(ro, e1, e2, 1, 0, 0, 0) = 46.0;
                        im.data_(ro, e1, e2, 2, 0, 0, 0) = 122.0;
                    }
                    else if (!first_half_e2 && first_half_e1 && !first_half_ro)
                    {
                        im.data_(ro, e1, e2, 0, 0, 0, 0) = 150.0;
                        im.data_(ro, e1, e2, 1, 0, 0, 0) = 156.0;
                        im.data_(ro, e1, e2, 2, 0, 0, 0) = 50.0;
                    }
                    else if (!first_half_e2 && !first_half_e1 && first_half_ro)
                    {
                        im.data_(ro, e1, e2, 0, 0, 0, 0) = 209.0;
                        im.data_(ro, e1, e2, 1, 0, 0, 0) = 139.0;
                        im.data_(ro, e1, e2, 2, 0, 0, 0) = 108.0;
                    }
                    else
                    {
                        im.data_(ro, e1, e2, 0, 0, 0, 0) = 128.0;
                        im.data_(ro, e1, e2, 1, 0, 0, 0) = 128.0;
                        im.data_(ro, e1, e2, 2, 0, 0, 0) = 128.0;
                    }
                }
            }
        }
    }

    void GenericReconOpenReconDemoGadget::create_rgb_2D_image(IsmrmrdImageArray& im)
    {
        size_t RO = 192;
        size_t E1 = 320;
        size_t E2 = 1;
        size_t CHA = 3; // for rgb
        size_t N = 12;
        size_t S = 1;
        size_t LOC = 1;

        im.data_.create(RO, E1, E2, CHA, N, S, LOC);
        Gadgetron::clear(im.data_);
        im.headers_.create(N, S, LOC);
        im.meta_.resize(N*S*LOC);

        size_t ro, e1, n;
        for (n=0; n<N; n++)
        {
            for (e1=0; e1<E1; e1++)
            {
                bool first_half_e1 = (e1<E1/2) ? true : false;
                for (ro=0; ro<RO; ro++)
                {
                    bool first_half_ro = (ro<RO/2) ? true : false;

                    if (first_half_e1 && first_half_ro) 
                    {
                        im.data_(ro, e1, 0, 0, n, 0, 0) = ((n+1)*255.0)/N;
                    }
                    else if (first_half_e1 && !first_half_ro)
                    {
                        im.data_(ro, e1, 0, 1, n, 0, 0) = ((n+1)*255.0)/N;
                    }
                    else if (!first_half_e1 && first_half_ro)
                    {
                        im.data_(ro, e1, 0, 2, n, 0, 0) = ((n+1)*255.0)/N;
                    }
                    else
                    {
                        im.data_(ro, e1, 0, 1, n, 0, 0) = ((n+1)*255.0)/N;
                        im.data_(ro, e1, 0, 2, n, 0, 0) = ((n+1)*255.0)/N;
                    }
                }
            }
        }
    }

    void set_ismrmrd_meta_values(ISMRMRD::MetaContainer& attrib, const std::string& name, const std::vector<double>& v)
    {
        if (v.empty())
        {
            GWARN_STREAM("set_ismrmrd_meta_values, input vector is empty ... " << name);
            return;
        }

        attrib.set(name.c_str(), v[0]);

        size_t ii;
        for (ii = 1; ii < v.size(); ii++)
        {
            attrib.append(name.c_str(), v[ii]);
        }
    }

    void GenericReconOpenReconDemoGadget::create_image_with_contours(IsmrmrdImageArray& im)
    {
        size_t RO = 160;
        size_t E1 = 120;
        size_t E2 = 1;
        size_t CHA = 1;
        size_t N = 8;
        size_t S = 1;
        size_t LOC = 1;

        im.data_.create(RO, E1, E2, CHA, N, S, LOC);
        Gadgetron::clear(im.data_);
        im.headers_.create(N, S, LOC);
        im.meta_.resize(N*S*LOC);

        size_t n, ro, e1;
        for (n=0; n<N; n++)
        {
            for (e1=0; e1<E1; e1++)
            {
                for (ro=0; ro<RO; ro++)
                {
                    im.data_(ro, e1, 0, 0, n, 0, 0) = 1.0 + (n*1024.0)/(N-1);
                }
            }

            // generate and attach contours
            std::vector<double> contour;

            // color of this contour, rgb
            contour.push_back(1.0);
            contour.push_back(0.0);
            contour.push_back(0.0);
            // line width
            contour.push_back(3.0);

            // let's draw circle for N=100 points
            double radius = 0.9 * std::min(RO/2, E1/2);
            double curr_radius = (0.1 + n/N) * radius;
            for (size_t p=0; p<360; p+=5)
            {
                double p_ro = RO/2 + curr_radius * std::cos(p * M_PI / 180.0); 
                double p_e1 = E1/2 + curr_radius * std::sin(p * M_PI / 180.0); 
                contour.push_back(p_ro);
                contour.push_back(p_e1);
            }

            std::stringstream ostr_c;
            ostr_c << "GT_ROI_" << n;
            set_ismrmrd_meta_values(im.meta_[n], ostr_c.str(), contour);

            // generate and attach landmarks
            std::vector<double> landmark;
            landmark.push_back((n*1.0)/(N-1)); // red
            landmark.push_back(1.0); // green
            landmark.push_back((n*0.8)/(N-1)); // blue
            landmark.push_back(2*n+1); // size
            landmark.push_back(n % 4); // type
            landmark.push_back(RO/(n+2)); // location in pixel, along the 1st dimension
            landmark.push_back(E1/(n+2)); // location in pixel, along the 2nd dimension
            landmark.push_back(0); // location in pixel, along the 3rd dimension

            std::stringstream ostr_l;
            ostr_l << "GT_Landmark_" << n;
            set_ismrmrd_meta_values(im.meta_[n], ostr_l.str(), landmark);

            // consider to bypass further process, e.g. distortion correction
            im.meta_[n].set("GADGETRON_DIRECT_IMAGE_SEND", 1.0);
        }
    }

    void GenericReconOpenReconDemoGadget::create_reports(IsmrmrdImageArray& reports)
    {
        size_t RO = 512;
        size_t E1 = 512;
        size_t E2 = 1;
        size_t CHA = 1;
        size_t N = 1;
        size_t S = 1;
        size_t LOC = 1;

        reports.data_.create(RO, E1, E2, CHA, N, S, LOC);
        Gadgetron::clear(reports.data_);
        reports.headers_.create(N, S, LOC);
        reports.meta_.resize(N*S*LOC);

        // plot a colored square with color LUT
        size_t ro, e1;
        for (e1=128; e1<384; e1++)
        {
            for (ro=128; ro<384; ro++)
            {
                reports.data_(ro, e1, 0, 0, 0, 0, 0) = e1 - 128;
            }
        }

        reports.meta_[0].set(GADGETRON_IMAGE_WINDOWCENTER, (long)(100));
        reports.meta_[0].set(GADGETRON_IMAGE_WINDOWWIDTH, (long)(200));

        reports.meta_[0].set("Correct_image_orientation", 1.0);
        reports.meta_[0].set("GADGETRON_DIRECT_IMAGE_SEND", 1.0);

        reports.meta_[0].set(GADGETRON_IMAGE_COLORMAP, "MicroDeltaHotMetal.pal");
    }

    GADGET_FACTORY_DECLARE(GenericReconOpenReconDemoGadget)
}