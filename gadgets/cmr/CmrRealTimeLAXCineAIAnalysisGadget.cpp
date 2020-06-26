
#include "CmrRealTimeLAXCineAIAnalysisGadget.h"
#include <cstdio>
#include "python_toolbox.h"
#include "hoNDImage_util.h"
#include "GadgetronTimer.h"
#include "mri_core_def.h"
#include "ismrmrd/meta.h"
#include "hoNDBSpline.h"
#include "ismrmrd/xml.h"
#include "cmr_time_stamp.h"
#include "hoNDKLT.h"
#include "cmr_file_and_directory_handling.h"
#include "cmr_image_container_util.h"
#include "cmr_ismrmrd_util.h"

namespace Gadgetron { 

    CmrRealTimeLAXCineAIAnalysisGadget::CmrRealTimeLAXCineAIAnalysisGadget() : BaseClass()
    {
        try
        {
            Gadgetron::initialize_python();

            // boost::filesystem::path gadgetron_python_path
            //     = this->context.paths.gadgetron_home / "share" / "gadgetron" / "python";

            // Gadgetron::add_python_path(gadgetron_python_path.generic_string());

            // this->gt_home_ = gadgetron_python_path.generic_string();
            // GDEBUG_STREAM("Set up python path : " << this->gt_home_);

            char* gt_home = std::getenv("GADGETRON_HOME");
            std::string path_name;
            if (gt_home != NULL)
            {
                size_t pos = std::string(gt_home).rfind("gadgetron");
                gt_home[pos - 1] = '\0';
                path_name = std::string(gt_home) + std::string("/share/gadgetron/python");

                if (Gadgetron::add_python_path(path_name) == GADGET_FAIL)
                {
                    GERROR_STREAM("Failed to add path: " << path_name << " to python ... ");
                }

                this->gt_home_ = path_name;
                GDEBUG_STREAM("Set up python path : " << this->gt_home_);
            }
        }
        catch (...)
        {
            this->gt_home_.clear();
            GERROR_STREAM("Exception happened when adding  path to python ... ");
        }
    }

    CmrRealTimeLAXCineAIAnalysisGadget::~CmrRealTimeLAXCineAIAnalysisGadget()
    {
    }

    int CmrRealTimeLAXCineAIAnalysisGadget::process_config(ACE_Message_Block* mb)
    {
        GADGET_CHECK_RETURN(BaseClass::process_config(mb) == GADGET_OK, GADGET_FAIL);

        ISMRMRD::IsmrmrdHeader h;
        ISMRMRD::deserialize(mb->rd_ptr(), h);

        if (h.encoding.size() == 0)
        {
            GDEBUG("Missing encoding section");
            return GADGET_FAIL;
        }

        ISMRMRD::EncodingSpace e_space = h.encoding[0].encodedSpace;
        ISMRMRD::EncodingSpace r_space = h.encoding[0].reconSpace;
        ISMRMRD::EncodingLimits e_limits = h.encoding[0].encodingLimits;

        meas_max_idx_.kspace_encode_step_1 = (uint16_t)(h.encoding[0].reconSpace.matrixSize.y);

        meas_max_idx_.set = (e_limits.set && (e_limits.set->maximum > 0)) ? e_limits.set->maximum : 0;
        meas_max_idx_.phase = (e_limits.phase && (e_limits.phase->maximum > 0)) ? e_limits.phase->maximum : 0;

        meas_max_idx_.kspace_encode_step_2 = (uint16_t)(h.encoding[0].reconSpace.matrixSize.z);

        meas_max_idx_.contrast = (e_limits.contrast && (e_limits.contrast->maximum > 0)) ? e_limits.contrast->maximum : 0;

        meas_max_idx_.slice = (e_limits.slice && (e_limits.slice->maximum > 0)) ? e_limits.slice->maximum : 0;

        meas_max_idx_.repetition = e_limits.repetition ? e_limits.repetition->maximum : 0;

        meas_max_idx_.average = e_limits.average ? e_limits.average->maximum : 0;

        meas_max_idx_.segment = 0;

        GDEBUG_STREAM("meas_max_idx_.slice is " << meas_max_idx_.slice);

        return GADGET_OK;
    }

    int CmrRealTimeLAXCineAIAnalysisGadget::process(Gadgetron::GadgetContainerMessage< IsmrmrdImageArray >* m1)
    {
        if (perform_timing.value()) { gt_timer_local_.start("CmrRealTimeLAXCineAIAnalysisGadget::process"); }

        GDEBUG_CONDITION_STREAM(verbose.value(), "CmrRealTimeLAXCineAIAnalysisGadget::process(...) starts ... ");

        // -------------------------------------------------------------

        process_called_times_++;

        // -------------------------------------------------------------

        IsmrmrdImageArray* data = m1->getObjectPtr();

        // print out data info
        if (verbose.value())
        {
            GDEBUG_STREAM("----> CmrRealTimeLAXCineAIAnalysisGadget::process(...) has been called " << process_called_times_ << " times ...");
            std::stringstream os;
            data->data_.print(os);
            GDEBUG_STREAM(os.str());
        }

        // -------------------------------------------------------------

        size_t encoding = (size_t)data->meta_[0].as_long("encoding", 0);
        GADGET_CHECK_RETURN(encoding < num_encoding_spaces_, GADGET_FAIL);

        size_t RO = data->data_.get_size(0);
        size_t E1 = data->data_.get_size(1);
        size_t E2 = data->data_.get_size(2);
        size_t CHA = data->data_.get_size(3);
        size_t N = data->data_.get_size(4);
        size_t S = data->data_.get_size(5);
        size_t PHS = data->data_.get_size(6);

        std::stringstream os;
        os << "_encoding_" << encoding << "_processing_" << process_called_times_;
        std::string str = os.str();

        // -------------------------------------------------------------

        if (!this->debug_folder_full_path_.empty())
        {
            gt_exporter_.export_array_complex(data->data_, this->debug_folder_full_path_ + "data" + str);
        }

        // -------------------------------------------------------------

        // copy the incoming image array
        IsmrmrdImageArray lax = *data;

        // send out the incoming images first
        if (this->next()->putq(m1) == -1)
        {
            GERROR("CmrRealTimeLAXCineAIAnalysisGadget::process, passing incoming image array on to next gadget");
            m1->release();
            return GADGET_FAIL;
        }

        // call the AI analysis
        IsmrmrdImageArray lax_ai, report;
        int status = this->perform_LAX_detection_AI(lax, lax_ai);

        // send out AI analysis results
        if(status==GADGET_OK)
        {
            Gadgetron::GadgetContainerMessage< IsmrmrdImageArray >* cm1 = new Gadgetron::GadgetContainerMessage< IsmrmrdImageArray >();
            *cm1->getObjectPtr() = lax_ai;
            if (this->next()->putq(cm1) == -1)
            {
                GERROR("CmrRealTimeLAXCineAIAnalysisGadget::process, passing lax_ai image array on to next gadget");
                cm1->release();
            }
        }

        // -------------------------------------------------------------

        if (perform_timing.value()) { gt_timer_local_.stop(); }

        return GADGET_OK;
    }

    int CmrRealTimeLAXCineAIAnalysisGadget::close(unsigned long flags)
    {
        GDEBUG_CONDITION_STREAM(true, "CmrRealTimeLAXCineAIAnalysisGadget - close(flags) : " << flags);
        if (BaseClass::close(flags) != GADGET_OK) return GADGET_FAIL;
        return GADGET_OK;
    }

    void CmrRealTimeLAXCineAIAnalysisGadget::convert_array_to_image_container(IsmrmrdImageArray& lax, hoNDImageContainer2D < hoMRImage<float, 2> >& lax_container)
    {
        size_t RO = lax.data_.get_size(0);
        size_t E1 = lax.data_.get_size(1);
        size_t PHS = lax.meta_.size();

        std::vector<size_t> dim2D(2);
        dim2D[0] = RO;
        dim2D[1] = E1;

        std::vector<size_t> cols(1, PHS);
        lax_container.create(cols);

        for (size_t phs = 0; phs < PHS; phs++)
        {
            lax_container(0, phs).create(dim2D);
            lax_container(0, phs).header_ = lax.headers_(phs);
            lax_container(0, phs).attrib_ = lax.meta_[phs];

            std::complex<float>* pData = lax.data_.begin() + phs * RO * E1;

            for (size_t k = 0; k < RO * E1; k++)
            {
                lax_container(0, phs)(k) = std::abs(pData[k]);
            }
        }
    }

    void CmrRealTimeLAXCineAIAnalysisGadget::convert_image_container_to_array(hoNDImageContainer2D < hoMRImage<float, 2> >& lax_container, IsmrmrdImageArray& lax)
    {
        size_t RO = lax_container(0, 0).get_size(0);
        size_t E1 = lax_container(0, 0).get_size(1);
        std::vector<size_t> cols = lax_container.cols();

        size_t PHS = cols[0];

        lax.data_.create(RO, E1, 1, 1, PHS, 1, 1);
        lax.headers_.create(PHS, 1, 1);
        lax.meta_.resize(PHS);

        for (size_t phs=0; phs<PHS; phs++)
        {
            float* pData = NULL;

            pData = lax_container(0, phs).begin();
            memcpy(&lax.headers_(phs, 0, 0), &lax_container(0, phs).header_, sizeof(ISMRMRD::ImageHeader));
            lax.meta_[phs] = lax_container(0, phs).attrib_;

            std::complex<float>* pLaxData = &lax.data_(0, 0, 0, 0, phs, 0, 0);
            for (size_t k = 0; k < RO * E1; k++)
            {
                pLaxData[k] = pData[k];
            }
        }
    }

    void plot_cross_on_image(hoMRImage<float, 2>& im, float px, float py)
    {
        float max_v;
        size_t max_v_ind;
        Gadgetron::maxAbsolute(im, max_v, max_v_ind);

        size_t RO = im.get_size(0);
        size_t E1 = im.get_size(1);

        int s_ro = (int(py)) - 2;
        if (s_ro < 0) s_ro = 0;

        int e_ro = s_ro + 9;
        if (e_ro >= RO) e_ro = RO - 1;

        int s_e1 = (int(px)) - 2;
        if (s_e1 < 0) s_e1 = 0;

        int e_e1 = s_e1 + 9;
        if (e_e1 >= E1) e_e1 = E1 - 1;

        size_t k;
        for (k=s_ro; k<=e_ro; k++)
        {
            im(k, (s_e1 + e_e1) / 2) = max_v + 1;
        }

        for (k = s_e1; k <= e_e1; k++)
        {
            im((s_ro+e_ro)/2, k) = max_v + 1;
        }
    }

    void CmrRealTimeLAXCineAIAnalysisGadget::plot_landmarks_on_images(hoNDImageContainer2D < hoMRImage<float, 2> > & lax_container, const hoNDArray<float>& pts)
    {
        size_t RO = lax_container(0, 0).get_size(0);
        size_t E1 = lax_container(0, 0).get_size(1);
        std::vector<size_t> cols = lax_container.cols();

        size_t PHS = cols[0];

        GADGET_CHECK_THROW(pts.get_size(2) <= PHS);

        for (size_t phs = 0; phs < pts.get_size(2); phs++)
        {
            if(pts(0, 0, phs) >= 0)
                plot_cross_on_image(lax_container(0, phs), pts(0, 0, phs), pts(0, 1, phs));

            if (pts(1, 0, phs) >= 0)
                plot_cross_on_image(lax_container(0, phs), pts(1, 0, phs), pts(1, 1, phs));

            if (pts(2, 0, phs) >= 0)
                plot_cross_on_image(lax_container(0, phs), pts(2, 0, phs), pts(2, 1, phs));
        }
    }

    void CmrRealTimeLAXCineAIAnalysisGadget::attach_info_to_report(hoNDImageContainer2D < hoMRImage<float, 2> > & lax_container, const hoNDArray<float>& pts)
    {
        size_t RO = lax_container(0, 0).get_size(0);
        size_t E1 = lax_container(0, 0).get_size(1);
        std::vector<size_t> cols = lax_container.cols();

        size_t PHS = cols[0];

        GADGET_CHECK_THROW(pts.get_size(2)<=PHS);

        for (size_t phs = 0; phs < pts.get_size(2); phs++)
        {
            std::vector<double> pt(2);
            pt[0] = pts(0, 0, phs);
            pt[1] = pts(0, 1, phs);
            Gadgetron::set_ismrmrd_meta_values(lax_container(0, phs).attrib_, "Gadgetron_Anterior_PT", pt);

            pt[0] = pts(1, 0, phs);
            pt[1] = pts(1, 1, phs);
            Gadgetron::set_ismrmrd_meta_values(lax_container(0, phs).attrib_, "Gadgetron_Posterior_PT", pt);

            pt[0] = pts(2, 0, phs);
            pt[1] = pts(2, 1, phs);
            Gadgetron::set_ismrmrd_meta_values(lax_container(0, phs).attrib_, "Gadgetron_Apical_PT", pt);
        }
    }

    int CmrRealTimeLAXCineAIAnalysisGadget::perform_LAX_detection_AI(IsmrmrdImageArray& lax, IsmrmrdImageArray& lax_ai)
    {
        try
        {
            Gadgetron::GadgetronTimer gt_timer(false);

            if (!this->debug_folder_full_path_.empty())
            {
                gt_exporter_.export_array_complex(lax.data_, this->debug_folder_full_path_ + "CmrRealTimeLAXCineAIAnalysisGadget_LAX");
            }

            size_t RO = lax.data_.get_size(0);
            size_t E1 = lax.data_.get_size(1);
            size_t PHS = lax.meta_.size();

            ImageContainerMagType lax_container;

            std::vector<size_t> dim2D(2);
            dim2D[0] = RO;
            dim2D[1] = E1;

            this->convert_array_to_image_container(lax, lax_container);

            if (!this->debug_folder_full_path_.empty())
            {
                hoNDArray<float> output;
                Gadgetron::convert_container_to_4D_array(lax_container, output);
                output.squeeze();

                std::stringstream str;
                str << "lax_container";
                gt_exporter_.export_array(output, this->debug_folder_full_path_ + "/" + str.str());
            }

            // -------------------------------------------
            // resample to 1mm2
            ImageContainerMagType lax_highres;
            lax_highres.copyFrom(lax_container);

            double FOV_RO = lax_container(0, 0).header_.field_of_view[0];
            double FOV_E1 = lax_container(0, 0).header_.field_of_view[1];

            double pixel_size = this->pixel_size_send.value();

            size_t new_RO = (size_t)(FOV_RO / pixel_size + 0.5);
            size_t new_E1 = (size_t)(FOV_E1 / pixel_size + 0.5);

            std::vector<size_t> dim2D_out(2);
            dim2D_out[0] = new_RO;
            dim2D_out[1] = new_E1;

            GDEBUG_STREAM("RT CINE LAX, new image size " << " [ " << new_RO << " " << new_E1 << "]");

            Gadgetron::resample_image_container(lax_container, dim2D_out, lax_highres);

            if (!this->debug_folder_full_path_.empty())
            {
                hoNDArray<float> output;
                Gadgetron::convert_container_to_4D_array(lax_highres, output);
                output.squeeze();

                std::stringstream str;
                str << "lax_highres";
                gt_exporter_.export_array(output, this->debug_folder_full_path_ + "/" + str.str());
            }

            // -------------------------------------------
            // convert to dicom orientation

            ImageContainerMagType lax_highres_dicom;
            lax_highres_dicom.copyFrom(lax_highres);

            if (!this->debug_folder_full_path_.empty())
            {
                hoNDArray<float> output;
                Gadgetron::convert_container_to_4D_array(lax_highres_dicom, output);
                output.squeeze();

                std::stringstream str;
                str << "lax_highres_dicom";
                gt_exporter_.export_array(output, this->debug_folder_full_path_ + "/" + str.str());
            }

            RO = lax_highres_dicom(0, 0).get_size(0);
            E1 = lax_highres_dicom(0, 0).get_size(1);

            // -------------------------------------------

            size_t PHS_detected = PHS;

            hoNDArray<float> lax_images;
            lax_images.create(RO, E1, PHS_detected);
            for (size_t phs = 0; phs < PHS_detected; phs++)
            {
                memcpy(&lax_images(0, 0, phs), lax_highres_dicom(0, phs).begin(), lax_highres_dicom(0, phs).get_number_of_bytes());
            }

            if (!this->debug_folder_full_path_.empty())
            {
                std::stringstream str;
                str << "lax_images_for_AI";
                gt_exporter_.export_array(lax_images, this->debug_folder_full_path_ + "/" + str.str());
            }

            lax_images.print(std::cout);

            // ---------------------------------------------------------
            // call cmr landmark detection
            hoNDArray<float> pts, probs;

            if(!this->gt_home_.empty())
            {
                GDEBUG_STREAM("=============================================");

                // load model
                PythonFunction<boost::python::object> load_model_cmr_landmark_detection("gadgetron_cmr_landmark_detection", "load_model_cmr_landmark_detection");
                boost::python::object model = load_model_cmr_landmark_detection(this->gt_home_, this->lax_landmark_detection_model.value());
                bp::incref(model.ptr());

                // apply model
                {
                    PythonFunction< hoNDArray<float>, hoNDArray<float> > perform_cmr_landmark_detection("gadgetron_cmr_landmark_detection", "perform_cmr_landmark_detection");
                    std::tie(pts, probs) = perform_cmr_landmark_detection(lax_images, model, 1.0, 8, 0.1, this->oper_RO.value(), this->oper_E1.value());
                }

                pts.print(std::cout);  
                probs.print(std::cout);

                GDEBUG_STREAM("=============================================");

                if (!this->debug_folder_full_path_.empty())
                {
                    std::stringstream str;
                    str << "lax_images_landmark_pts";
                    gt_exporter_.export_array(pts, this->debug_folder_full_path_ + "/" + str.str());
                }

                if (!this->debug_folder_full_path_.empty())
                {
                    std::stringstream str;
                    str << "lax_images_landmark_probs";
                    gt_exporter_.export_array(probs, this->debug_folder_full_path_ + "/" + str.str());
                }

                // -------------------------------------------

                GDEBUG_STREAM("attach landmarks to images ...");
                this->plot_landmarks_on_images(lax_highres_dicom, pts);

                if (!this->debug_folder_full_path_.empty())
                {
                    std::vector<std::vector<std::string> > attribs;
                    Gadgetron::serialize_contrainer_meta_attributes(lax_highres_dicom, attribs);
                    Gadgetron::write_contrainer_meta_attributes(attribs, this->debug_folder_full_path_ + "/lax_highres_dicom_attribs");
                }

                if (!this->debug_folder_full_path_.empty())
                {
                    hoNDArray<float> output;
                    Gadgetron::convert_container_to_4D_array(lax_highres_dicom, output);
                    output.squeeze();

                    std::stringstream str;
                    str << "lax_highres_dicom_with_pts";
                    gt_exporter_.export_array(output, this->debug_folder_full_path_ + "/" + str.str());
                }

                for (size_t phs = 0; phs < PHS; phs++)
                {
                    lax_highres_dicom(0, phs).header_.image_series_index *= 100;

                    lax_highres_dicom(0, phs).attrib_.append(GADGETRON_IMAGEPROCESSINGHISTORY, GADGETRON_AI);
                    lax_highres_dicom(0, phs).attrib_.append(GADGETRON_SEQUENCEDESCRIPTION, GADGETRON_AI);

                    Gadgetron::set_attrib_from_ismrmrd_header(lax_highres_dicom(0, phs).header_, lax_highres_dicom(0, phs).attrib_);
                }

                GDEBUG_STREAM("attach landmarks to attributes ...");
                this->attach_info_to_report(lax_highres_dicom, pts);
            }

            GDEBUG_STREAM("prepare outputs ...");
            this->convert_image_container_to_array(lax_highres_dicom, lax_ai);

            if (!this->debug_folder_full_path_.empty())
            {
                std::stringstream str;
                str << "lax_ai";
                gt_exporter_.export_array_complex(lax_ai.data_, this->debug_folder_full_path_ + "/" + str.str());
            }
        }
        catch (...)
        {
            GERROR_STREAM("Exceptions happened in CmrRealTimeLAXCineAIAnalysisGadget::perform_LAX_detection_AI(...) ... ");
            return GADGET_FAIL;
        }

        return GADGET_OK;
    }

    // ----------------------------------------------------------------------------------------

    GADGET_FACTORY_DECLARE(CmrRealTimeLAXCineAIAnalysisGadget)
}
