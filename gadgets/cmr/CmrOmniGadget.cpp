
#include "CmrOmniGadget.h"
#include "hoNDImage_util.h"
#include "hoNDArray_reductions.h"
#include <boost/algorithm/string.hpp>
#include <sstream>

namespace Gadgetron {

    CmrOmniGadget::CmrOmniGadget() : BaseClass()
    {
        num_of_PD_images_ = 0;
    }

    CmrOmniGadget::~CmrOmniGadget()
    {
    }

    int CmrOmniGadget::process_config(ACE_Message_Block* mb)
    {
        GADGET_CHECK_RETURN(BaseClass::process_config(mb) == GADGET_OK, GADGET_FAIL);

        try
        {
            std::string gadgetron_home = this->context.paths.gadgetron_home.generic_string();
            boost::filesystem::path gadgetron_python_path = this->context.paths.gadgetron_home / "share" / "gadgetron" / "python";

            Gadgetron::initialize_python();
            Gadgetron::add_python_path(gadgetron_python_path.generic_string());
            this->gt_home_ = gadgetron_python_path.generic_string();

            boost::filesystem::path cmr_python_path = this->context.paths.gadgetron_home / "share" / "gadgetron" / "python" / "cmr_ml";
            Gadgetron::add_python_path(cmr_python_path.generic_string());

            boost::filesystem::path model_path = this->context.paths.gadgetron_home / "share" / "gadgetron" / "python" / "cmr_ml" / "models";
            Gadgetron::add_python_path(model_path.generic_string());
            this->model_dir_ = model_path.generic_string();

            GDEBUG_STREAM("Set up python path using context : " << this->gt_home_);
            GDEBUG_STREAM("Set up model path using context : " << this->model_dir_);
        }
        catch (...)
        {
            GERROR_STREAM("Exception happened when adding  path to python ... ");
            return GADGET_FAIL;
        }

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

        if (!this->prepare_AI()) { return GADGET_FAIL; }

        if (h.userParameters)
        {
            for (std::vector<ISMRMRD::UserParameterLong>::const_iterator i = h.userParameters->userParameterLong.begin(); i != h.userParameters->userParameterLong.end(); ++i)
            {
                if (this->file_type.value() == "Perfusion")
                {
                    if (std::strcmp(i->name.c_str(), "NumOfProtonDensityImages") == 0)
                    {
                        num_of_PD_images_ = i->value;                        
                    }
                }
            }
        }
        GDEBUG_STREAM("NumOfProtonDensityImages is " << num_of_PD_images_);
        // -------------------------------------------------

        return GADGET_OK;
    }
    
    bool CmrOmniGadget::prepare_AI()
    {
        std::string gt_home;
        char* v = std::getenv("GADGETRON_HOME");
        if (v == NULL)
        {
#ifdef _WIN32
            gt_home = "D:/gtuser/mrprogs/install/";
#else
            gt_home = "/usr/local/";
#endif // _WIN32
        }
        else
        {
            gt_home = std::string(v);
        }

        gt_home.append("/share/gadgetron/python/");
        GDEBUG_STREAM("Gadgetron python directory : " << gt_home);

        Gadgetron::initialize_python();

        GDEBUG_STREAM("Omninet model file : " << this->model.value());

        try
        {
            if (!this->gt_home_.empty())
            {
                std::string model_name = this->model.value();
                if (!boost::algorithm::ends_with(model_name, ".pts"))
                {
                    model_name += ".pts";
                }

                GDEBUG_STREAM("Load OmniNet model : " << model_name);

                if (this->perform_timing.value()) { gt_timer_.start("OmniNet model loading ... "); }
                {
                    GILLock lg;
                    PythonFunction<boost::python::object> load_model_for_inference("omninet", "load_model_for_inference");
                    model_ = load_model_for_inference(this->model_dir_, model_name);
                    bp::incref(model_.ptr());
                }
                if (this->perform_timing.value()) { gt_timer_.stop(); }

                GDEBUG_STREAM("Load OmniNet model completed ... ");
            }
        }
        catch (...)
        {
            GERROR_STREAM("Loading OmniNet model failed ... ");
            return false;
        }

        return true;
    }

    int CmrOmniGadget::process(Gadgetron::GadgetContainerMessage< IsmrmrdReconData >* m1)
    {
        if (perform_timing.value()) { gt_timer_local_.start("CmrOmniGadget::process"); }

        process_called_times_++;

        IsmrmrdReconData* recon_bit_ = m1->getObjectPtr();
        if (recon_bit_->rbit_.size() > num_encoding_spaces_)
        {
            GWARN_STREAM("Incoming recon_bit has more encoding spaces than the protocol : " << recon_bit_->rbit_.size() << " instead of " << num_encoding_spaces_);
        }

        // for every encoding space
        for (size_t e = 0; e < recon_bit_->rbit_.size(); e++)
        {
            std::stringstream os;
            os << "_encoding_" << e;

            GDEBUG_CONDITION_STREAM(verbose.value(), "Calling " << process_called_times_ << " , encoding space : " << e);
            GDEBUG_CONDITION_STREAM(verbose.value(), "======================================================================");

            // ---------------------------------------------------------------
            // export incoming data
            if (!debug_folder_full_path_.empty())
            {
                gt_exporter_.export_array_complex(recon_bit_->rbit_[e].data_.data_, debug_folder_full_path_ + "data" + os.str());
            }

            if (!debug_folder_full_path_.empty())
            {
                gt_exporter_.export_array_complex(recon_bit_->rbit_[e].ref_->data_, debug_folder_full_path_ + "ref" + os.str());
            }

            // ---------------------------------------------------------------

            if (recon_bit_->rbit_[e].data_.data_.get_number_of_elements() > 0)
            {
                // ---------------------------------------------------------------

                if (perform_timing.value()) { gt_timer_.start("CmrOmniGadget::perform_omni_recon"); }
                this->perform_omni_recon(recon_bit_->rbit_[e], e);
                if (perform_timing.value()) { gt_timer_.stop(); }

                // ---------------------------------------------------------------

                if (perform_timing.value()) { gt_timer_.start("CmrOmniGadget::compute_image_header"); }
                this->compute_image_header(recon_bit_->rbit_[e], res_, e);
                if (perform_timing.value()) { gt_timer_.stop(); }

                // ---------------------------------------------------------------

                if (!debug_folder_full_path_.empty())
                {
                    this->gt_exporter_.export_array_complex(res_.data_, debug_folder_full_path_ + "recon_res" + os.str());
                }

                if (perform_timing.value()) { gt_timer_.start("CmrOmniGadget::send_out_image_array, res"); }
                this->send_out_image_array(res_, e, image_series.value() + ((int)e + 111), file_type.value() == "Retro" ? GADGETRON_IMAGE_RETRO : GADGETRON_IMAGE_REGULAR);
                if (perform_timing.value()) { gt_timer_.stop(); }
            }
        }

        m1->release();

        if (perform_timing.value()) { gt_timer_local_.stop(); }

        return GADGET_OK;
    }

    void CmrOmniGadget::perform_omni_recon(IsmrmrdReconBit& recon_bit, size_t encoding)
    {
        try
        {
            size_t RO  = recon_bit.data_.data_.get_size(0);
            size_t E1  = recon_bit.data_.data_.get_size(1);
            size_t E2  = recon_bit.data_.data_.get_size(2);
            size_t CHA = recon_bit.data_.data_.get_size(3);
            size_t N   = recon_bit.data_.data_.get_size(4);
            size_t S   = recon_bit.data_.data_.get_size(5);
            size_t SLC = recon_bit.data_.data_.get_size(6);

            std::stringstream os;
            os << "_encoding_" << encoding;

            GADGET_CHECK_THROW(E2==1);

            GDEBUG_STREAM("OmniNet recon for encoding space " << encoding << " with matrix size : [" << RO << "," << E1 << "," << E2 << "] , #coils : " << CHA << " , #ave : " << N << " , #sets : " << S << " , #slices : " << SLC);
            GDEBUG_STREAM("data is " << Gadgetron::nrm2(recon_bit.data_.data_) << " , ref is " << Gadgetron::nrm2(recon_bit.ref_->data_));

            Gadgetron::GadgetronTimer timer(false);

            // compute coil map
            if (perform_timing.value()) { timer.start("compute coil map ... "); }

            hoNDArray< std::complex<float> > ref_calib, ref_coil_map;
            this->make_ref_coil_map(*recon_bit.ref_, recon_bit.data_.data_.get_dimensions(), ref_calib, ref_coil_map, encoding);

            if (!debug_folder_full_path_.empty()) 
            { 
                gt_exporter_.export_array_complex(ref_calib, debug_folder_full_path_ + "ref_calib" + os.str()); 
                gt_exporter_.export_array_complex(ref_coil_map, debug_folder_full_path_ + "ref_coil_map" + os.str()); 
            }

            hoNDArray< std::complex<float> > coil_map;
            this->perform_coil_map_estimation(ref_coil_map, coil_map, encoding);

            if (!debug_folder_full_path_.empty()) 
            { 
                gt_exporter_.export_array_complex(coil_map, debug_folder_full_path_ + "coil_map" + os.str()); 
            }
            if (perform_timing.value()) { timer.stop(); }

            res_.data_.create(RO, E1, E2, 1, N, S, SLC);

            if (perform_timing.value()) { timer.start("compute omni net model ... "); }

            if ((S==1) && (N>=1))
            {
                hoNDArray< std::complex<float> > kspace;
                kspace.create(RO, E1, CHA, N, SLC);
                Gadgetron::clear(kspace);

                size_t ro, e1, cha, slc, n;
                for (slc=0; slc<SLC; slc++)
                {
                    for (n=0; n<N; n++)
                    {
                        for (cha=0; cha<CHA; cha++)
                        {
                            for (e1=0; e1<E1; e1++)
                            {
                                for (ro=0; ro<RO; ro++)
                                {
                                    kspace(ro, e1, cha, n, slc) = recon_bit.data_.data_(ro, e1, 0, cha, n, 0, slc);
                                }
                            }
                        }
                    }
                }

                if (!debug_folder_full_path_.empty()) 
                { 
                    gt_exporter_.export_array_complex(kspace, debug_folder_full_path_ + "kspace" + os.str()); 
                }

                hoNDArray< std::complex<float> > coil_map_model;
                coil_map_model.create(RO, E1, CHA, SLC);
                for (slc=0; slc<SLC; slc++)
                {
                    for (cha=0; cha<CHA; cha++)
                    {
                        for (e1=0; e1<E1; e1++)
                        {
                            for (ro=0; ro<RO; ro++)
                            {
                                coil_map_model(ro, e1, cha, slc) = coil_map(ro, e1, 0, cha, 0, 0, slc);
                            }
                        }
                    }
                }
                if (!debug_folder_full_path_.empty()) 
                { 
                    gt_exporter_.export_array_complex(coil_map_model, debug_folder_full_path_ + "coil_map_model" + os.str()); 
                }

                if (perform_timing.value()) { timer.start("compute omni net model ... "); }

                hoNDArray< std::complex<float> > res_kspace;
                {
                    GILLock lg;
                    PythonFunction< hoNDArray<std::complex<float>> > apply_omninet("omninet", "apply_omninet");
                    res_kspace = apply_omninet(kspace, coil_map_model, this->model_, file_type.value(), this->num_of_PD_images_);
                }
                if (perform_timing.value()) { timer.stop(); }

                if (!debug_folder_full_path_.empty())
                {
                    gt_exporter_.export_array_complex(res_kspace, debug_folder_full_path_ + "res_kspace" + os.str());
                }

                for (slc=0; slc<SLC; slc++)
                {
                    for (n=0; n<N; n++)
                    {
                        for (cha=0; cha<CHA; cha++)
                        {
                            for (e1=0; e1<E1; e1++)
                            {
                                for (ro=0; ro<RO; ro++)
                                {
                                    res_.data_(ro, e1, 0, cha, n, 0, slc) = res_kspace(ro, e1, cha, n, slc);
                                }
                            }
                        }
                    }
                }
            }

            if ((S>1) && (N==1))
            {
                hoNDArray< std::complex<float> > kspace;
                kspace.create(RO, E1, CHA, S, SLC);
                Gadgetron::clear(kspace);

                size_t ro, e1, cha, slc, s;
                for (slc=0; slc<SLC; slc++)
                {
                    for (s=0; s<S; s++)
                    {
                        for (cha=0; cha<CHA; cha++)
                        {
                            for (e1=0; e1<E1; e1++)
                            {
                                for (ro=0; ro<RO; ro++)
                                {
                                    kspace(ro, e1, cha, s, slc) = recon_bit.data_.data_(ro, e1, 0, cha, 0, s, slc);
                                }
                            }
                        }
                    }
                }

                if (!debug_folder_full_path_.empty()) 
                { 
                    gt_exporter_.export_array_complex(kspace, debug_folder_full_path_ + "kspace" + os.str()); 
                }

                hoNDArray< std::complex<float> > coil_map_model;
                coil_map_model.create(RO, E1, CHA, SLC);
                for (slc=0; slc<SLC; slc++)
                {
                    for (cha=0; cha<CHA; cha++)
                    {
                        for (e1=0; e1<E1; e1++)
                        {
                            for (ro=0; ro<RO; ro++)
                            {
                                coil_map_model(ro, e1, cha, slc) = coil_map(ro, e1, 0, cha, 0, 0, slc);
                            }
                        }
                    }
                }
                if (!debug_folder_full_path_.empty()) 
                { 
                    gt_exporter_.export_array_complex(coil_map_model, debug_folder_full_path_ + "coil_map_model" + os.str()); 
                }

                if (perform_timing.value()) { timer.start("compute omni net model ... "); }

                hoNDArray< std::complex<float> > res_kspace;
                {
                    GILLock lg;
                    PythonFunction< hoNDArray<std::complex<float>> > apply_omninet("omninet", "apply_omninet");
                    res_kspace = apply_omninet(kspace, coil_map_model, this->model_, file_type.value(), this->num_of_PD_images_);
                }
                if (perform_timing.value()) { timer.stop(); }

                if (!debug_folder_full_path_.empty())
                {
                    gt_exporter_.export_array_complex(res_kspace, debug_folder_full_path_ + "res_kspace" + os.str());
                }

                for (slc=0; slc<SLC; slc++)
                {
                    for (s=0; s<S; s++)
                    {
                        for (cha=0; cha<CHA; cha++)
                        {
                            for (e1=0; e1<E1; e1++)
                            {
                                for (ro=0; ro<RO; ro++)
                                {
                                    res_.data_(ro, e1, 0, cha, 0, s, slc) = res_kspace(ro, e1, cha, s, slc);
                                }
                            }
                        }
                    }
                }
            }

            res_.headers_.create(N, S, SLC);;
            res_.meta_.resize(N*S*SLC);
        }
        catch (...)
        {
            GADGET_THROW("Errors happened in CmrOmniGadget::perform_omni_recon(...) ... ");
        }
    }

    GADGET_FACTORY_DECLARE(CmrOmniGadget)
}