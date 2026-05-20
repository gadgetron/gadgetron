
#include "CmrPSIRNetGadget.h"
#include "hoNDImage_util.h"
#include <boost/algorithm/string.hpp>
#include <sstream>

namespace Gadgetron {

    CmrPSIRNetGadget::CmrPSIRNetGadget() : BaseClass()
    {
    }

    CmrPSIRNetGadget::~CmrPSIRNetGadget()
    {
    }

    int CmrPSIRNetGadget::process_config(ACE_Message_Block* mb)
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

        if (h.sequenceParameters.is_present())
        {
            if (h.sequenceParameters.get().TI.is_present())
            {
                TI_ = h.sequenceParameters.get().TI.get();
            }
        }
        else
        {
            GWARN_STREAM("Inversion time does not exist in the seq protocols ... ");
        }

        if (!this->prepare_AI()) { return GADGET_FAIL; }

        // -------------------------------------------------

        return GADGET_OK;
    }
    
    bool CmrPSIRNetGadget::prepare_AI()
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

        GDEBUG_STREAM("PSIRNet model file : " << this->model.value());

        try
        {
            if (!this->gt_home_.empty())
            {
                std::string model_name = this->model.value();
                if (!boost::algorithm::ends_with(this->model.value(), "pts"))
                {
                    std::string model_name = this->model.value() + ".pts";
                }

                GDEBUG_STREAM("Load PSIRNet model : " << model_name);

                if (this->perform_timing.value()) { gt_timer_.start("PSIRNet model loading ... "); }
                {
                    GILLock lg;
                    PythonFunction<boost::python::object> load_model("psirnet", "load_model_for_inference");
                    model_ = load_model(this->model_dir_, model_name);
                    bp::incref(model_.ptr());
                }
                if (this->perform_timing.value()) { gt_timer_.stop(); }

                GDEBUG_STREAM("Load PSIRNet model completed ... ");
            }
        }
        catch (...)
        {
            GERROR_STREAM("Loading PSIRNet model failed ... ");
            return false;
        }

        return true;
    }

    int CmrPSIRNetGadget::process(Gadgetron::GadgetContainerMessage< IsmrmrdReconData >* m1)
    {
        if (perform_timing.value()) { gt_timer_local_.start("CmrPSIRNetGadget::process"); }

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

                if (perform_timing.value()) { gt_timer_.start("CmrPSIRNetGadget::perform_psir"); }
                this->perform_psir(recon_bit_->rbit_[e], e);
                if (perform_timing.value()) { gt_timer_.stop(); }

                // ---------------------------------------------------------------

                if (perform_timing.value()) { gt_timer_.start("CmrPSIRNetGadget::compute_image_header, psir images"); }
                this->compute_image_header(recon_bit_->rbit_[e], res_psir_, e);
                if (perform_timing.value()) { gt_timer_.stop(); }

                // ---------------------------------------------------------------
                // adjust headers for psir and mag IR
                this->compute_image_header_psir_magir(res_psir_, res_magir_, e);
                // ---------------------------------------------------------------

                if (!debug_folder_full_path_.empty())
                {
                    this->gt_exporter_.export_array_complex(res_psir_.data_, debug_folder_full_path_ + "recon_res" + os.str());
                    this->gt_exporter_.export_array_complex(res_magir_.data_, debug_folder_full_path_ + "recon_res_magir" + os.str());
                }

                if (perform_timing.value()) { gt_timer_.start("CmrPSIRNetGadget::send_out_image_array, magir"); }
                this->send_out_image_array(res_magir_, e, image_series.value() + ((int)e + 109), GADGETRON_IMAGE_MAGIR);
                if (perform_timing.value()) { gt_timer_.stop(); }

                if (perform_timing.value()) { gt_timer_.start("CmrPSIRNetGadget::send_out_image_array, psir"); }
                this->send_out_image_array(res_psir_, e, image_series.value() + ((int)e + 111), GADGETRON_IMAGE_PSIR);
                if (perform_timing.value()) { gt_timer_.stop(); }
            }
        }

        m1->release();

        if (perform_timing.value()) { gt_timer_local_.stop(); }

        return GADGET_OK;
    }

    void CmrPSIRNetGadget::perform_psir(IsmrmrdReconBit& recon_bit, size_t encoding)
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

            GADGET_CHECK_THROW(E2==1);
            GADGET_CHECK_THROW(S==2);

            GDEBUG_STREAM("PSIRNet recon for encoding space " << encoding << " with matrix size : [" << RO << "," << E1 << "," << E2 << "] , #coils : " << CHA << " , #ave : " << N << " , #sets : " << S << " , #slices : " << SLC);
            GDEBUG_STREAM("data is " << Gadgetron::nrm2(recon_bit.data_.data_) << " , ref is " << Gadgetron::nrm2(recon_bit.ref_->data_));

            Gadgetron::GadgetronTimer timer(false);

            // set up outputs
            res_psir_.data_.create(RO, E1, E2, 1, N, 1, SLC);
            res_magir_.data_.create(RO, E1, E2, 1, N, 1, SLC);

            // get IR and PD k-space data
            hoNDArray< std::complex<float> > kspace_ir;
            kspace_ir.create(RO, E1, CHA, N*SLC);
            Gadgetron::clear(kspace_ir);

            hoNDArray< std::complex<float> > kspace_pd;
            kspace_pd.create(RO, E1, CHA, N*SLC);
            Gadgetron::clear(kspace_pd);

            size_t ro, e1, cha, n, s, slc;
            for (slc=0; slc<SLC; slc++)
            {
                for (n=0; n<N; n++)
                {
                    //memcpy(&kspace_ir(0, 0, 0, n + slc*N), &recon_bit.data_.data_(0, 0, 0, 0, n, 0, slc), RO*E1*CHA*sizeof(std::complex<float>));
                    //memcpy(&kspace_pd(0, 0, 0, n + slc*N), &recon_bit.data_.data_(0, 0, 0, 0, n, 1, slc), RO*E1*CHA*sizeof(std::complex<float>));

                    for (cha=0; cha<CHA; cha++)
                    {
                        for (e1=0; e1<E1; e1++)
                        {
                            for (ro=0; ro<RO; ro++)
                            {             
                                kspace_ir(ro, e1, cha, n + slc*N) = recon_bit.data_.data_(ro, e1, 0, cha, n, 0, slc);
                                kspace_pd(ro, e1, cha, n + slc*N) = recon_bit.data_.data_(ro, e1, 0, cha, n, 1, slc);
                            }
                        }
                    }
                }
            }

            GDEBUG_STREAM("kspace_ir is " << Gadgetron::nrm2(kspace_ir) << " , kspace_pd is " << Gadgetron::nrm2(kspace_pd));

            std::stringstream os;
            os << "_encoding_" << encoding;

            if (!debug_folder_full_path_.empty()) 
            { 
                gt_exporter_.export_array_complex(kspace_ir, debug_folder_full_path_ + "kspace_ir" + os.str()); 
            }

            if (!debug_folder_full_path_.empty()) 
            { 
                gt_exporter_.export_array_complex(kspace_pd, debug_folder_full_path_ + "kspace_pd" + os.str()); 
            }

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

            // call the model
            hoNDArray< std::complex<float> > coil_map_model;
            coil_map_model.create(RO, E1, CHA, SLC);
            Gadgetron::clear(coil_map_model);
            for (slc=0; slc<SLC; slc++)
            {
                // memcpy(&coil_map_model(0, 0, 0, slc), &coil_map(0, 0, 0, slc), RO*E1*CHA*sizeof(std::complex<float>));
                for (cha=0; cha<CHA; cha++)
                {
                    for (e1=0; e1<E1; e1++)
                    {
                        for (ro=0; ro<RO; ro++)
                        {             
                            coil_map_model(ro, e1, cha, slc) = coil_map(ro, e1, 0, cha, slc);
                        }
                    }
                }
            }

            GDEBUG_STREAM("coil_map_model is " << Gadgetron::nrm2(coil_map_model));

            if (!debug_folder_full_path_.empty()) 
            { 
                gt_exporter_.export_array_complex(coil_map_model, debug_folder_full_path_ + "coil_map_model" + os.str()); 
            }

            if (perform_timing.value()) { timer.start("compute psir ... "); }
            hoNDArray< float > psir; // RO, E1, 1, N*SLC
            {
                GILLock lg;
                PythonFunction< hoNDArray< float > > apply_psirnet("psirnet", "apply_psirnet");
                psir = apply_psirnet(kspace_ir, kspace_pd, coil_map_model, this->model_);
            }
            if (perform_timing.value()) { timer.stop(); }

            if (!debug_folder_full_path_.empty()) 
            { 
                gt_exporter_.export_array(psir, debug_folder_full_path_ + "psir" + os.str()); 
            }

            // set the results
            res_psir_.headers_.create(N, 1, SLC);
            res_psir_.meta_.resize(N*SLC);

            res_magir_.headers_.create(N, 1, SLC);
            res_magir_.meta_.resize(N*SLC);

            float scale_factor = this->scale_factor_after_SCC.value();

            for (slc=0; slc<SLC; slc++)
            {
                for (n=0; n<N; n++)
                {             
                    for (e1=0; e1<E1; e1++)
                    {
                        for (ro=0; ro<RO; ro++)
                        { 
                            res_psir_.data_(ro, e1, 0, 0, n, 0, slc) = psir(ro, e1, 0, n + slc*N) * scale_factor + this->offset_factor_after_SCC.value();
                        }
                    }
                }
            }

            // compute magnitude image
            Gadgetron::abs(res_psir_.data_, res_magir_.data_);

            if (!debug_folder_full_path_.empty()) { gt_exporter_.export_array_complex(res_psir_.data_, debug_folder_full_path_ + "psir_" + os.str()); }
            if (!debug_folder_full_path_.empty()) { gt_exporter_.export_array_complex(res_magir_.data_, debug_folder_full_path_ + "magir_" + os.str()); }
        }
        catch (...)
        {
            GADGET_THROW("Errors happened in CmrPSIRNetGadget::perform_psir(...) ... ");
        }
    }

    bool CmrPSIRNetGadget::calculate_window_level(hoNDArray<std::complex<float>>& magIRImage, hoNDArray<std::complex<float>>& PSIRImage, float& window_center, float& window_width)
    {
        try
        {
            // since gfactor is not taken into account here, we need higher thresholding ratio
            float thres =  2;

            size_t N = PSIRImage.get_number_of_elements();

            std::vector<size_t> dim;
            PSIRImage.get_dimensions(dim);
            hoNDArray<float> mask(dim);
            Gadgetron::clear(mask);

            // get the foreground
            long long n;
            size_t numOfPixelInMask = 0;
            for ( n=0; n<N; n++ )
            {
                if ( magIRImage(n).real() > thres )
                {
                    mask(n) = 1;
                    numOfPixelInMask++;
                }
            }

            if ( numOfPixelInMask == 0 ) return true;

            // get the median of foreground
            std::vector<float> valueInMask(numOfPixelInMask, 0);
            size_t ind(0);
            for ( n=0; n<N; n++ )
            {
                if ( mask(n) == 1 )
                {
                    valueInMask[ind++] = PSIRImage(n).real();
                }
            }

            std::sort(valueInMask.begin(), valueInMask.end());
            float medianValueInMask = valueInMask[numOfPixelInMask/2];
            float windowing_high_end_percentile = 0.95;
            float w_high = valueInMask[ (size_t)(numOfPixelInMask *windowing_high_end_percentile) ];

            // get the second median and normal level, should be the center of myocardium
            hoNDArray<float> thresdhold(dim);
            Gadgetron::clear(thresdhold);

            for ( n=0; n<N; n++ )
            {
                if ( (PSIRImage(n).real()<medianValueInMask) && (mask(n) == 1) )
                {
                    thresdhold(n) = 1;
                }
            }            

            float normal_level = 0;

            numOfPixelInMask = 0;
            for ( n=0; n<N; n++ )
            {
                if ( thresdhold(n) == 1 )
                {
                    normal_level += PSIRImage(n).real();
                    numOfPixelInMask++;
                }
            }
            normal_level /= numOfPixelInMask;

            std::vector<float> valueInMaskMyo;
            valueInMaskMyo.resize(numOfPixelInMask, 0);

            ind = 0;
            for ( n=0; n<N; n++ )
            {
                if ( thresdhold(n) == 1 )
                {
                    valueInMaskMyo[ind++] = PSIRImage(n).real();
                }
            }
            std::nth_element(valueInMaskMyo.begin(), valueInMaskMyo.begin() + valueInMaskMyo.size() / 2, valueInMaskMyo.end());
            medianValueInMask = valueInMaskMyo[numOfPixelInMask / 2];

            // get the windowing setting

            float range = 1.1 * w_high - normal_level;
            float min_level = w_high - range;

            window_center = min_level + range / 2;
            window_width = range;
        }
        catch(...)
        {
            GERROR_STREAM("Error happened in calculate_window_level(hoNDArray<std::complex<float>>& magPDFiltered, hoNDArray<std::complex<float>>& PSIRImage, float& window_center, float& window_width) ... ");
            return false;
        }

        return true;
    }

    int CmrPSIRNetGadget::compute_image_header_psir_magir(IsmrmrdImageArray& res_psir, IsmrmrdImageArray& res_magir, size_t encoding)
    {
        try
        {
            size_t RO = res_psir.data_.get_size(0);
            size_t E1 = res_psir.data_.get_size(1);
            size_t E2 = res_psir.data_.get_size(2);
            size_t CHA = res_psir.data_.get_size(3);
            size_t N = res_psir.data_.get_size(4);
            size_t S = res_psir.data_.get_size(5);
            size_t SLC = res_psir.data_.get_size(6);

            res_magir.headers_ = res_psir.headers_;
            res_magir.meta_ = res_psir.meta_;

            hoNDArray<std::complex<float>> magIRImage;
            magIRImage.create(RO, E1);

            hoNDArray<std::complex<float>> PSIRImage;
            PSIRImage.create(RO, E1);

            // loop through the headers and meta and set fields
            size_t n, s, slc;
            for (slc=0; slc<SLC; slc++)
            {
                for (n=0; n<N; n++)
                {
                    // copy the header and meta information from the PSIR result to the MAGIR result
                    res_magir.headers_(n, 0, slc) = res_psir.headers_(n, 0, slc);
                    res_magir.meta_[n + slc*N] = res_psir.meta_[n + slc*N];

                    res_psir.meta_[n + slc*N].set(GADGETRON_IMAGE_SCALE_RATIO, 1.0);
                    res_psir.meta_[n + slc*N].set(GADGETRON_IMAGECOMMENT, GADGETRON_IMAGE_PSIR);
                    res_psir.meta_[n + slc*N].set(GADGETRON_SEQUENCEDESCRIPTION, GADGETRON_IMAGE_PSIR);
                    res_psir.meta_[n + slc*N].set(GADGETRON_DATA_ROLE, GADGETRON_IMAGE_PSIR);
                    if(!TI_.empty()) res_psir.meta_[n + slc*N].set(GADGETRON_IMAGE_INVERSIONTIME, TI_[0]);

                    res_magir.meta_[n + slc*N].set(GADGETRON_IMAGE_SCALE_RATIO, 1.0);
                    res_magir.meta_[n + slc*N].set(GADGETRON_IMAGECOMMENT, GADGETRON_IMAGE_MAGIR);
                    res_magir.meta_[n + slc*N].set(GADGETRON_SEQUENCEDESCRIPTION, GADGETRON_IMAGE_MAGIR);
                    res_magir.meta_[n + slc*N].set(GADGETRON_DATA_ROLE, GADGETRON_IMAGE_MAGIR);
                    if(!TI_.empty()) res_magir.meta_[n + slc*N].set(GADGETRON_IMAGE_INVERSIONTIME, TI_[0]);

                    // compute window level for the PSIR image
                    memcpy(PSIRImage.begin(), &res_psir.data_(0,0,0,0,n,0,slc), sizeof(std::complex<float>)*RO*E1);
                    memcpy(magIRImage.begin(), &res_magir.data_(0,0,0,0,n,0,slc), sizeof(std::complex<float>)*RO*E1);

                    float window_center = this->offset_factor_after_SCC.value();
                    float window_width = 1200;
                    if (!calculate_window_level(magIRImage, PSIRImage, window_center, window_width))
                    {
                        GERROR_STREAM("Failed to calculate window level for PSIR image");
                    }
                    else
                    {
                        GDEBUG_STREAM("Calculated window level for PSIR image " << n << " in slice " << slc << " : window center = " << window_center << " , window width = " << window_width);
                        res_psir.meta_[n + slc*N].set(GADGETRON_IMAGE_WINDOWCENTER, window_center);
                        res_psir.meta_[n + slc*N].set(GADGETRON_IMAGE_WINDOWWIDTH, window_width);
                    }
                }
            }
        }
        catch (...)
        {
            GERROR_STREAM("Errors in GenericReconGadget::compute_image_header_psir_magir(...) ... ");
            return GADGET_FAIL;
        }

        return GADGET_OK;
    }

    GADGET_FACTORY_DECLARE(CmrPSIRNetGadget)
}
