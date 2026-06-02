
#include "CmrPSIRNetGadget.h"
#include "hoNDImage_util.h"
#include "hoNDArray_reductions.h"
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
                if (!boost::algorithm::ends_with(model_name, ".pts"))
                {
                    model_name += ".pts";
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

            // PSIRNet is a single-shot model: each call takes one IR + one PD k-space
            // (see psirnet/src/models.py:PSIRNet). We invoke it independently for every
            // (average, slice) pair so that downstream tools / cardiologists can pick
            // any single-shot reconstruction retrospectively (e.g. compare against a
            // MOCO multi-average reference). The N dimension is preserved end-to-end.
            //
            // To keep the Python wrapper agnostic to the meaning of the batch dim, we
            // pack (n, slc) into a flat batch index `b = n + slc*N` for the model inputs
            // and unpack the same way on the outputs.
            res_psir_.data_.create(RO, E1, E2, 1, N, 1, SLC);
            res_magir_.data_.create(RO, E1, E2, 1, N, 1, SLC);

            const size_t B = N * SLC;

            // get IR and PD k-space data, flattened along batch = n + slc*N
            hoNDArray< std::complex<float> > kspace_ir;
            kspace_ir.create(RO, E1, CHA, B);
            Gadgetron::clear(kspace_ir);

            hoNDArray< std::complex<float> > kspace_pd;
            kspace_pd.create(RO, E1, CHA, B);
            Gadgetron::clear(kspace_pd);

            size_t ro, e1, cha, slc, n;
            for (slc=0; slc<SLC; slc++)
            {
                for (n=0; n<N; n++)
                {
                    const size_t b = n + slc*N;
                    for (cha=0; cha<CHA; cha++)
                    {
                        for (e1=0; e1<E1; e1++)
                        {
                            for (ro=0; ro<RO; ro++)
                            {
                                kspace_ir(ro, e1, cha, b) = recon_bit.data_.data_(ro, e1, 0, cha, n, 0, slc);
                                kspace_pd(ro, e1, cha, b) = recon_bit.data_.data_(ro, e1, 0, cha, n, 1, slc);
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
            // The coil map is one per slice (independent of average), so we tile it
            // along the same flat batch dim used for the k-space inputs.
            hoNDArray< std::complex<float> > coil_map_model;
            coil_map_model.create(RO, E1, CHA, B);
            Gadgetron::clear(coil_map_model);
            for (slc=0; slc<SLC; slc++)
            {
                for (n=0; n<N; n++)
                {
                    const size_t b = n + slc*N;
                    for (cha=0; cha<CHA; cha++)
                    {
                        for (e1=0; e1<E1; e1++)
                        {
                            for (ro=0; ro<RO; ro++)
                            {
                                coil_map_model(ro, e1, cha, b) = coil_map(ro, e1, 0, cha, slc);
                            }
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
            hoNDArray< float > psir; // (RO, E1, 1, B) real-valued PSIR
            {
                GILLock lg;
                PythonFunction< hoNDArray<float> > apply_psirnet("psirnet", "apply_psirnet");
                psir = apply_psirnet(kspace_ir, kspace_pd, coil_map_model, this->model_);
            }
            if (perform_timing.value()) { timer.stop(); }

            if (!debug_folder_full_path_.empty())
            {
                gt_exporter_.export_array(psir, debug_folder_full_path_ + "psir" + os.str());
            }

            // set the results: one PSIR + one MagIR image per (slice, average)
            res_psir_.headers_.create(N, 1, SLC);
            res_psir_.meta_.resize(N * SLC);

            res_magir_.headers_.create(N, 1, SLC);
            res_magir_.meta_.resize(N * SLC);

            // Unpack the flat-batch real-valued psir into the (RO,E1,E2,CHA,N,S,SLC)
            // complex IsmrmrdImageArrays, applying the scale factor inline.
            // PSIR stays signed (sign carried by image_type=REAL downstream);
            // MagIR is |PSIR|. The additive psir_offset is applied later in
            // compute_image_header_psir_magir, on res_psir_ only.
            const float scale_factor = this->scale_factor_after_SCC.value();
            for (slc=0; slc<SLC; slc++)
            {
                for (n=0; n<N; n++)
                {
                    const size_t b = n + slc*N;
                    for (e1=0; e1<E1; e1++)
                    {
                        for (ro=0; ro<RO; ro++)
                        {
                            const float v = psir(ro, e1, 0, b) * scale_factor;
                            res_psir_.data_(ro, e1, 0, 0, n, 0, slc)  = std::complex<float>(v, 0.0f);
                            res_magir_.data_(ro, e1, 0, 0, n, 0, slc) = std::complex<float>(std::abs(v), 0.0f);
                        }
                    }
                }
            }

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
            float thres =  10;

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

            hoNDArray<float> magIRImage_float;
            magIRImage_float.create(RO, E1);

            hoNDArray<std::complex<float>> PSIRImage;
            PSIRImage.create(RO, E1);

            // PSIR is a real-valued, signed image. The scanner only consumes unsigned shorts,
            // so we shift the PSIR pixel values into a strictly non-negative range at the very
            // end of this function. Until then, everything downstream of `perform_psir` must
            // see the signed values. Two things are required to make the rest of the chain do
            // the right thing with PSIR:
            //   1. `ComplexToFloatGadget` branches on `header.image_type`: MAGNITUDE -> abs(),
            //      REAL -> real(). The base-class `compute_image_header(...)` hard-codes
            //      MAGNITUDE, which would strip the sign. We tag PSIR as REAL here.
            //   2. We record the additive offset in `GADGETRON_IMAGE_SCALE_OFFSET` meta so the
            //      receiving viewer/DICOM stage can recover the signed value.
            // MagIR is left as MAGNITUDE (it is |IR|, non-negative) and is not shifted.
            // Sourced from the `psir_offset` gadget property so the actual additive shift, the
            // SCALE_OFFSET meta, and the window-centre adjustment all stay in lock-step.
            const double psir_offset = this->psir_offset.value();

            // loop through the headers and meta and set fields
            size_t n, s, slc;
            for (slc=0; slc<SLC; slc++)
            {
                for (n=0; n<N; n++)
                {
                    // copy the header and meta information from the PSIR result to the MAGIR result
                    res_magir.headers_(n, 0, slc) = res_psir.headers_(n, 0, slc);
                    res_magir.meta_[n + slc*N] = res_psir.meta_[n + slc*N];

                    // Override the MAGNITUDE default planted by GenericReconGadget::compute_image_header.
                    res_psir.headers_(n, 0, slc).image_type  = ISMRMRD::ISMRMRD_IMTYPE_REAL;
                    res_magir.headers_(n, 0, slc).image_type = ISMRMRD::ISMRMRD_IMTYPE_MAGNITUDE;

                    res_psir.meta_[n + slc*N].set(GADGETRON_IMAGE_SCALE_RATIO, 1.0);
                    res_psir.meta_[n + slc*N].set(GADGETRON_IMAGECOMMENT, GADGETRON_IMAGE_PSIR);
                    res_psir.meta_[n + slc*N].set(GADGETRON_SEQUENCEDESCRIPTION, GADGETRON_IMAGE_PSIR);
                    res_psir.meta_[n + slc*N].set(GADGETRON_DATA_ROLE, GADGETRON_IMAGE_PSIR);
                    res_psir.meta_[n + slc*N].set(GADGETRON_IMAGE_SCALE_OFFSET, (double)psir_offset);
                    if(!TI_.empty()) res_psir.meta_[n + slc*N].set(GADGETRON_IMAGE_INVERSIONTIME, TI_[0]);

                    res_magir.meta_[n + slc*N].set(GADGETRON_IMAGE_SCALE_RATIO, 1.0);
                    res_magir.meta_[n + slc*N].set(GADGETRON_IMAGECOMMENT, GADGETRON_IMAGE_MAGIR);
                    res_magir.meta_[n + slc*N].set(GADGETRON_SEQUENCEDESCRIPTION, GADGETRON_IMAGE_MAGIR);
                    res_magir.meta_[n + slc*N].set(GADGETRON_DATA_ROLE, GADGETRON_IMAGE_MAGIR);
                    if(!TI_.empty()) res_magir.meta_[n + slc*N].set(GADGETRON_IMAGE_INVERSIONTIME, TI_[0]);

                    // compute window level for the PSIR image on the signed (un-shifted) values.
                    memcpy(PSIRImage.begin(), &res_psir.data_(0,0,0,0,n,0,slc), sizeof(std::complex<float>)*RO*E1);
                    memcpy(magIRImage.begin(), &res_magir.data_(0,0,0,0,n,0,slc), sizeof(std::complex<float>)*RO*E1);

                    float window_center = 2048;
                    float window_width = 1200;
                    if (!calculate_window_level(magIRImage, PSIRImage, window_center, window_width))
                    {
                        GERROR_STREAM("Failed to calculate window level for PSIR image");
                    }
                    else
                    {
                        GDEBUG_STREAM("Calculated window level for PSIR image " << n << " in slice " << slc << " : window center = " << window_center << " , window width = " << window_width);
                        // The pixel data is shifted by `psir_offset` below, so the window centre
                        // must shift by the same amount to remain valid for the displayed pixels.
                        res_psir.meta_[n + slc*N].set(GADGETRON_IMAGE_WINDOWCENTER, (double)(window_center + psir_offset));
                        res_psir.meta_[n + slc*N].set(GADGETRON_IMAGE_WINDOWWIDTH, (double)window_width);
                    }

                    // compute magIR windowing
                    Gadgetron::abs(magIRImage, magIRImage_float);
                    float fraction = 0.95;
                    float high_end = Gadgetron::percentile(magIRImage_float, fraction);
                    window_center = (high_end - 0.02 * this->scale_factor_after_SCC.value()) / 2;
                    res_magir.meta_[n + slc*N].set(GADGETRON_IMAGE_WINDOWCENTER, (double)window_center);
                    res_magir.meta_[n + slc*N].set(GADGETRON_IMAGE_WINDOWWIDTH, (double)window_width);
                    GDEBUG_STREAM("Calculated window level for MagIR image " << n << " in slice " << slc << " : window center = " << window_center << " , window width = " << window_width);
                }
            }

            // Final positive-shift so the unsigned-short conversion downstream preserves the sign.
            // Only the PSIR image is shifted; MagIR is already non-negative.
            res_psir.data_ += static_cast<double>(psir_offset);
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