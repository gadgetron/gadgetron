/** \file   gtPlusPhaseSensitiveRecon.h
    \brief  Implement phase sensitive reconstruction
            The input is a 2D image container
    \author Hui Xue
*/

#pragma once

#include "gtPlusISMRMRDReconUtil.h"

namespace Gadgetron { namespace gtPlus {

    template<typename ValueType, unsigned int D> 
    class gtPlusPhaseSensitiveRecon
    {
    public:

        typedef gtPlusPhaseSensitiveRecon<ValueType, D> Self;

        typedef Gadgetron::hoNDImage<ValueType, D> ImageType;
        typedef Gadgetron::hoNDImage<ValueType, 2> Image2DType;
        typedef Gadgetron::hoNDImage<ValueType, 3> Image3DType;

        typedef ValueType T;
        typedef ValueType element_type;
        typedef ValueType value_type;

        typedef typename realType<T>::Type real_value_type;
        typedef Gadgetron::hoNDImage<real_value_type, D> ImageMagType;

        typedef Gadgetron::hoNDImageContainer2D<ImageType> ImageContinerType;

        gtPlusPhaseSensitiveRecon();
        virtual ~gtPlusPhaseSensitiveRecon();

        /// set the default parameters
        virtual bool setDefaultParameters();

        /// compute the phase sensitive recon
        virtual bool performPSIRRecon(const ImageContinerType& input);

        /// print the class information
        virtual void print(std::ostream& os) const;

        // ----------------------------------
        // parameters
        // ----------------------------------

        /// input image container containing images for processing
        ImageContinerType input_;

        /// row as the PD images
        unsigned int row_PD_;

        /// column used for PSIR windowing computation
        unsigned int col_compute_PSIR_windowing_;

        /// whether to perform PSIR
        bool perform_PSIR_;

        /// whether to perform surface coil correction using PD
        bool perform_SCC_;

        /// scaling factor after SCC
        real_value_type scale_factor_after_SCC_;

        /// whether to compute PSIR windowing
        bool compute_PSIR_windowing_;

        /// median filter window width
        unsigned int filter_width_[6];

        /// intensity scale factor for input images
        float intensity_scale_factor_;

        /// output parameters
        /// PSIR, magIR and magPD results
        ImageContinerType PSIR_;

        /// input image container containing images for processing
        ImageContinerType magIR_;

        /// input image container containing images for processing
        ImageContinerType magPD_;

        /// windowing 
        float window_center_;
        float window_width_;

        /// verbose mode
        bool verbose_;

        // ----------------------------------
        // debug and timing
        // ----------------------------------
        // clock for timing
        Gadgetron::GadgetronTimer gt_timer1_;
        Gadgetron::GadgetronTimer gt_timer2_;
        Gadgetron::GadgetronTimer gt_timer3_;

        bool performTiming_;

        // exporter
        Gadgetron::gtPlus::gtPlusIOAnalyze gt_exporter_;

        // debug folder
        std::string debugFolder_;

    protected:

        bool initialize(const ImageContinerType& input);

        /// compute a proper window level/width for PSIR images
        bool calculate_window_level(ImageType& magPDFiltered, ImageType& PSIRImage, ImageType& PSIRImageBeforeSCC, float& window_center, float& window_width);
    };

    template<typename ValueType, unsigned int D> 
    gtPlusPhaseSensitiveRecon<ValueType, D>::gtPlusPhaseSensitiveRecon() : performTiming_(false)
    {
        gt_timer1_.set_timing_in_destruction(false);
        gt_timer2_.set_timing_in_destruction(false);
        gt_timer3_.set_timing_in_destruction(false);

        GADGET_CHECK_THROW(this->setDefaultParameters());
    }

    template<typename ValueType, unsigned int D> 
    gtPlusPhaseSensitiveRecon<ValueType, D>::
    ~gtPlusPhaseSensitiveRecon()
    {
    }

    template<typename ValueType, unsigned int D> 
    bool gtPlusPhaseSensitiveRecon<ValueType, D>::setDefaultParameters()
    {
        row_PD_ = 1;
        col_compute_PSIR_windowing_ = 0;
        perform_PSIR_ = true;
        perform_SCC_ = true;
        scale_factor_after_SCC_ = 100;
        compute_PSIR_windowing_ = true;

        filter_width_[0] = 7;
        filter_width_[1] = 7;
        filter_width_[2] = 7;
        filter_width_[3] = 7;
        filter_width_[4] = 7;
        filter_width_[5] = 7;

        intensity_scale_factor_ = 1;

        window_center_ = -1;
        window_width_ = -1;

        verbose_ = false;

        return true;
    }

    template<typename ValueType, unsigned int D> 
    bool gtPlusPhaseSensitiveRecon<ValueType, D>::
    initialize(const ImageContinerType& input)
    {
        try
        {
            size_t row = input.rows();
            std::vector<size_t> cols = input.cols();

            // allocate results
            std::vector<size_t> dim(row-1, cols[0]);
            GADGET_CHECK_RETURN_FALSE(PSIR_.create(dim));
            GADGET_CHECK_RETURN_FALSE(magIR_.create(dim));

            dim.clear();
            dim.resize(1);
            dim[0] = cols[0];
            GADGET_CHECK_RETURN_FALSE(magPD_.create(dim));

            input_ = input;
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Error happened in gtPlusPhaseSensitiveRecon<ValueType, D>::initialize(const TargetContinerType& targetContainer) ... ");
            return false;
        }

        return true;
    }

    template<typename ValueType, unsigned int D> 
    bool gtPlusPhaseSensitiveRecon<ValueType, D>::performPSIRRecon(const ImageContinerType& input)
    {
        try
        {
            GADGET_CHECK_RETURN_FALSE(this->initialize(input));

            size_t row = input_.rows();
            std::vector<size_t> cols = input_.cols();

            size_t RO = input_(0, 0).get_size(0);
            size_t E1 = input_(0, 0).get_size(1);

            size_t r, c;

            ImageType PDMag(RO, E1), PDPhase(RO, E1);

            for ( c=0; c<cols[0]; c++ )
            {
                ImageType& PDImage = input_(row_PD_, c);

                Gadgetron::absolute(PDImage, PDMag);
                Gadgetron::addEpsilon(PDMag);

                ImageType& magPD = magPD_(0, c);
                magPD = PDMag;

                if ( perform_PSIR_ )
                {
                    Gadgetron::divide(PDImage, PDMag, PDPhase);

                    size_t rInd = 0;
                    for ( r=0; r<row; r++ )
                    {
                        if ( r == row_PD_ ) continue;

                        ImageType& IRImage = input_(r, c);
                        ImageType& PSIRImage = PSIR_(rInd, c);

                        PSIRImage = IRImage;

                        Gadgetron::multiplyConj(PSIRImage, PDPhase, PSIRImage);
                        Gadgetron::complex_to_real(PSIRImage, PSIRImage);

                        ImageType& magIRImage = magIR_(rInd, c);
                        magIRImage = PSIRImage;

                        Gadgetron::absolute(magIRImage, magIRImage);

                        rInd++;
                    }
                }
                else
                {
                    size_t rInd = 0;
                    for ( r=0; r<row; r++ )
                    {
                        if ( r == row_PD_ ) continue;

                        ImageType& IRImage = input_(r, c);

                        ImageType& magIRImage = magIR_(rInd, c);
                        magIRImage = IRImage;

                        Gadgetron::absolute(IRImage, magIRImage);

                        rInd++;
                    }
                }
            }

            if ( !debugFolder_.empty() )
            {
                hoNDArray<T> out;

                for ( r=0; r<row; r++ )
                {
                    input_.to_NDArray(r, out);
                    std::ostringstream ostr;
                    ostr << "Input_row" << r;
                    GADGET_EXPORT_ARRAY_COMPLEX(debugFolder_, gt_exporter_, out, ostr.str());
                }

                if ( perform_PSIR_ )
                {
                    for ( r=0; r<row-1; r++ )
                    {
                        PSIR_.to_NDArray(r, out);
                        std::ostringstream ostr;
                        ostr << "PSIR_row" << r;
                        GADGET_EXPORT_ARRAY_COMPLEX(debugFolder_, gt_exporter_, out, ostr.str());
                    }
                }

                for ( r=0; r<row-1; r++ )
                {
                    magIR_.to_NDArray(r, out);
                    std::ostringstream ostr;
                    ostr << "MagIR_row" << r;
                    GADGET_EXPORT_ARRAY_COMPLEX(debugFolder_, gt_exporter_, out, ostr.str());
                }

                magPD_.to_NDArray(0, out);
                std::ostringstream ostr;
                ostr << "magPD";
                GADGET_EXPORT_ARRAY_COMPLEX(debugFolder_, gt_exporter_, out, ostr.str());
            }

            if ( perform_SCC_ )
            {
                ImageType& magPD = magPD_(0, col_compute_PSIR_windowing_);

                // filter the magPD using median
                ImageMagType magPDReal(magPD.get_dimensions());
                Gadgetron::complex_to_real(magPD, magPDReal);

                ImageMagType magPDFiltered(magPDReal);
                Gadgetron::filterMedian(magPDReal, filter_width_, magPDFiltered);

                ImageType magPDFilteredCx(magPDFiltered.get_dimensions());
                Gadgetron::real_to_complex(magPDFiltered, magPDFilteredCx);

                ImageType PSIRImageBeforeSCC = PSIR_(0, col_compute_PSIR_windowing_);

                for ( c=0; c<cols[0]; c++ )
                {
                    for ( r=0; r<row-1; r++ )
                    {
                        if ( perform_PSIR_ )
                        {
                            ImageType& PSIRImage = PSIR_(r, c);
                            Gadgetron::divide(PSIRImage, magPDFilteredCx, PSIRImage);
                            Gadgetron::scal(scale_factor_after_SCC_, PSIRImage);
                        }

                        ImageType& magIRImage = magIR_(r, c);
                        Gadgetron::divide(magIRImage, magPDFilteredCx, magIRImage);
                        Gadgetron::scal(scale_factor_after_SCC_, magIRImage);
                    }
                }

                // compute the windowing
                if ( perform_PSIR_ && compute_PSIR_windowing_ )
                {
                    ImageType& PSIRImage = PSIR_(0, col_compute_PSIR_windowing_);
                    ImageType& magIRImage = magIR_(0, col_compute_PSIR_windowing_);

                    this->calculate_window_level(magPDFilteredCx, PSIRImage, PSIRImageBeforeSCC, window_center_, window_width_);
                }

                if ( !debugFolder_.empty() )
                {
                    hoNDArray<T> out;

                    if ( perform_PSIR_ )
                    {
                        for ( r=0; r<row-1; r++ )
                        {
                            PSIR_.to_NDArray(r, out);
                            std::ostringstream ostr;
                            ostr << "PSIR_Norm_row" << r;
                            GADGET_EXPORT_ARRAY_COMPLEX(debugFolder_, gt_exporter_, out, ostr.str());
                        }
                    }

                    for ( r=0; r<row-1; r++ )
                    {
                        magIR_.to_NDArray(r, out);
                        std::ostringstream ostr;
                        ostr << "MagIR_Norm_row" << r;
                        GADGET_EXPORT_ARRAY_COMPLEX(debugFolder_, gt_exporter_, out, ostr.str());
                    }

                    std::ostringstream ostr;
                    ostr << "magPD_used_for_Norm";
                    GADGET_EXPORT_IMAGE_COMPLEX(debugFolder_, gt_exporter_, magPD, ostr.str());

                    std::ostringstream ostr2;
                    ostr2 << "magPD_used_for_Norm_filtered";
                    GADGET_EXPORT_IMAGE_COMPLEX(debugFolder_, gt_exporter_, magPDFilteredCx, ostr.str());
                }
            }
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Error happened in gtPlusPhaseSensitiveRecon<ValueType, D>::performPSIRRecon(...) ... ");
            return false;
        }

        return true;
    }

    template<typename ValueType, unsigned int D> 
    bool gtPlusPhaseSensitiveRecon<ValueType, D>::calculate_window_level(ImageType& magPDFiltered, ImageType& PSIRImage, ImageType& PSIRImageBeforeSCC, float& window_center, float& window_width)
    {
        try
        {
            if ( !PSIRImage.attrib_.attribute2_.get(GTPLUS_IMAGE_SCALE_RATIO, 0, intensity_scale_factor_) )
            {
                intensity_scale_factor_ = 8;
            }

            real_value_type thres = intensity_scale_factor_ * 2;

            size_t N = PSIRImage.get_number_of_elements();

            std::vector<size_t> dim = PSIRImage.get_dimensions();
            ImageMagType mask(dim);
            Gadgetron::clear(mask);

            long long n;
            size_t numOfPixelInMask = 0;
            for ( n=0; n<N; n++ )
            {
                if ( magPDFiltered(n).real() > thres )
                {
                    mask(n) = 1;
                    numOfPixelInMask++;
                }
            }

            if ( numOfPixelInMask == 0 ) return true;

            std::vector<real_value_type> valueInMask(numOfPixelInMask, 0);
            size_t ind;
            for ( n=0; n<N; n++ )
            {
                if ( mask(n) == 1 )
                {
                    valueInMask[ind++] = PSIRImage(n).real();
                }
            }

            std::sort(valueInMask.begin(), valueInMask.end());
            real_value_type medianValueInMask = valueInMask[numOfPixelInMask/2];

            ImageMagType thresdhold(dim);
            Gadgetron::clear(thresdhold);

            for ( n=0; n<N; n++ )
            {
                if ( (PSIRImage(n).real()<medianValueInMask) && (mask(n) == 1) )
                {
                    thresdhold(n) = 1;
                }
            }

            real_value_type normal_level = 0;

            ind = 0;
            for ( n=0; n<N; n++ )
            {
                if ( thresdhold(n) == 1 )
                {
                    normal_level += PSIRImage(n).real();
                    ind++;
                }
            }
            normal_level /= ind;

            // make a 1D histogram
            size_t numOfBins = 64;
            hoNDArray<double> hist1D(numOfBins), cdf1D(numOfBins);
            hoNDArray<double> histBins(numOfBins);

            real_value_type minValue = valueInMask[0];
            real_value_type maxValue = valueInMask[numOfPixelInMask-1];

            real_value_type range_t = real_value_type(1.0)/(maxValue - minValue + std::numeric_limits<real_value_type>::epsilon());

            for ( n=0; n<numOfPixelInMask; n++ )
            {
                size_t indT = static_cast<size_t>( range_t*(valueInMask[n]-minValue)*(numOfBins-1) + 0.5 );
                hist1D(indT)++;
            }

            for ( n=0; n<numOfBins; n++ )
            {
                hist1D(n) /= numOfPixelInMask;
            }

            for ( n=0; n<numOfBins; n++ )
            {
                cdf1D(n) = 0;

                size_t m;
                for ( m=0; m<n; m++ )
                {
                    cdf1D(n) += hist1D(n);
                }
            }

            for ( n=0; n<numOfBins; n++ )
            {
                if ( cdf1D(n) > 0.97 ) break;
            }

            real_value_type f1dnormreal_val97 = n * (maxValue - minValue)/numOfBins;

            real_value_type range = 1.1 * f1dnormreal_val97 - normal_level;
            real_value_type min_level = f1dnormreal_val97 - range;

            window_center = min_level + range/2;
            window_width = range;
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Error happened in gtPlusPhaseSensitiveRecon<ValueType, D>::calculate_window_level(ImageType& magPDFiltered, ImageType& PSIRImage, ImageType& PSIRImageBeforeSCC, float& window_center, float& window_width) ... ");
            return false;
        }

        return true;
    }

    template<typename ValueType, unsigned int D> 
    void gtPlusPhaseSensitiveRecon<ValueType, D>::print(std::ostream& os) const
    {
        using namespace std;

        os << "-------------- GTPlus PSIR Reconstruction -------------" << endl;

        os << "Image dimension is : " << D << endl;

        std::string elemTypeName = std::string(typeid(ValueType).name());
        os << "Image data type is : " << elemTypeName << std::endl;

        os << "Row as the PD images is : " << row_PD_ << std::endl;
        os << "Column used for PSIR windowing computation is : " << col_compute_PSIR_windowing_ << std::endl;
        os << "-------------------------------------------------------" << endl;
    }

}}
