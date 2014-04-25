/** \file   gtPlusMoCoAveraging.h
    \brief  Implement 2D and 3D moco averaging
            The input is a 2D image container
            for every row, the moco+ave will be performed
            the different row can have the same or different reference frame
    \author Hui Xue
*/

#pragma once

#include "gtPlusISMRMRDReconUtil.h"
#include "hoImageRegContainer2DRegistration.h"

namespace Gadgetron { namespace gtPlus {

    struct compObj
    {
        compObj() {}
        ~compObj() {}

        bool operator()(const std::pair<double, size_t>& m1, const std::pair<double, size_t>& m2) const
        {
            return (m1.first <= m2.first);
        }
    };

    template<typename ValueType, typename CoordType, unsigned int D> 
    class gtPlusMoCoAveraging
    {
    public:

        typedef gtPlusMoCoAveraging<ValueType, CoordType, D> Self;

        typedef Gadgetron::hoNDImage<ValueType, D> ImageType;
        typedef Gadgetron::hoNDImage<ValueType, 2> Image2DType;
        typedef Gadgetron::hoNDImage<ValueType, 3> Image3DType;

        typedef ValueType T;
        typedef ValueType element_type;
        typedef ValueType value_type;

        typedef CoordType coord_type;

        typedef Gadgetron::hoImageRegContainer2DRegistration<ValueType, CoordType, D, D> RegContainer2DType;
        typedef typename RegContainer2DType::TransformationDeformationFieldType TransformationDeformationFieldType;
        typedef typename RegContainer2DType::DeformationFieldType DeformationFieldType;

        typedef Gadgetron::hoNDImageContainer2D<ImageType> ImageContinerType;

        typedef Gadgetron::hoNDImageContainer2D<DeformationFieldType> DeformationFieldContinerType;

        typedef std::vector< std::vector< std::pair<double, size_t> > > MoCoQualityType;

        gtPlusMoCoAveraging();
        virtual ~gtPlusMoCoAveraging();

        /// set the default parameters
        virtual bool setDefaultParameters();

        /// compute moco and moco quality, but not averaging
        virtual bool computeMoCoAveraging(const ImageContinerType& input)
        {
            GADGET_CHECK_RETURN_FALSE(this->initialize(input));
            GADGET_CHECK_RETURN_FALSE(this->pickReference());
            GADGET_CHECK_RETURN_FALSE(this->performMoCo());
            GADGET_CHECK_RETURN_FALSE(this->computeMoCoQuality());
            GADGET_CHECK_RETURN_FALSE(this->performAveraging(register_.warped_container_, averaged_));

            return true;
        }

        /// apply computed moco averaging on an input container and compute averaged images
        template<typename ValueType2> 
        bool applyMoCoAveraging(const Gadgetron::hoNDImageContainer2D< Gadgetron::hoNDImage<ValueType2, D> >& input, 
                Gadgetron::hoNDImageContainer2D< Gadgetron::hoNDImage<ValueType2, D> >& warpped, 
                Gadgetron::hoNDImageContainer2D< Gadgetron::hoNDImage<ValueType2, D> >& averaged)
        {
            GADGET_CHECK_RETURN_FALSE(register_.warpContainer2D(input, input, register_.deformation_field_, warpped));
            GADGET_CHECK_RETURN_FALSE(this->performAveraging(warpped, averaged));
            return true;
        }

        /// print the class information
        virtual void print(std::ostream& os) const;

        /// perform averaging
        template<typename ValueType2> 
        bool performAveraging( const Gadgetron::hoNDImageContainer2D< Gadgetron::hoNDImage<ValueType2, D> >& input, 
                Gadgetron::hoNDImageContainer2D< Gadgetron::hoNDImage<ValueType2, D> >& averaged )
        {
            try
            {
                size_t row = input.rows();
                std::vector<size_t> cols = input.cols();

                std::vector<size_t> dim(1);
                dim[0] = row;
                GADGET_CHECK_RETURN_FALSE(averaged.create(dim));

                size_t r;
                for ( r=0; r<row; r++ )
                {
                    long long col = (long long)cols[r];
                    long long usedImageForAve = col * percentage_kept_for_averaging_;

                    size_t ref = reference_[r];

                    averaged(0, r) = input(r, ref);
                    Gadgetron::clear(averaged(0, r));

                    if ( usedImageForAve <= 1 ) continue;

                    long long ii;
                    for ( ii=0; ii<usedImageForAve; ii++ )
                    {
                        Gadgetron::add(averaged(0, r), input(r, moco_quality_measure_[r][ii].second ), averaged(0, r));
                    }

                    Gadgetron::scal( (typename realType<ValueType2>::Type)( 1.0/std::sqrt( (double)usedImageForAve ) ),  averaged(0, r));
                }
            }
            catch(...)
            {
                GADGET_ERROR_MSG("Error happened in gtPlusMoCoAveraging<ValueType, CoordType, D>::performAveraging(...) ... ");
                return false;
            }

            return true;
        }

        /// perform the cross-series moco on the averaged_ container
        virtual bool performCrossRowMoCo(size_t refRow);

        template<typename ValueType2> 
        bool applyCrossRowMoCo(const Gadgetron::hoNDImageContainer2D< Gadgetron::hoNDImage<ValueType2, D> >& input, 
                Gadgetron::hoNDImageContainer2D< Gadgetron::hoNDImage<ValueType2, D> >& warpped)
        {
            GADGET_CHECK_RETURN_FALSE(register_cross_row_.warpContainer2D(input, input, register_cross_row_.deformation_field_, warpped));
            return true;
        }

        // ----------------------------------
        // parameters
        // ----------------------------------

        /// input image container containing images for processing
        ImageContinerType input_;

        /// the image container register
        RegContainer2DType register_;

        /// register for cross-row moco
        RegContainer2DType register_cross_row_;

        /// whether to perform cross row registration
        bool cross_row_reg_;

        /// all rows have the same reference
        bool cross_row_reference_;

        /// method to compute moco quality
        /// "SSD", "Deformation"
        std::string moco_quality_strategy_;

        /// averaging percentage
        /// e.g., 0.5 means half of the data is used for averaging
        double percentage_kept_for_averaging_;

        /// whether to perform soft averaging
        bool soft_averaging_;

        /// whether need the mocoed images from input
        /// sometimes, only the deformation filed is needed
        bool moco_input_needed_;

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

        // ----------------------------------
        // registration results
        // ----------------------------------

        /// mocoed images and deformation fields are stored in the register_

        /// averaged images, [numOfRow 1]
        ImageContinerType averaged_;

        /// refrence picked for every row
        std::vector<unsigned int> reference_;

        /// moco quality measurements for every row
        /// the moco quality of a frame is measured by a double value
        MoCoQualityType moco_quality_measure_;

    protected:

        bool initialize(const ImageContinerType& input);

        /// pick the reference
        virtual bool pickReference();

        /// perform moco averaging
        virtual bool performMoCo();

        /// compute moco quality measure
        bool computeMoCoQuality();

        /// compute moco quality measure using SSD
        bool computeMoCoQualitySSD(const ImageContinerType& warpped, const std::vector<unsigned int>& ref, MoCoQualityType& mocoQ);

        /// compute moco quality measure using deformation field
        bool computeMoCoQualityDeformation(const DeformationFieldContinerType* deform, const std::vector<unsigned int>& ref, MoCoQualityType& mocoQ);
    };

    template<typename ValueType, typename CoordType, unsigned int D> 
    gtPlusMoCoAveraging<ValueType, CoordType, D>::gtPlusMoCoAveraging() : performTiming_(false)
    {
        gt_timer1_.set_timing_in_destruction(false);
        gt_timer2_.set_timing_in_destruction(false);
        gt_timer3_.set_timing_in_destruction(false);

        GADGET_CHECK_THROW(this->setDefaultParameters());
    }

    template<typename ValueType, typename CoordType, unsigned int D> 
    gtPlusMoCoAveraging<ValueType, CoordType, D>::
    ~gtPlusMoCoAveraging()
    {
    }

    template<typename ValueType, typename CoordType, unsigned int D> 
    bool gtPlusMoCoAveraging<ValueType, CoordType, D>::setDefaultParameters()
    {
        GADGET_CHECK_RETURN_FALSE(register_.setDefaultParameters());

        cross_row_reg_ = false;
        cross_row_reference_ = false;
        moco_quality_strategy_ = "SSD";
        percentage_kept_for_averaging_ = 0.5;
        soft_averaging_ = true;;
        moco_input_needed_ = true;

        verbose_ = false;

        return true;
    }

    template<typename ValueType, typename CoordType, unsigned int D> 
    bool gtPlusMoCoAveraging<ValueType, CoordType, D>::
    initialize(const ImageContinerType& input)
    {
        try
        {
            size_t row = input.rows();
            std::vector<size_t> dim(1, row);

            // allocate averaged images
            GADGET_CHECK_RETURN_FALSE(averaged_.create(dim));

            reference_.resize(row, 0);

            moco_quality_measure_.resize(row);

            std::vector<size_t> cols = input.cols();

            size_t r, c;
            for ( r=0; r<row; r++ )
            {
                moco_quality_measure_[r].resize(cols[r]);
                for ( c=0; c<cols[r]; c++ )
                {
                    moco_quality_measure_[r][c].first = 0;
                    moco_quality_measure_[r][c].second = c;
                }
            }

            input_ = input;
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Error happened in gtPlusMoCoAveraging<ValueType, CoordType, D>::initialize(const TargetContinerType& targetContainer) ... ");
            return false;
        }

        return true;
    }

    template<typename ValueType, typename CoordType, unsigned int D> 
    bool gtPlusMoCoAveraging<ValueType, CoordType, D>::pickReference()
    {
        try
        {
            size_t row = input_.rows();
            std::vector<size_t> cols = input_.cols();

            size_t r;
            for ( r=0; r<row; r++ )
            {
                long long col = (long long)cols[r];

                hoMatrixReal<ValueType> SSD(col, col);
                Gadgetron::clear(SSD);

                long long m, n;
                #pragma omp parallel default(none) private(m, n) shared(r, col, SSD)
                {
                    ImageType diff(input_(r, 0).get_dimensions());

                    #pragma omp for
                    for ( m=0; m<col; m++ )
                    {
                        ValueType v(0);
                        for ( n=m+1; n<col; n++ )
                        {
                            Gadgetron::subtract(input_(r, m), input_(r, n), diff);

                            Gadgetron::norm2(diff, v);
                            SSD(m, n) = v;
                            SSD(n, m) = v;
                        }
                    }

                    diff.clear();
                }

                // sort for every column
                GADGET_CHECK_RETURN_FALSE(SSD.sort_ascending_along_row());

                // pick the middel row
                hoMatrixReal<ValueType> minimalSSDCol;
                GADGET_CHECK_RETURN_FALSE(SSD.subMatrix(minimalSSDCol, col/2, col/2, 0, col-1));

                // find minimal SSD
                ValueType minSSD;
                size_t ind;
                Gadgetron::minAbsolute(minimalSSDCol, minSSD, ind);

                reference_[r] = ind;
            }

            if ( cross_row_reference_ && row>1 )
            {
                std::vector<unsigned int> refComputed(reference_);

                std::sort(refComputed.begin(), refComputed.end());
                size_t commonRef = refComputed[row/2];

                for ( r=0; r<row; r++ )
                {
                    if ( commonRef < cols[r] )
                    {
                        reference_[r] = commonRef;
                    }
                }
            }

            if ( verbose_ )
            {
                GADGET_MSG("gtPlusMoCoAveraging - reference : ");
                for ( r=0; r<row; r++ )
                {
                    GADGET_MSG("row " << r << " - " << reference_[r]);
                }
            }
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Error happened in gtPlusMoCoAveraging<ValueType, CoordType, D>::pickReference(...) ... ");
            return false;
        }

        return true;
    }

    template<typename ValueType, typename CoordType, unsigned int D> 
    bool gtPlusMoCoAveraging<ValueType, CoordType, D>::performMoCo()
    {
        try
        {
            register_.debugFolder_ = debugFolder_;
            register_.verbose_ = verbose_;

            if ( verbose_ )
            {
                register_.print(std::cout);
            }

            bool warped = moco_input_needed_;

            GADGET_CHECK_PERFORM(performTiming_, gt_timer3_.start("gtPlusMoCoAveraging - perform moco ... "));

            if ( register_.container_reg_mode_ == GT_IMAGE_REG_CONTAINER_FIXED_REFERENCE )
            {
                register_.registerOverContainer2DFixedReference(input_, reference_, warped);
            }

            if ( register_.container_reg_mode_ == GT_IMAGE_REG_CONTAINER_PROGRESSIVE )
            {
                warped = true;
                register_.registerOverContainer2DProgressive(input_, reference_);
            }

            GADGET_CHECK_PERFORM(performTiming_, gt_timer3_.stop());

            if ( !debugFolder_.empty() && warped )
            {
                hoNDArray<T> out;

                size_t row = input_.rows();

                size_t r;
                for ( r=0; r<row; r++ )
                {
                    register_.warped_container_.to_NDArray(r, out);

                    std::ostringstream ostr;
                    ostr << "MOCO_row" << r;

                    GADGET_EXPORT_ARRAY(debugFolder_, gt_exporter_, out, ostr.str());
                }
            }
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Error happened in gtPlusMoCoAveraging<ValueType, CoordType, D>::performMoCo(...) ... ");
            return false;
        }

        return true;
    }

    template<typename ValueType, typename CoordType, unsigned int D> 
    bool gtPlusMoCoAveraging<ValueType, CoordType, D>::performCrossRowMoCo(size_t refRow)
    {
        try
        {
            size_t col = averaged_.cols(0);
            GADGET_CHECK_RETURN_FALSE(refRow < col);

            register_cross_row_.debugFolder_ = debugFolder_;
            register_cross_row_.verbose_ = verbose_;

            if ( verbose_ )
            {
                register_cross_row_.print(std::cout);
            }

            bool warped = moco_input_needed_;

            GADGET_CHECK_PERFORM(performTiming_, gt_timer3_.start("gtPlusMoCoAveraging - perform cross-row moco ... "));

            std::vector<unsigned int> refCrossRow(1, refRow);

            if ( register_cross_row_.container_reg_mode_ == GT_IMAGE_REG_CONTAINER_FIXED_REFERENCE )
            {
                register_cross_row_.registerOverContainer2DFixedReference(averaged_, refCrossRow, warped);
            }

            if ( register_.container_reg_mode_ == GT_IMAGE_REG_CONTAINER_PROGRESSIVE )
            {
                warped = true;
                register_cross_row_.registerOverContainer2DProgressive(averaged_, refCrossRow);
            }

            GADGET_CHECK_PERFORM(performTiming_, gt_timer3_.stop());

            if ( !debugFolder_.empty() && warped )
            {
                hoNDArray<T> out;
                register_cross_row_.warped_container_.to_NDArray(0, out);

                std::ostringstream ostr;
                ostr << "MOCO_cross_row";

                GADGET_EXPORT_ARRAY(debugFolder_, gt_exporter_, out, ostr.str());
            }
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Error happened in gtPlusMoCoAveraging<ValueType, CoordType, D>::performCrossRowMoCo(size_t refRow) ... ");
            return false;
        }

        return true;
    }

    template<typename ValueType, typename CoordType, unsigned int D> 
    bool gtPlusMoCoAveraging<ValueType, CoordType, D>::computeMoCoQualitySSD(const ImageContinerType& warpped, const std::vector<unsigned int>& reference, MoCoQualityType& mocoQ)
    {
        try
        {
            size_t row = warpped.rows();
            std::vector<size_t> cols = warpped.cols();

            long long r;

            for ( r=0; r<(long long)row; r++ )
            {
                size_t ref = reference[r];
                mocoQ[r][ref].first = 0;

                const ImageType& refImage = warpped(r, ref);

                long long c;
                {
                    #pragma omp parallel default(none) private(c) shared(refImage, cols, r, warpped, mocoQ, ref) if (D==2)
                    {
                        ImageType diff(refImage);
                        ValueType v;

                        #pragma omp for 
                        for ( c=0; c<(long long)cols[r]; c++ )
                        {
                            if ( c == ref ) continue;

                            const ImageType& warppedImage = warpped(r, c);

                            Gadgetron::subtract(refImage, warppedImage, diff);
                            Gadgetron::norm2(diff, v);
                            mocoQ[r][c].first = v;
                        }

                        diff.clear();
                    }
                }

                std::sort(mocoQ[r].begin(), mocoQ[r].end(), compObj() );
            }
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Error happened in gtPlusMoCoAveraging<ValueType, CoordType, D>::computeMoCoQualitySSD(...) ... ");
            return false;
        }

        return true;
    }

    template<typename ValueType, typename CoordType, unsigned int D> 
    bool gtPlusMoCoAveraging<ValueType, CoordType, D>::computeMoCoQualityDeformation(const DeformationFieldContinerType* deform, const std::vector<unsigned int>& reference, MoCoQualityType& mocoQ)
    {
        try
        {
            size_t row = deform[0].rows();
            std::vector<size_t> cols = deform[0].cols();

            long long r;

            for ( r=0; r<(long long)row; r++ )
            {
                size_t ref = reference[r];
                mocoQ[r][ref].first = 0;

                long long c;
                {
                    #pragma omp parallel default(none) private(c) shared(cols, r, deform, mocoQ, ref) if (D==2)
                    {
                        TransformationDeformationFieldType deformTransform;
                        DeformationFieldType* deformField[D];

                        hoNDArray<ValueType> jac2D;
                        ValueType meanDeform, maxDeform, meanLogJac2D, maxLogJac2D;

                        #pragma omp for 
                        for ( c=0; c<(long long)cols[r]; c++ )
                        {
                            if ( c == ref ) continue;

                            unsigned int ii;
                            for ( ii=0; ii<D; ii++ )
                            {
                                deformField[ii] = const_cast<DeformationFieldType*>( &(deform[ii](r, c)) );
                            }

                            deformTransform.jacobianPosition(jac2D, deformField, 1);
                            deformTransform.analyzeJacobianAndDeformation(jac2D, deformField, meanDeform, maxDeform, meanLogJac2D, maxLogJac2D, 1);

                            mocoQ[r][ref].first = meanDeform;
                        }
                    }
                }

                std::sort(mocoQ[r].begin(), mocoQ[r].end(), compObj() );
            }
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Error happened in gtPlusMoCoAveraging<ValueType, CoordType, D>::computeMoCoQualitySSD(...) ... ");
            return false;
        }

        return true;
    }

    template<typename ValueType, typename CoordType, unsigned int D> 
    bool gtPlusMoCoAveraging<ValueType, CoordType, D>::computeMoCoQuality()
    {
        try
        {
            if ( moco_quality_strategy_ == "SSD" )
            {
                GADGET_CHECK_RETURN_FALSE(computeMoCoQualitySSD(register_.warped_container_, reference_, moco_quality_measure_));
            }
            else if ( moco_quality_strategy_ == "Deformation" )
            {
                GADGET_CHECK_RETURN_FALSE(computeMoCoQualityDeformation(register_.deformation_field_, reference_, moco_quality_measure_));
            }
            else
            {
                GADGET_ERROR_MSG("Incorrect moco quality measurement mode : " << moco_quality_strategy_);
                return false;
            }
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Error happened in gtPlusMoCoAveraging<ValueType, CoordType, D>::computeMoCoQuality(...) ... ");
            return false;
        }

        return true;
    }

    template<typename ValueType, typename CoordType, unsigned int D> 
    void gtPlusMoCoAveraging<ValueType, CoordType, D>::print(std::ostream& os) const
    {
        using namespace std;

        os << "-------------- GTPlus MoCo Averaging -------------" << endl;

        os << "Image dimension is : " << D << endl;

        std::string elemTypeName = std::string(typeid(ValueType).name());
        os << "Image data type is : " << elemTypeName << std::endl;

        elemTypeName = std::string(typeid(CoordType).name());
        os << "Transformation coordinate data type is : " << elemTypeName << std::endl;

        os << "Whether to perform cross row registration is : " << cross_row_reg_ << std::endl;
        os << "Method to compute moco quality is : " << moco_quality_strategy_ << std::endl;
        os << "Averaging percentage is : " << percentage_kept_for_averaging_ << std::endl;
        os << "Whether to perform soft averaging is : " << soft_averaging_ << std::endl;
        os << "--------------------------------------------------" << endl;
        os << "Info of register is : " << endl;
        os << "--------------------------------------------------" << endl;
        register_.print(os);
        os << "--------------------------------------------------" << endl;
    }

}}
