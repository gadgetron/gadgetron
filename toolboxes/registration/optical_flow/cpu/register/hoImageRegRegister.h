/** \file   hoImageRegRegister.h
    \brief  Define the class to perform image registration in gadgetron
    \author Hui Xue
*/

#ifndef hoImageRegRegister_H_
#define hoImageRegRegister_H_

#pragma once

#include "hoNDArray.h"
#include "hoNDImage.h"
#include "hoNDInterpolator.h"
#include "hoNDBoundaryHandler.h"
#include "hoMatrix.h"
#include "hoNDArray_utils.h"
#include "hoNDArray_elemwise.h"
#include "hoNDImage_util.h"

// transformation
#include "hoImageRegTransformation.h"
#include "hoImageRegParametricTransformation.h"
#include "hoImageRegTransformation.h"
#include "hoImageRegHomogenousTransformation.h"
#include "hoImageRegRigid2DTransformation.h"
#include "hoImageRegRigid3DTransformation.h"

// warper
#include "hoImageRegWarper.h"

// solver
#include "hoImageRegDeformationFieldSolver.h"
#include "hoImageRegParametricSolver.h"
#include "hoImageRegDeformationFieldBidirectionalSolver.h"
#include "hoImageRegParametricDownHillSolver.h"
#include "hoImageRegParametricGradientDescentSolver.h"

// dissimilarity
#include "hoImageRegDissimilaritySSD.h"
#include "hoImageRegDissimilarityLocalCCR.h"
#include "hoImageRegDissimilarityMutualInformation.h"
#include "hoImageRegDissimilarityNormalizedMutualInformation.h"

namespace Gadgetron {

    /// perform the image registration using pyramid scheme
    template<typename TargetType, typename SourceType, typename CoordType> 
    class hoImageRegRegister
    {
    public:

        typedef hoImageRegRegister<TargetType, SourceType, CoordType> Self;

        typedef typename TargetType::value_type ValueType;
        enum { DIn = TargetType::NDIM };
        enum { DOut = SourceType::NDIM };

        typedef hoNDImage<ValueType, 2> Target2DType;
        typedef Target2DType Source2DType;

        typedef hoNDImage<ValueType, 3> Target3DType;
        typedef Target2DType Source3DType;

        typedef ValueType T;
        typedef ValueType element_type;
        typedef ValueType value_type;

        typedef CoordType coord_type;

        /// boundary handler and interpolator for target image
        typedef hoNDBoundaryHandler<TargetType> BoundaryHandlerTargetType;
        typedef hoNDBoundaryHandlerFixedValue<TargetType> BoundaryHandlerTargetFixedValueType;
        typedef hoNDBoundaryHandlerBorderValue<TargetType> BoundaryHandlerTargetBorderValueType;
        typedef hoNDBoundaryHandlerPeriodic<TargetType> BoundaryHandlerTargetPeriodicType;
        typedef hoNDBoundaryHandlerMirror<TargetType> BoundaryHandlerTargetMirrorType;

        typedef hoNDInterpolator<TargetType> InterpTargetType;
        typedef hoNDInterpolatorLinear<TargetType> InterpTargetLinearType;
        typedef hoNDInterpolatorNearestNeighbor<TargetType> InterpTargetNearestNeighborType;
        typedef hoNDInterpolatorBSpline<TargetType, DIn> InterpTargetBSplineType;

        /// boundary handler and interpolator for source image
        typedef hoNDBoundaryHandler<SourceType> BoundaryHandlerSourceType;
        typedef hoNDBoundaryHandlerFixedValue<SourceType> BoundaryHandlerSourceFixedValueType;
        typedef hoNDBoundaryHandlerBorderValue<SourceType> BoundaryHandlerSourceBorderValueType;
        typedef hoNDBoundaryHandlerPeriodic<SourceType> BoundaryHandlerSourcePeriodicType;
        typedef hoNDBoundaryHandlerMirror<SourceType> BoundaryHandlerSourceMirrorType;

        typedef hoNDInterpolator<SourceType> InterpSourceType;
        typedef hoNDInterpolatorLinear<SourceType> InterpSourceLinearType;
        typedef hoNDInterpolatorNearestNeighbor<SourceType> InterpSourceNearestNeighborType;
        typedef hoNDInterpolatorBSpline<SourceType, DIn> InterpSourceBSplineType;

        /// warper type
        typedef hoImageRegWarper<TargetType, SourceType, CoordType> WarperType;

        /// image dissimilarity type
        typedef hoImageRegDissimilarity<SourceType> DissimilarityType;

        hoImageRegRegister(unsigned int resolution_pyramid_levels=3, ValueType bg_value=ValueType(0));
        virtual ~hoImageRegRegister();

        /// initialize the registration
        /// should be called after all images and parameters of registration are set
        virtual bool initialize();

        /// set target and source, create the multi-resolution pyramid and set up the interpolators
        virtual void setTarget(TargetType& target);
        virtual void setSource(SourceType& source);

        /// create dissimilarity measures
        DissimilarityType* createDissimilarity(GT_IMAGE_DISSIMILARITY v, unsigned int level);

        /// perform the registration
        virtual bool performRegistration() = 0;

        /// print the class information
        virtual void printContent(std::ostream& os) const;
        virtual void print(std::ostream& os) const;

        /// parameters

        /// whether to perform world coordinate registration
        bool use_world_coordinates_;

        /// number of resolution pyramid levels
        unsigned int resolution_pyramid_levels_;

        /// use fast pyramid creation by dividing the image size by 2
        /// if the use_world_coordinates_ == true and resolution_pyramid_divided_by_2_ == true, , resolution_pyramid_downsample_ratio_
        /// and resolution_pyramid_blurring_sigma_ will be ignored
        bool resolution_pyramid_divided_by_2_;

        /// downsample ratio of the resolution pyramid for every dimension and every level
        /// e.g. ratio=2, downsample by 100%
        std::vector< std::vector<float> > resolution_pyramid_downsample_ratio_;

        /// extra gaussian blurring can be applied on every resolution pyramid
        /// if use_world_coordinates_=true, sigma is in the unit of world coordinate
        /// otherwise, it is in the unit of image pixel
        std::vector< std::vector<float> > resolution_pyramid_blurring_sigma_;

        /// boundary handler and interpolator type for warper, for every resolution level, different interpolator can be used
        std::vector<GT_BOUNDARY_CONDITION> boundary_handler_type_warper_;
        std::vector<GT_IMAGE_INTERPOLATOR> interp_type_warper_;

        /// boundary handler and interpolator type for pyramid construction
        GT_BOUNDARY_CONDITION boundary_handler_type_pyramid_construction_;
        GT_IMAGE_INTERPOLATOR interp_type_pyramid_construction_;

        /// image dissimilarity
        /// for different pyramid level, different dissimilarity can be used
        std::vector<GT_IMAGE_DISSIMILARITY> dissimilarity_type_;

        /// solver for every pyramid level
        std::vector<GT_IMAGE_REG_SOLVER> solver_type_;

        ///// whether to set the origin of target/source to image center
        //bool orgin_at_image_center_;

        /// parameters for dissimilarity measures, for every paramid level
        /// LocalCCR
        std::vector<std::vector<ValueType> > dissimilarity_LocalCCR_sigmaArg_;

        /// Histogram based
        std::vector<unsigned int> dissimilarity_hist_num_bin_target_;
        std::vector<unsigned int> dissimilarity_hist_num_bin_warpped_;
        bool dissimilarity_hist_pv_interpolation_;
        std::vector<size_t> dissimilarity_hist_step_size_ignore_pixel_;

        /// Mutual information
        std::vector<ValueType> dissimilarity_MI_betaArg_;

        // ----------------------------------
        // debug and timing
        // ----------------------------------
        // clock for timing
        Gadgetron::GadgetronTimer gt_timer1_;
        Gadgetron::GadgetronTimer gt_timer2_;
        Gadgetron::GadgetronTimer gt_timer3_;

        bool performTiming_;

        // exporter
        Gadgetron::ImageIOAnalyze gt_exporter_;

        // debug folder
        std::string debugFolder_;

    protected:

        TargetType* target_;
        SourceType* source_;

        /// back ground values, used to mark regions in the target image which will not be warped
        ValueType bg_value_;

        /// store the multi-resolution images for every pyramid level
        std::vector<TargetType> target_pyramid_;
        std::vector<TargetType> source_pyramid_;

        /// store the boundary handler and interpolator for warpers
        std::vector<BoundaryHandlerTargetType*> target_bh_warper_;
        std::vector<InterpTargetType*> target_interp_warper_;

        std::vector<BoundaryHandlerSourceType*> source_bh_warper_;
        std::vector<InterpSourceType*> source_interp_warper_;

        /// store the boundary handler and interpolator for pyramid construction
        BoundaryHandlerTargetType* target_bh_pyramid_construction_;
        InterpTargetType* target_interp_pyramid_construction_;

        BoundaryHandlerSourceType* source_bh_pyramid_construction_;
        InterpSourceType* source_interp_pyramid_construction_;

        /// store warpers for ever pyramid level
        std::vector<WarperType> warper_pyramid_;

        /// store the image dissimilarity for every pyramid level
        std::vector<DissimilarityType*> dissimilarity_pyramid_;

        /// store warpers for ever pyramid level
        std::vector<WarperType> warper_pyramid_inverse_;

        /// store the image dissimilarity for every pyramid level
        std::vector<DissimilarityType*> dissimilarity_pyramid_inverse_;
    };

    template<typename TargetType, typename SourceType, typename CoordType> 
    hoImageRegRegister<TargetType, SourceType, CoordType>::
    hoImageRegRegister(unsigned int resolution_pyramid_levels, ValueType bg_value) 
    : target_(NULL), source_(NULL), bg_value_(bg_value), performTiming_(false)
    {
        gt_timer1_.set_timing_in_destruction(false);
        gt_timer2_.set_timing_in_destruction(false);
        gt_timer3_.set_timing_in_destruction(false);

        use_world_coordinates_ = true;

        resolution_pyramid_levels_ = resolution_pyramid_levels;

        resolution_pyramid_divided_by_2_ = true;

        resolution_pyramid_downsample_ratio_.resize(resolution_pyramid_levels_-1, std::vector<float>(DIn, 2.0f) );

        resolution_pyramid_blurring_sigma_.resize(resolution_pyramid_levels_, std::vector<float>(DIn, 0.0f));

        boundary_handler_type_warper_.resize(resolution_pyramid_levels_, GT_BOUNDARY_CONDITION_FIXEDVALUE);
        interp_type_warper_.resize(resolution_pyramid_levels_, GT_IMAGE_INTERPOLATOR_LINEAR);

        boundary_handler_type_pyramid_construction_ = GT_BOUNDARY_CONDITION_BORDERVALUE;
        interp_type_pyramid_construction_ = GT_IMAGE_INTERPOLATOR_LINEAR;

        dissimilarity_type_.resize(resolution_pyramid_levels_, GT_IMAGE_DISSIMILARITY_NMI);

        solver_type_.resize(resolution_pyramid_levels_, GT_IMAGE_REG_SOLVER_DOWNHILL);

        target_bh_warper_.resize(resolution_pyramid_levels_, NULL);
        target_interp_warper_.resize(resolution_pyramid_levels_, NULL);

        source_bh_warper_.resize(resolution_pyramid_levels_, NULL);
        source_interp_warper_.resize(resolution_pyramid_levels_, NULL);

        target_bh_pyramid_construction_ = NULL;
        target_interp_pyramid_construction_ = NULL;

        source_bh_pyramid_construction_ = NULL;
        source_interp_pyramid_construction_ = NULL;

        dissimilarity_pyramid_.resize(resolution_pyramid_levels_, NULL);
        dissimilarity_pyramid_inverse_.resize(resolution_pyramid_levels_, NULL);

        dissimilarity_LocalCCR_sigmaArg_.resize(resolution_pyramid_levels_, std::vector<ValueType>(DOut, 2.0) );

        dissimilarity_hist_num_bin_target_.resize(resolution_pyramid_levels_, 64);
        dissimilarity_hist_num_bin_warpped_.resize(resolution_pyramid_levels_, 64);
        dissimilarity_hist_pv_interpolation_ = false;
        dissimilarity_hist_step_size_ignore_pixel_.resize(resolution_pyramid_levels_, 1);

        dissimilarity_MI_betaArg_.resize(resolution_pyramid_levels_, 2.0);
    }

    template<typename TargetType, typename SourceType, typename CoordType> 
    hoImageRegRegister<TargetType, SourceType, CoordType>::~hoImageRegRegister()
    {
        unsigned int ii;
        for ( ii=0; ii<resolution_pyramid_levels_; ii++ )
        {
            delete target_bh_warper_[ii];
            delete target_interp_warper_[ii];

            delete source_bh_warper_[ii];
            delete source_interp_warper_[ii];

            delete dissimilarity_pyramid_[ii];
            delete dissimilarity_pyramid_inverse_[ii];
        }

        delete target_bh_pyramid_construction_;
        delete target_interp_pyramid_construction_;

        delete source_bh_pyramid_construction_;
        delete source_interp_pyramid_construction_;
    }

    template<typename TargetType, typename SourceType, typename CoordType> 
    hoImageRegDissimilarity<SourceType>* hoImageRegRegister<TargetType, SourceType, CoordType>::createDissimilarity(GT_IMAGE_DISSIMILARITY v, unsigned int level)
    {
        hoImageRegDissimilarity<SourceType>* res = NULL;

        unsigned int ii;

        switch (v)
        {
            case GT_IMAGE_DISSIMILARITY_SSD:
                res = new hoImageRegDissimilaritySSD<SourceType>();
                break;

            case GT_IMAGE_DISSIMILARITY_LocalCCR:
            {
                hoImageRegDissimilarityLocalCCR<SourceType>* ptr = new hoImageRegDissimilarityLocalCCR<SourceType>();
                for ( ii=0; ii<DOut; ii++ )
                {
                    ptr->sigmaArg_[ii] = dissimilarity_LocalCCR_sigmaArg_[level][ii];
                }

                res = ptr;
            }
                break;

            case GT_IMAGE_DISSIMILARITY_MI:
            {
                hoImageRegDissimilarityMutualInformation<SourceType>* ptr = new hoImageRegDissimilarityMutualInformation<SourceType>();

                ptr->betaArg_[0] = dissimilarity_MI_betaArg_[level];
                ptr->betaArg_[1] = dissimilarity_MI_betaArg_[level];
                ptr->num_bin_target_ = dissimilarity_hist_num_bin_target_[level];
                ptr->num_bin_warpped_ = dissimilarity_hist_num_bin_warpped_[level];
                ptr->pv_interpolation_ = dissimilarity_hist_pv_interpolation_;
                ptr->step_size_ignore_pixel_ = dissimilarity_hist_step_size_ignore_pixel_[level];

                res = ptr;
            }
                break;

            case GT_IMAGE_DISSIMILARITY_NMI:
            {
                hoImageRegDissimilarityNormalizedMutualInformation<SourceType>* ptr = new hoImageRegDissimilarityNormalizedMutualInformation<SourceType>();

                ptr->num_bin_target_ = dissimilarity_hist_num_bin_target_[level];
                ptr->num_bin_warpped_ = dissimilarity_hist_num_bin_warpped_[level];
                ptr->pv_interpolation_ = dissimilarity_hist_pv_interpolation_;
                ptr->step_size_ignore_pixel_ = dissimilarity_hist_step_size_ignore_pixel_[level];

                res = ptr;
            }
                break;

            default:
                GERROR_STREAM("Unrecognized image dissimilarity type : " << v);
        }

        res->setBackgroundValue(bg_value_);

        return res;
    }

    template<typename TargetType, typename SourceType, typename CoordType> 
    bool hoImageRegRegister<TargetType, SourceType, CoordType>::initialize()
    {
        try
        {
            GADGET_CHECK_RETURN_FALSE(target_!=NULL);
            GADGET_CHECK_RETURN_FALSE(source_!=NULL);

            GADGET_CHECK_RETURN_FALSE(resolution_pyramid_downsample_ratio_.size()==resolution_pyramid_levels_-1);
            GADGET_CHECK_RETURN_FALSE(resolution_pyramid_blurring_sigma_.size()==resolution_pyramid_levels_);

            GADGET_CHECK_RETURN_FALSE(boundary_handler_type_warper_.size()==resolution_pyramid_levels_);
            GADGET_CHECK_RETURN_FALSE(interp_type_warper_.size()==resolution_pyramid_levels_);

            GADGET_CHECK_RETURN_FALSE(dissimilarity_type_.size()==resolution_pyramid_levels_);
            GADGET_CHECK_RETURN_FALSE(solver_type_.size()==resolution_pyramid_levels_);

            target_pyramid_.resize(resolution_pyramid_levels_);
            source_pyramid_.resize(resolution_pyramid_levels_);

            target_pyramid_[0] = *target_;
            source_pyramid_[0] = *source_;

            target_bh_pyramid_construction_ = createBoundaryHandler<TargetType>(boundary_handler_type_pyramid_construction_);
            target_interp_pyramid_construction_ = createInterpolator<TargetType, DOut>(interp_type_pyramid_construction_);
            target_interp_pyramid_construction_->setBoundaryHandler(*target_bh_pyramid_construction_);

            source_bh_pyramid_construction_ = createBoundaryHandler<SourceType>(boundary_handler_type_pyramid_construction_);
            source_interp_pyramid_construction_ = createInterpolator<SourceType, DIn>(interp_type_pyramid_construction_);
            source_interp_pyramid_construction_->setBoundaryHandler(*source_bh_pyramid_construction_);

            /// allocate all objects
            unsigned int ii, jj;
            for ( ii=0; ii<resolution_pyramid_levels_-1; ii++ )
            {
                // create pyramid
                target_bh_pyramid_construction_->setArray(target_pyramid_[ii]);
                target_interp_pyramid_construction_->setArray(target_pyramid_[ii]);

                if ( use_world_coordinates_ )
                {
                    if ( resolution_pyramid_divided_by_2_ )
                    {
                        Gadgetron::downsampleImageBy2WithAveraging(target_pyramid_[ii], *target_bh_pyramid_construction_, target_pyramid_[ii+1]);
                    }
                    else
                    {
                        std::vector<float> ratio = resolution_pyramid_downsample_ratio_[ii];
                        Gadgetron::downsampleImage(target_pyramid_[ii], *target_interp_pyramid_construction_, target_pyramid_[ii+1], &ratio[0]);

                        std::vector<float> sigma = resolution_pyramid_blurring_sigma_[ii+1];
                        for ( jj=0; jj<DOut; jj++ )
                        {
                            sigma[jj] /= target_pyramid_[ii+1].get_pixel_size(jj); // world to pixel
                        }

                        Gadgetron::filterGaussian(target_pyramid_[ii+1], &sigma[0]);
                    }
                }
                else
                {
                    std::vector<float> ratio = resolution_pyramid_downsample_ratio_[ii];

                    bool downsampledBy2 = true;
                    for ( jj=0; jj<DOut; jj++ )
                    {
                        if ( std::abs(ratio[jj]-2.0f) > FLT_EPSILON )
                        {
                            downsampledBy2 = false;
                            break;
                        }
                    }

                    if ( downsampledBy2 )
                    {
                        Gadgetron::downsampleImageBy2WithAveraging(target_pyramid_[ii], *target_bh_pyramid_construction_, target_pyramid_[ii+1]);
                        // Gadgetron::downsampleImage(target_pyramid_[ii], *target_interp_pyramid_construction_, target_pyramid_[ii+1], &ratio[0]);
                    }
                    else
                    {
                        Gadgetron::downsampleImage(target_pyramid_[ii], *target_interp_pyramid_construction_, target_pyramid_[ii+1], &ratio[0]);
                        std::vector<float> sigma = resolution_pyramid_blurring_sigma_[ii+1];
                        Gadgetron::filterGaussian(target_pyramid_[ii+1], &sigma[0]);
                    }
                }

                // source

                source_bh_pyramid_construction_->setArray(source_pyramid_[ii]);
                source_interp_pyramid_construction_->setArray(source_pyramid_[ii]);

                if ( use_world_coordinates_ )
                {
                    if ( resolution_pyramid_divided_by_2_ )
                    {
                        Gadgetron::downsampleImageBy2WithAveraging(source_pyramid_[ii], *source_bh_pyramid_construction_, source_pyramid_[ii+1]);
                    }
                    else
                    {
                        std::vector<float> ratio = resolution_pyramid_downsample_ratio_[ii];
                        Gadgetron::downsampleImage(source_pyramid_[ii], *source_interp_pyramid_construction_, source_pyramid_[ii+1], &ratio[0]);

                        std::vector<float> sigma = resolution_pyramid_blurring_sigma_[ii+1];
                        for ( jj=0; jj<DOut; jj++ )
                        {
                            sigma[jj] /= source_pyramid_[ii+1].get_pixel_size(jj); // world to pixel
                        }

                        Gadgetron::filterGaussian(source_pyramid_[ii+1], &sigma[0]);
                    }
                }
                else
                {
                    std::vector<float> ratio = resolution_pyramid_downsample_ratio_[ii];

                    bool downsampledBy2 = true;
                    for ( jj=0; jj<DOut; jj++ )
                    {
                        if ( std::abs(ratio[jj]-2.0f) > FLT_EPSILON )
                        {
                            downsampledBy2 = false;
                            break;
                        }
                    }

                    if ( downsampledBy2 )
                    {
                        Gadgetron::downsampleImageBy2WithAveraging(source_pyramid_[ii], *source_bh_pyramid_construction_, source_pyramid_[ii+1]);
                        //Gadgetron::downsampleImage(source_pyramid_[ii], *source_interp_pyramid_construction_, source_pyramid_[ii+1], &ratio[0]);
                    }
                    else
                    {
                        Gadgetron::downsampleImage(source_pyramid_[ii], *source_interp_pyramid_construction_, source_pyramid_[ii+1], &ratio[0]);
                        std::vector<float> sigma = resolution_pyramid_blurring_sigma_[ii+1];
                        Gadgetron::filterGaussian(source_pyramid_[ii+1], &sigma[0]);
                    }
                }
            }

            for ( ii=0; ii<resolution_pyramid_levels_; ii++ )
            {
                target_bh_warper_[ii] = createBoundaryHandler<TargetType>(boundary_handler_type_warper_[ii]);
                target_bh_warper_[ii]->setArray(target_pyramid_[ii]);

                target_interp_warper_[ii] = createInterpolator<TargetType, DOut>(interp_type_warper_[ii]);
                target_interp_warper_[ii]->setArray(target_pyramid_[ii]);
                target_interp_warper_[ii]->setBoundaryHandler(*target_bh_warper_[ii]);

                source_bh_warper_[ii] = createBoundaryHandler<SourceType>(boundary_handler_type_warper_[ii]);
                source_bh_warper_[ii]->setArray(source_pyramid_[ii]);

                source_interp_warper_[ii] = createInterpolator<SourceType, DIn>(interp_type_warper_[ii]);
                source_interp_warper_[ii]->setArray(source_pyramid_[ii]);
                source_interp_warper_[ii]->setBoundaryHandler(*source_bh_warper_[ii]);

                dissimilarity_pyramid_[ii] = createDissimilarity(dissimilarity_type_[ii], ii);
                dissimilarity_pyramid_[ii]->initialize(target_pyramid_[ii]);
                dissimilarity_pyramid_[ii]->debugFolder_ = this->debugFolder_;

                dissimilarity_pyramid_inverse_[ii] = createDissimilarity(dissimilarity_type_[ii], ii);
                dissimilarity_pyramid_inverse_[ii]->initialize(source_pyramid_[ii]);
                dissimilarity_pyramid_inverse_[ii]->debugFolder_ = this->debugFolder_;
            }

            if ( !debugFolder_.empty() )
            {
                for ( ii=0; ii<resolution_pyramid_levels_; ii++ )
                {
                    std::ostringstream ostr_t;
                    ostr_t << "target_" << ii;

                    gt_exporter_.export_image(target_pyramid_[ii], debugFolder_+ostr_t.str());

                    std::ostringstream ostr_s;
                    ostr_s << "source_" << ii;

                    gt_exporter_.export_image(source_pyramid_[ii], debugFolder_+ostr_s.str());
                }
            }
        }
        catch(...)
        {
            GERROR_STREAM("Errors happened in hoImageRegRegister<TargetType, SourceType, CoordType>::initialize() ... ");
        }

        return true;
    }

    template<typename TargetType, typename SourceType, typename CoordType> 
    inline void hoImageRegRegister<TargetType, SourceType, CoordType>::setTarget(TargetType& target)
    {
        target_ = &target;
    }

    template<typename TargetType, typename SourceType, typename CoordType> 
    inline void hoImageRegRegister<TargetType, SourceType, CoordType>::setSource(SourceType& source)
    {
        source_ = &source;
    }

    template<typename TargetType, typename SourceType, typename CoordType> 
    void hoImageRegRegister<TargetType, SourceType, CoordType>::printContent(std::ostream& os) const
    {
        using namespace std;
        os << "Input dimension is : " << DIn << endl;
        os << "Output dimension is : " << DOut << endl;

        std::string elemTypeName = std::string(typeid(ValueType).name());
        os << "Image data type is : " << elemTypeName << std::endl;

        elemTypeName = std::string(typeid(CoordType).name());
        os << "Transformation coordinate data type is : " << elemTypeName << std::endl;

        os << "Whether to perform world coordinate registration is : " << use_world_coordinates_ << std::endl;
        os << "Number of resolution pyramid levels is : " << resolution_pyramid_levels_ << std::endl;

        os << "------------" << std::endl;
        os << "Downsample ratio of the resolution pyramid for every dimension and every level is : " << std::endl;

        unsigned int ii, jj;
        for ( ii=0; ii<resolution_pyramid_levels_-1; ii++ )
        {
            os << "Level " << ii << " [ ";
            for ( jj=0; jj<resolution_pyramid_downsample_ratio_[ii].size(); jj++ )
            {
                os << resolution_pyramid_downsample_ratio_[ii][jj] << " ";
            }
            os << " ] " << std::endl;
        }

        os << "------------" << std::endl;
        os << "Gaussian blurring sigma for every dimension and every level is : " << std::endl;
        for ( ii=0; ii<resolution_pyramid_levels_; ii++ )
        {
            os << "Level " << ii << " [ ";
            for ( jj=0; jj<resolution_pyramid_blurring_sigma_[ii].size(); jj++ )
            {
                os << resolution_pyramid_blurring_sigma_[ii][jj] << " ";
            }
            os << " ] " << std::endl;
        }

        os << "------------" << std::endl;
        os << "Boundary handler and interpolator type for warper is : " << std::endl;
        for ( ii=0; ii<resolution_pyramid_levels_; ii++ )
        {
            os << " Level " << ii << " - " 
                << getBoundaryHandlerName(boundary_handler_type_warper_[ii]) 
                << " - " << getInterpolatorName(interp_type_warper_[ii]) << std::endl;
        }

        os << "------------" << std::endl;
        os << "Boundary handler and interpolator type for pyramid construction is : " << std::endl;
        for ( ii=0; ii<resolution_pyramid_levels_; ii++ )
        {
            os << " Level " << ii << " - " 
                << getBoundaryHandlerName(boundary_handler_type_pyramid_construction_) 
                << " - " << getInterpolatorName(interp_type_pyramid_construction_) << std::endl;
        }

        os << "------------" << std::endl;
        os << "Image dissimilarity is : " << std::endl;
        for ( ii=0; ii<resolution_pyramid_levels_; ii++ )
        {
            os << " Level " << ii << " - " 
                << getDissimilarityName(dissimilarity_type_[ii]) << std::endl;
        }

        os << "------------" << std::endl;
        os << "Image registration solver is : " << std::endl;
        for ( ii=0; ii<resolution_pyramid_levels_; ii++ )
        {
            os << " Level " << ii << " - " 
                << getImageRegSolverName(solver_type_[ii]) << std::endl;
        }
    }

    template<typename TargetType, typename SourceType, typename CoordType> 
    void hoImageRegRegister<TargetType, SourceType, CoordType>::print(std::ostream& os) const
    {
        using namespace std;
        os << "--------------Gagdgetron image register -------------" << endl;
        this->printContent(os);
        os << "-----------------------------------------------------" << std::endl;
    }
}
#endif // hoImageRegRegister_H_
