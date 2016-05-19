/** \file   cmr_motion_correction.cpp
    \brief  Implement functionalities to perform some cardiac MR motion correction
    \author Hui Xue
*/

#include "cmr_motion_correction.h"

namespace Gadgetron { 

template <typename T> 
void find_key_frame_SSD_2DT(const Gadgetron::hoNDArray<T>& input, size_t& key_frame)
{
    try
    {
        size_t RO = input.get_size(0);
        size_t E1 = input.get_size(1);

        long long col = input.get_size(2);

        hoMatrixReal<T> SSD(col, col);
        Gadgetron::clear(SSD);

        T* pInput = const_cast<T*>(input.begin());

        long long m, n;
#pragma omp parallel default(none) private(m, n) shared(pInput, RO, E1, col, SSD)
        {
            hoNDArray<T> diff(RO, E1);

#pragma omp for
            for ( m=0; m<col; m++ )
            {
                hoNDArray<T> a(RO, E1, pInput+m*RO*E1);

                T v(0);
                for ( n=m+1; n<col; n++ )
                {
                    hoNDArray<T> b(RO, E1, pInput+n*RO*E1);

                    Gadgetron::subtract(a, b, diff);

                    Gadgetron::norm2(diff, v);
                    SSD(m, n) = v;
                    SSD(n, m) = v;
                }
            }

            diff.clear();
        }

        // sort for every column
        SSD.sort_ascending_along_row();

        // pick the middel row
        hoMatrixReal<T> minimalSSDCol;
        SSD.subMatrix(minimalSSDCol, col/2, col/2, 0, col-1);

        // find minimal SSD
        T minSSD(0);
        Gadgetron::minAbsolute(minimalSSDCol, minSSD, key_frame);

        if ( key_frame >= (size_t)col ) key_frame=(size_t)col-1;
    }
    catch(...)
    {
        GADGET_THROW("Exceptions happened in find_key_frame_SSD(...) ... ");
    }
}

template EXPORTCMR void find_key_frame_SSD_2DT(const Gadgetron::hoNDArray<float>& input, size_t& key_frame);
template EXPORTCMR void find_key_frame_SSD_2DT(const Gadgetron::hoNDArray<double>& input, size_t& key_frame);

// ------------------------------------------------------------------------

template <typename T> 
void perform_moco_fixed_key_frame_2DT(const Gadgetron::hoNDArray<T>& input, size_t key_frame, T reg_strength, std::vector<unsigned int> iters, bool bidirectional_moco, Gadgetron::hoImageRegContainer2DRegistration<T, float, 2, 2>& reg)
{
    try
    {
        size_t RO = input.get_size(0);
        size_t E1 = input.get_size(1);
        size_t N = input.get_size(2);

        GADGET_CHECK_THROW(key_frame<N);

        Gadgetron::hoNDImageContainer2D< hoNDImage<T, 2> > im;
        std::vector<size_t> cols(1, N);

        std::vector<size_t> dim(3);
        dim[0] = RO;
        dim[1] = E1;
        dim[2] = N;

        im.create(const_cast<T*>(input.begin()), dim);

        size_t div_num = 3;
        float dissimilarity_thres = 1e-6;
        size_t inverse_deform_enforce_iter = 10;
        float inverse_deform_enforce_weight = 0.5;
        size_t level = iters.size();

        reg.setDefaultParameters( (unsigned int)level, false);

        reg.container_reg_mode_ = GT_IMAGE_REG_CONTAINER_FIXED_REFERENCE;
        reg.bg_value_ = -1;

        reg.container_reg_transformation_ = (bidirectional_moco ? GT_IMAGE_REG_TRANSFORMATION_DEFORMATION_FIELD_BIDIRECTIONAL : GT_IMAGE_REG_TRANSFORMATION_DEFORMATION_FIELD);
        reg.max_iter_num_pyramid_level_ = iters;

        reg.boundary_handler_type_warper_.clear();
        reg.boundary_handler_type_warper_.resize(level, GT_BOUNDARY_CONDITION_BORDERVALUE);

        reg.interp_type_warper_.clear();
        reg.interp_type_warper_.resize(level, GT_IMAGE_INTERPOLATOR_LINEAR);

        reg.regularization_hilbert_strength_pyramid_level_.clear();
        reg.regularization_hilbert_strength_pyramid_level_.resize(level);
        for ( size_t ii=0; ii<level; ii++ )
        {
            reg.regularization_hilbert_strength_pyramid_level_[ii].resize(2, reg_strength);
        }

        reg.dissimilarity_type_ = GT_IMAGE_DISSIMILARITY_LocalCCR;
        reg.dissimilarity_thres_pyramid_level_.clear();

        reg.dissimilarity_thres_pyramid_level_.resize(level, dissimilarity_thres);

        reg.inverse_deform_enforce_iter_pyramid_level_.clear();
        reg.inverse_deform_enforce_iter_pyramid_level_.resize(level, inverse_deform_enforce_iter);

        reg.inverse_deform_enforce_weight_pyramid_level_.clear();
        reg.inverse_deform_enforce_weight_pyramid_level_.resize(level, inverse_deform_enforce_weight);

        reg.div_num_pyramid_level_.clear();
        reg.div_num_pyramid_level_.resize(level, div_num);

        std::vector<unsigned int> referenceFrame(1, key_frame);
        reg.registerOverContainer2DFixedReference(im, referenceFrame, false, false);
    }
    catch(...)
    {
        GADGET_THROW("Exceptions happened in perform_moco_fixed_key_frame_2DT(...) ... ");
    }
}

template EXPORTCMR void perform_moco_fixed_key_frame_2DT(const Gadgetron::hoNDArray<float>& input, size_t key_frame, float reg_strength, std::vector<unsigned int> iters, bool bidirectional_moco, Gadgetron::hoImageRegContainer2DRegistration<float, float, 2, 2>& reg);

// ------------------------------------------------------------------------

template <typename T> 
void perform_moco_pair_wise_frame_2DT(const Gadgetron::hoNDArray<T>& target, const Gadgetron::hoNDArray<T>& source, 
    size_t key_frame, T reg_strength, std::vector<unsigned int> iters, bool bidirectional_moco, Gadgetron::hoImageRegContainer2DRegistration<T, float, 2, 2>& reg)
{
    try
    {
        size_t RO = target.get_size(0);
        size_t E1 = target.get_size(1);
        size_t N = target.get_size(2);

        GADGET_CHECK_THROW(source.get_size(0)==RO);
        GADGET_CHECK_THROW(source.get_size(1)==E1);
        GADGET_CHECK_THROW(source.get_size(2)==N);

        Gadgetron::hoNDImageContainer2D< hoNDImage<T, 2> > im_target, im_source;
        std::vector<size_t> cols(1, N);

        std::vector<size_t> dim(3);
        dim[0] = RO;
        dim[1] = E1;
        dim[2] = N;

        im_target.create(const_cast<T*>(target.begin()), dim);
        im_source.create(const_cast<T*>(source.begin()), dim);

        size_t div_num = 3;
        float dissimilarity_thres = 1e-6;
        size_t inverse_deform_enforce_iter = 10;
        float inverse_deform_enforce_weight = 0.5;
        size_t level = iters.size();

        reg.setDefaultParameters( (unsigned int)level, false);

        reg.container_reg_mode_ = GT_IMAGE_REG_CONTAINER_PAIR_WISE;
        reg.bg_value_ = -1;

        reg.container_reg_transformation_ = (bidirectional_moco ? GT_IMAGE_REG_TRANSFORMATION_DEFORMATION_FIELD_BIDIRECTIONAL : GT_IMAGE_REG_TRANSFORMATION_DEFORMATION_FIELD);
        reg.max_iter_num_pyramid_level_ = iters;

        reg.boundary_handler_type_warper_.clear();
        reg.boundary_handler_type_warper_.resize(level, GT_BOUNDARY_CONDITION_BORDERVALUE);

        reg.interp_type_warper_.clear();
        reg.interp_type_warper_.resize(level, GT_IMAGE_INTERPOLATOR_LINEAR);

        reg.regularization_hilbert_strength_pyramid_level_.clear();
        reg.regularization_hilbert_strength_pyramid_level_.resize(level);
        for ( size_t ii=0; ii<level; ii++ )
        {
            reg.regularization_hilbert_strength_pyramid_level_[ii].resize(2, reg_strength);
        }

        reg.dissimilarity_type_ = GT_IMAGE_DISSIMILARITY_LocalCCR;
        reg.dissimilarity_thres_pyramid_level_.clear();

        reg.dissimilarity_thres_pyramid_level_.resize(level, dissimilarity_thres);

        reg.inverse_deform_enforce_iter_pyramid_level_.clear();
        reg.inverse_deform_enforce_iter_pyramid_level_.resize(level, inverse_deform_enforce_iter);

        reg.inverse_deform_enforce_weight_pyramid_level_.clear();
        reg.inverse_deform_enforce_weight_pyramid_level_.resize(level, inverse_deform_enforce_weight);

        reg.div_num_pyramid_level_.clear();
        reg.div_num_pyramid_level_.resize(level, div_num);

        std::vector<unsigned int> referenceFrame(1, key_frame);
        reg.registerOverContainer2DPairWise(im_target, im_source, false, false);
    }
    catch(...)
    {
        GADGET_THROW("Exceptions happened in perform_moco_pair_wise_frame_2DT(...) ... ");
    }
}

template EXPORTCMR void perform_moco_pair_wise_frame_2DT(const Gadgetron::hoNDArray<float>& target, const Gadgetron::hoNDArray<float>& source, size_t key_frame, float reg_strength, std::vector<unsigned int> iters, bool bidirectional_moco, Gadgetron::hoImageRegContainer2DRegistration<float, float, 2, 2>& reg);

}
