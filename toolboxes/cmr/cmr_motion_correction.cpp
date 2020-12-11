/** \file   cmr_motion_correction.cpp
    \brief  Implement functionalities to perform some cardiac MR motion correction
    \author Hui Xue
*/

#include "cmr_motion_correction.h"

#include "ismrmrd/ismrmrd.h"
#include "ismrmrd/meta.h"

#include "boost/tuple/tuple.hpp"
#include "boost/tuple/tuple_comparison.hpp"
#include "boost/tuple/tuple_io.hpp"

#include "ho2DArray.h"
#include "ho3DArray.h"
#include "ho4DArray.h"
#include "ho5DArray.h"
#include "ho6DArray.h"
#include "ho7DArray.h"
#include "hoNDArray_linalg.h"
#include "hoNDFFT.h"
#include "hoNDKLT.h"
#include "hoNDArray_utils.h"
#include "hoNDArray_elemwise.h"
#include "hoNDImage_util.h"
#include "hoMRImage.h"
#include "hoNDArray_reductions.h"
#include "hoNDArray_linalg.h"
#include "hoNDFFT.h"

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

                    v = Gadgetron::nrm2(diff);
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
void find_key_frame_SSD_2DT(const Gadgetron::hoNDArray<T>& input, size_t& key_frame, std::vector< std::pair<double, size_t> >& moco_quality)
{
    try
    {
        size_t RO = input.get_size(0);
        size_t E1 = input.get_size(1);

        long long col = input.get_size(2);

        hoMatrixReal<T> SSD(col, col);
        Gadgetron::clear(SSD);

        T* pInput = const_cast<T*>(input.begin());

        moco_quality.resize(col);

        long long m, n;
#pragma omp parallel default(none) private(m, n) shared(pInput, RO, E1, col, SSD)
        {
            hoNDArray<T> diff(RO, E1);

#pragma omp for
            for (m = 0; m < col; m++)
            {
                hoNDArray<T> a(RO, E1, pInput + m*RO*E1);

                T v(0);
                for (n = m + 1; n < col; n++)
                {
                    hoNDArray<T> b(RO, E1, pInput + n*RO*E1);

                    Gadgetron::subtract(a, b, diff);

                    v = Gadgetron::nrm2(diff);
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
        SSD.subMatrix(minimalSSDCol, col / 2, col / 2, 0, col - 1);

        // find minimal SSD
        T minSSD(0);
        Gadgetron::minAbsolute(minimalSSDCol, minSSD, key_frame);

        if (key_frame >= (size_t)col) key_frame = (size_t)col - 1;

        hoNDArray<T> diff(RO, E1);
        hoNDArray<T> a(RO, E1, pInput + key_frame*RO*E1);

        for (m = 0; m < col; m++)
        {
            /*if(m==key_frame)
            moco_quality[m].first = 0;
            else
            {
            hoNDArray<T> b(RO, E1, pInput + m*RO*E1);
            Gadgetron::subtract(a, b, diff);
            Gadgetron::norm2(diff, v);
            moco_quality[m].first = v;
            }*/

            moco_quality[m].first = minimalSSDCol(0, m);
            moco_quality[m].second = m;
        }

        std::sort(moco_quality.begin(), moco_quality.end(), cmr_moco_ave_compObj());
    }
    catch (...)
    {
        GADGET_THROW("Exceptions happened in find_key_frame_SSD_2DT(moco_quality) ... ");
    }
}

template EXPORTCMR void find_key_frame_SSD_2DT(const Gadgetron::hoNDArray<float>& input, size_t& key_frame, std::vector< std::pair<double, size_t> >& moco_quality);
template EXPORTCMR void find_key_frame_SSD_2DT(const Gadgetron::hoNDArray<double>& input, size_t& key_frame, std::vector< std::pair<double, size_t> >& moco_quality);

// ------------------------------------------------------------------------

template <typename T>
void compute_SSD_2DT(const Gadgetron::hoNDArray<T>& input, size_t key_frame, std::vector< std::pair<double, size_t> >& moco_quality)
{
    try
    {
        size_t RO = input.get_size(0);
        size_t E1 = input.get_size(1);

        long long col = input.get_size(2);
        T* pInput = const_cast<T*>(input.begin());

        hoNDArray<T> diff(RO, E1);
        hoNDArray<T> a(RO, E1, pInput + key_frame*RO*E1);

        moco_quality.resize(col);

        T v(0);
        size_t m;
        for (m = 0; m < col; m++)
        {
            if (m == key_frame)
                moco_quality[m].first = 0;
            else
            {
                hoNDArray<T> b(RO, E1, pInput + m*RO*E1);
                Gadgetron::subtract(a, b, diff);
                v = Gadgetron::nrm2(diff);
                moco_quality[m].first = v;
            }

            moco_quality[m].second = m;
        }

        std::sort(moco_quality.begin(), moco_quality.end(), cmr_moco_ave_compObj());
    }
    catch (...)
    {
        GADGET_THROW("Exceptions happened in compute_SSD_2DT(moco_quality) ... ");
    }
}

template EXPORTCMR void compute_SSD_2DT(const Gadgetron::hoNDArray<float>& input, size_t key_frame, std::vector< std::pair<double, size_t> >& moco_quality);
template EXPORTCMR void compute_SSD_2DT(const Gadgetron::hoNDArray<double>& input, size_t key_frame, std::vector< std::pair<double, size_t> >& moco_quality);

// ------------------------------------------------------------------------

template <typename T>
void compute_deformation_jacobian(const Gadgetron::hoNDArray<T>& dx, const Gadgetron::hoNDArray<T>& dy, std::vector<T>& mean_deform, std::vector<T>& max_deform, std::vector<T>& mean_log_jac, std::vector<T>& max_log_jac)
{
    try
    {
        hoImageRegDeformationField<T, 2> deformTransform;
        T meanDeform, maxDeform, meanLogJac, maxLogJac;
        hoNDArray<T> jac;
        hoNDImage<T, 2> d[2];
        hoNDImage<T, 2>* deform_field[2];

        size_t c;

        T minMeanDeform = (T)1e6;
        T minMeanLogJac = (T)1e6;
        size_t refSource = 0;

        size_t RO = dx.get_size(0);
        size_t E1 = dx.get_size(1);
        size_t N = dx.get_size(2);

        T* pDx = const_cast<T*>(dx.begin());
        T* pDy = const_cast<T*>(dy.begin());

        std::vector<size_t> dim(2);
        dim[0] = RO;
        dim[1] = E1;

        mean_deform.resize(N, 0);
        max_deform.resize(N, 0);
        mean_log_jac.resize(N, 0);
        max_log_jac.resize(N, 0);

        for (c = 0; c < N; c++)
        {
            d[0].create(dim, pDx + c*RO*E1);
            d[1].create(dim, pDy + c*RO*E1);

            deform_field[0] = &d[0];
            deform_field[1] = &d[1];

            deformTransform.jacobianPosition(jac, deform_field, 2);
            deformTransform.analyzeJacobianAndDeformation(jac, deform_field, meanDeform, maxDeform, meanLogJac, maxLogJac, 2);

            mean_deform[c] = meanDeform;
            max_deform[c] = maxDeform;
            mean_log_jac[c] = meanLogJac;
            max_log_jac[c] = maxLogJac;
        }
    }
    catch (...)
    {
        GADGET_THROW("Exceptions happened in compute_deformation_jacobian(...) ... ");
    }
}

template EXPORTCMR void compute_deformation_jacobian(const Gadgetron::hoNDArray<float>& dx, const Gadgetron::hoNDArray<float>& dy, std::vector<float>& mean_deform, std::vector<float>& max_deform, std::vector<float>& mean_log_jac, std::vector<float>& max_log_jac);

// ------------------------------------------------------------------------

template <typename T>
void find_key_frame_deformation_2DT(const Gadgetron::hoNDArray<T>& dx, const Gadgetron::hoNDArray<T>& dy, size_t& key_frame, std::vector< std::pair<double, size_t> >& moco_quality)
{
    try
    {
        std::vector<T> mean_deform, max_deform, mean_log_jac, max_log_jac;
        compute_deformation_jacobian(dx, dy, mean_deform, max_deform, mean_log_jac, max_log_jac);

        size_t N = dx.get_size(2);

        size_t c;
        for (c = 0; c < N; c++)
        {
            moco_quality[c].first = (double)mean_deform[c];
            moco_quality[c].second = c;
        }

        std::sort(moco_quality.begin(), moco_quality.end(), cmr_moco_ave_compObj());
    }
    catch (...)
    {
        GADGET_THROW("Exceptions happened in find_key_frame_deformation_2DT(...) ... ");
    }
}

template EXPORTCMR void find_key_frame_deformation_2DT(const Gadgetron::hoNDArray<float>& dx, const Gadgetron::hoNDArray<float>& dy, size_t& key_frame, std::vector< std::pair<double, size_t> >& moco_quality);

// ------------------------------------------------------------------------

template <typename T>
void perform_averaging(const Gadgetron::hoNDArray<T>& input, size_t& key_frame, const std::vector< std::pair<double, size_t> >& moco_quality, double percentage_kept_for_averaging, Gadgetron::hoNDArray<T>& ave)
{
    try
    {
        size_t RO = input.get_size(0);
        size_t E1 = input.get_size(1);
        size_t N = input.get_size(2);

        ave.create(RO, E1);

        size_t usedImageForAve = (size_t)(N * percentage_kept_for_averaging);
        if (usedImageForAve < 1)
        {
            usedImageForAve = 0;
        }

        memcpy(ave.begin(), input.begin() + moco_quality[0].second*RO*E1, sizeof(T)*RO*E1);

        T* pInput = const_cast<T*>(input.begin());

        size_t ii;
        for (ii = 1; ii < usedImageForAve; ii++)
        {
            hoNDArray<T> im2D(RO, E1, pInput + moco_quality[ii].second*RO*E1);
            Gadgetron::add(ave, im2D, ave);
        }

        Gadgetron::scal((typename realType<T>::Type)(1.0 / std::sqrt((double)usedImageForAve)), ave);
    }
    catch (...)
    {
        GADGET_THROW("Exceptions happened in perform_averaging(...) ... ");
    }
}

template EXPORTCMR void perform_averaging(const Gadgetron::hoNDArray<float>& input, size_t& key_frame, const std::vector< std::pair<double, size_t> >& moco_quality, double percentage_kept_for_averaging, Gadgetron::hoNDArray<float>& ave);

template EXPORTCMR void perform_averaging(const Gadgetron::hoNDArray<double>& input, size_t& key_frame, const std::vector< std::pair<double, size_t> >& moco_quality, double percentage_kept_for_averaging, Gadgetron::hoNDArray<double>& ave);

template EXPORTCMR void perform_averaging(const Gadgetron::hoNDArray< std::complex<float> >& input, size_t& key_frame, const std::vector< std::pair<double, size_t> >& moco_quality, double percentage_kept_for_averaging, Gadgetron::hoNDArray< std::complex<float> >& ave);

template EXPORTCMR void perform_averaging(const Gadgetron::hoNDArray< std::complex<double> >& input, size_t& key_frame, const std::vector< std::pair<double, size_t> >& moco_quality, double percentage_kept_for_averaging, Gadgetron::hoNDArray< std::complex<double> >& ave);

// ------------------------------------------------------------------------

template <typename T>
void perform_moco_fixed_key_frame_2DT(Gadgetron::hoNDImageContainer2D< hoNDImage<T, 2> >& input, const std::vector<unsigned int>& key_frame, T reg_strength, std::vector<unsigned int> iters, bool bidirectional_moco, bool warp_input, Gadgetron::hoImageRegContainer2DRegistration<Gadgetron::hoNDImage<T, 2>, Gadgetron::hoNDImage<T, 2>, double>& reg)
{
    try
    {
        size_t div_num = 3;
        float dissimilarity_thres = 1e-6;
        size_t inverse_deform_enforce_iter = 10;
        float inverse_deform_enforce_weight = 0.5;
        size_t level = iters.size();

        reg.setDefaultParameters((unsigned int)level, false);

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
        for (size_t ii = 0; ii<level; ii++)
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

        reg.registerOverContainer2DFixedReference(input, key_frame, warp_input, false);
    }
    catch (...)
    {
        GADGET_THROW("Exceptions happened in perform_moco_fixed_key_frame_2DT(container) ... ");
    }
}

template EXPORTCMR void perform_moco_fixed_key_frame_2DT(Gadgetron::hoNDImageContainer2D< hoNDImage<float, 2> >& input, const std::vector<unsigned int>& key_frame, float reg_strength, std::vector<unsigned int> iters, bool bidirectional_moco, bool warp_input, Gadgetron::hoImageRegContainer2DRegistration<Gadgetron::hoNDImage<float, 2>, Gadgetron::hoNDImage<float, 2>, double>& reg);

template EXPORTCMR void perform_moco_fixed_key_frame_2DT(Gadgetron::hoNDImageContainer2D< hoNDImage<double, 2> > & input, const std::vector<unsigned int>& key_frame, double reg_strength, std::vector<unsigned int> iters, bool bidirectional_moco, bool warp_input, Gadgetron::hoImageRegContainer2DRegistration<Gadgetron::hoNDImage<double, 2>, Gadgetron::hoNDImage<double, 2>, double>& reg);

// ------------------------------------------------------------------------

template <typename T>
void perform_moco_fixed_key_frame_2DT(const Gadgetron::hoNDArray<T>& input, size_t key_frame, bool warp_input, Gadgetron::hoImageRegContainer2DRegistration<Gadgetron::hoNDImage<T, 2>, Gadgetron::hoNDImage<T, 2>, double>& reg)
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

    std::vector<unsigned int> referenceFrame(1, key_frame);

    reg.registerOverContainer2DFixedReference(im, referenceFrame, warp_input, false);
}

template EXPORTCMR void perform_moco_fixed_key_frame_2DT(const Gadgetron::hoNDArray<float>& input, size_t key_frame, bool warp_input, Gadgetron::hoImageRegContainer2DRegistration<Gadgetron::hoNDImage<float, 2>, Gadgetron::hoNDImage<float, 2>, double>& reg);

template EXPORTCMR void perform_moco_fixed_key_frame_2DT(const Gadgetron::hoNDArray<double>& input, size_t key_frame, bool warp_input, Gadgetron::hoImageRegContainer2DRegistration<Gadgetron::hoNDImage<double, 2>, Gadgetron::hoNDImage<double, 2>, double>& reg);

// ------------------------------------------------------------------------

template <typename T>
void perform_moco_fixed_key_frame_2DT(const Gadgetron::hoNDArray<T>& input, size_t key_frame, T reg_strength, std::vector<unsigned int> iters, bool bidirectional_moco, bool warp_input, Gadgetron::hoImageRegContainer2DRegistration<Gadgetron::hoNDImage<T, 2>, Gadgetron::hoNDImage<T, 2>, double>& reg)
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

        std::vector<unsigned int> referenceFrame(1, key_frame);

        Gadgetron::perform_moco_fixed_key_frame_2DT(im, referenceFrame, reg_strength, iters, bidirectional_moco, warp_input, reg);
    }
    catch (...)
    {
        GADGET_THROW("Exceptions happened in perform_moco_fixed_key_frame_2DT(...) ... ");
    }
}

template EXPORTCMR void perform_moco_fixed_key_frame_2DT(const Gadgetron::hoNDArray<float>& input, size_t key_frame, float reg_strength, std::vector<unsigned int> iters, bool bidirectional_moco, bool warp_input, Gadgetron::hoImageRegContainer2DRegistration<Gadgetron::hoNDImage<float, 2>, Gadgetron::hoNDImage<float, 2>, double>& reg);

template EXPORTCMR void perform_moco_fixed_key_frame_2DT(const Gadgetron::hoNDArray<double>& input, size_t key_frame, double reg_strength, std::vector<unsigned int> iters, bool bidirectional_moco, bool warp_input, Gadgetron::hoImageRegContainer2DRegistration<Gadgetron::hoNDImage<double, 2>, Gadgetron::hoNDImage<double, 2>, double>& reg);

// ------------------------------------------------------------------------

template <typename T>
void perform_moco_pair_wise_frame_2DT(const Gadgetron::hoNDArray<T>& target, const Gadgetron::hoNDArray<T>& source, bool warp_input, Gadgetron::hoImageRegContainer2DRegistration<Gadgetron::hoNDImage<T, 2>, Gadgetron::hoNDImage<T, 2>, double>& reg)
{
    try
    {
        size_t RO = target.get_size(0);
        size_t E1 = target.get_size(1);
        size_t N = target.get_size(2);

        GADGET_CHECK_THROW(source.get_size(0) == RO);
        GADGET_CHECK_THROW(source.get_size(1) == E1);
        GADGET_CHECK_THROW(source.get_size(2) == N);

        Gadgetron::hoNDImageContainer2D< hoNDImage<T, 2> > im_target, im_source;
        std::vector<size_t> cols(1, N);

        std::vector<size_t> dim(3);
        dim[0] = RO;
        dim[1] = E1;
        dim[2] = N;

        im_target.create(const_cast<T*>(target.begin()), dim);
        im_source.create(const_cast<T*>(source.begin()), dim);

        reg.registerOverContainer2DPairWise(im_target, im_source, warp_input, false);
    }
    catch (...)
    {
        GADGET_THROW("Exceptions happened in perform_moco_pair_wise_frame_2DT(...) ... ");
    }
}

template EXPORTCMR void perform_moco_pair_wise_frame_2DT(const Gadgetron::hoNDArray<float>& target, const Gadgetron::hoNDArray<float>& source, bool warp_input, Gadgetron::hoImageRegContainer2DRegistration<Gadgetron::hoNDImage<float, 2>, Gadgetron::hoNDImage<float, 2>, double>& reg);

template EXPORTCMR void perform_moco_pair_wise_frame_2DT(const Gadgetron::hoNDArray<double>& target, const Gadgetron::hoNDArray<double>& source, bool warp_input, Gadgetron::hoImageRegContainer2DRegistration<Gadgetron::hoNDImage<double, 2>, Gadgetron::hoNDImage<double, 2>, double>& reg);

// ------------------------------------------------------------------------

template <typename T>
void perform_moco_pair_wise_frame_2DT(const Gadgetron::hoNDArray<T>& target, const Gadgetron::hoNDArray<T>& source,
    T reg_strength, std::vector<unsigned int> iters, bool bidirectional_moco, bool warp_input, Gadgetron::hoImageRegContainer2DRegistration<Gadgetron::hoNDImage<T, 2>, Gadgetron::hoNDImage<T, 2>, double>& reg)
{
    try
    {
        size_t RO = target.get_size(0);
        size_t E1 = target.get_size(1);
        size_t N = target.get_size(2);

        GADGET_CHECK_THROW(source.get_size(0) == RO);
        GADGET_CHECK_THROW(source.get_size(1) == E1);
        GADGET_CHECK_THROW(source.get_size(2) == N);

        size_t div_num = 3;
        float dissimilarity_thres = 1e-6;
        size_t inverse_deform_enforce_iter = 10;
        float inverse_deform_enforce_weight = 0.5;
        size_t level = iters.size();

        reg.setDefaultParameters((unsigned int)level, false);

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
        for (size_t ii = 0; ii<level; ii++)
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

        Gadgetron::perform_moco_pair_wise_frame_2DT(target, source, warp_input, reg);
    }
    catch (...)
    {
        GADGET_THROW("Exceptions happened in perform_moco_pair_wise_frame_2DT(...) ... ");
    }
}

template EXPORTCMR void perform_moco_pair_wise_frame_2DT(const Gadgetron::hoNDArray<float>& target, const Gadgetron::hoNDArray<float>& source, float reg_strength, std::vector<unsigned int> iters, bool bidirectional_moco, bool warp_input, Gadgetron::hoImageRegContainer2DRegistration<Gadgetron::hoNDImage<float, 2>, Gadgetron::hoNDImage<float, 2>, double>& reg);

template EXPORTCMR void perform_moco_pair_wise_frame_2DT(const Gadgetron::hoNDArray<double>& target, const Gadgetron::hoNDArray<double>& source, double reg_strength, std::vector<unsigned int> iters, bool bidirectional_moco, bool warp_input, Gadgetron::hoImageRegContainer2DRegistration<Gadgetron::hoNDImage<double, 2>, Gadgetron::hoNDImage<double, 2>, double>& reg);

// ------------------------------------------------------------------------

template <typename T>
void apply_deformation_field(const Gadgetron::hoNDArray<T>& target, const Gadgetron::hoNDArray<double>& dx, const Gadgetron::hoNDArray<double>& dy, Gadgetron::hoNDArray<T>& output, Gadgetron::GT_BOUNDARY_CONDITION bh)
{
    try
    {
        size_t RO = target.get_size(0);
        size_t E1 = target.get_size(1);
        size_t N = target.get_size(2);

        GADGET_CHECK_THROW(dx.get_size(0) == RO);
        GADGET_CHECK_THROW(dx.get_size(1) == E1);
        GADGET_CHECK_THROW(dx.get_size(2) == N);

        GADGET_CHECK_THROW(dy.get_size(0) == RO);
        GADGET_CHECK_THROW(dy.get_size(1) == E1);
        GADGET_CHECK_THROW(dy.get_size(2) == N);

        Gadgetron::hoNDImageContainer2D< hoNDImage<T, 2> > im_target, im_source;
        std::vector<size_t> cols(1, N);

        std::vector<size_t> dim(3);
        dim[0] = RO;
        dim[1] = E1;
        dim[2] = N;

        im_target.create(const_cast<T*>(target.begin()), dim);

        im_source.copyFrom(im_target);

        Gadgetron::hoImageRegContainer2DRegistration<hoNDImage<T, 2>, hoNDImage<T, 2>, double> reg;

        hoNDImageContainer2D< hoNDImage<double, 2> > deformation_field[2];
        deformation_field[0].create(const_cast<double*>(dx.begin()), dim);
        deformation_field[1].create(const_cast<double*>(dy.begin()), dim);

        reg.warpContainer2D(im_target, im_target, deformation_field, im_source, bh);

        im_source.to_NDArray(0, output);
    }
    catch (...)
    {
        GADGET_THROW("Exceptions happened in apply_deformation_field(...) ... ");
    }
}

template EXPORTCMR void apply_deformation_field(const Gadgetron::hoNDArray<float>& input, const Gadgetron::hoNDArray<double>& dx, const Gadgetron::hoNDArray<double>& dy, Gadgetron::hoNDArray<float>& output, Gadgetron::GT_BOUNDARY_CONDITION bh);

template EXPORTCMR void apply_deformation_field(const Gadgetron::hoNDArray<double>& input, const Gadgetron::hoNDArray<double>& dx, const Gadgetron::hoNDArray<double>& dy, Gadgetron::hoNDArray<double>& output, Gadgetron::GT_BOUNDARY_CONDITION bh);

// ------------------------------------------------------------------------

template <typename T>
void concatenate_deform_fields_2DT(const hoNDArray<T>& dx, const hoNDArray<T>& dy, size_t key_frame, hoNDArray<T>& dx_out, hoNDArray<T>& dy_out)
{
    try
    {
        typedef hoNDArray<T> ArrayType;
        typedef hoNDInterpolatorLinear<ArrayType> InterpType;
        typedef hoNDBoundaryHandlerFixedValue<ArrayType> BHType;

        size_t RO = dx.get_size(0);
        size_t E1 = dx.get_size(1);
        size_t N = dx.get_size(2);

        GADGET_CHECK_THROW(dx.dimensions_equal(&dy));
        GADGET_CHECK_THROW(key_frame >= 0);
        GADGET_CHECK_THROW(key_frame < N);

        dx_out = dx;
        dy_out = dy;

        std::vector<InterpType> interp_ro(N), interp_e1(N);
        std::vector<BHType> bhFixedValue_ro(N), bhFixedValue_e1(N);
        std::vector<ArrayType> curr_dx(N), curr_dy(N);

        for (size_t n = 0; n < N; n++)
        {
            curr_dx[n].create(RO, E1, const_cast<T*>(dx.begin() + n * RO * E1));
            curr_dy[n].create(RO, E1, const_cast<T*>(dy.begin() + n * RO * E1));

            interp_ro[n].setBoundaryHandler(bhFixedValue_ro[n]);
            interp_e1[n].setBoundaryHandler(bhFixedValue_e1[n]);

            interp_ro[n].setArray(curr_dx[n]);
            interp_e1[n].setArray(curr_dy[n]);

            bhFixedValue_ro[n].setArray(curr_dx[n]);
            bhFixedValue_e1[n].setArray(curr_dy[n]);
        }

        double p_ro, p_e1;

        long long i;

        #pragma omp parallel for default(none) private(i, p_ro, p_e1) shared(key_frame, N, dx_out, dy_out, RO, E1, interp_ro, interp_e1)
        for (i = key_frame + 2; i < N; i++)
        {
            long ro, e1, j;
            for (ro = 0; ro < RO; ro++)
            {
                for (e1 = 0; e1 < E1; e1++)
                {
                    // forward
                    p_ro = ro; p_e1 = e1;
                    for (j = key_frame + 1; j < i; j++)
                    {
                        p_ro += interp_ro[j](p_ro, p_e1);
                        p_e1 += interp_e1[j](p_ro, p_e1);
                    }

                    dx_out(ro, e1, i) = p_ro + interp_ro[j](p_ro, p_e1) - ro;
                    dy_out(ro, e1, i) = p_e1 + interp_e1[j](p_ro, p_e1) - e1;
                }
            }
        }

        #pragma omp parallel for default(none) private(i, p_ro, p_e1) shared(key_frame, N, dx_out, dy_out, RO, E1, interp_ro, interp_e1)
        for (i = key_frame - 2; i >= 0; i--)
        {
            long ro, e1, j;
            for (ro = 0; ro < RO; ro++)
            {
                for (e1 = 0; e1 < E1; e1++)
                {
                    p_ro = ro; p_e1 = e1;
                    for (j = key_frame - 1; j > i; j--)
                    {
                        p_ro += interp_ro[j](p_ro, p_e1);
                        p_e1 += interp_e1[j](p_ro, p_e1);
                    }

                    dx_out(ro, e1, i) = p_ro + interp_ro[j](p_ro, p_e1) - ro;
                    dy_out(ro, e1, i) = p_e1 + interp_e1[j](p_ro, p_e1) - e1;
                }
            }
        }

    }
    catch (...)
    {
        GADGET_THROW("Error happened in concatenate_deform_fields_2DT(...) ... ");
    }
}

template EXPORTCMR void concatenate_deform_fields_2DT(const hoNDArray<float>& dx, const hoNDArray<float>& dy, size_t key_frame, hoNDArray<float>& dx_out, hoNDArray<float>& dy_out);
template EXPORTCMR void concatenate_deform_fields_2DT(const hoNDArray<double>& dx, const hoNDArray<double>& dy, size_t key_frame, hoNDArray<double>& dx_out, hoNDArray<double>& dy_out);

// ------------------------------------------------------------------------

template <typename T>
void find_key_frame_use_deformation_cross_series(const Gadgetron::hoNDArray<T>& target, const Gadgetron::hoNDArray<T>& source, Gadgetron::hoImageRegContainer2DRegistration<Gadgetron::hoNDImage<T, 2>, Gadgetron::hoNDImage<T, 2>, double>& reg, size_t& key_frame, std::vector< std::pair<double, size_t> >& moco_quality, hoNDArray<double>& dx, hoNDArray<double>& dy)
{
    try
    {
        size_t RO = target.get_size(0);
        size_t E1 = target.get_size(1);
        size_t N = target.get_size(2);

        GADGET_CHECK_THROW(RO == source.get_size(0));
        GADGET_CHECK_THROW(E1 == source.get_size(1));
        size_t M = source.get_size(2);

        moco_quality.resize(N);

        dx.create(RO, E1, M, N);
        Gadgetron::clear(dx);

        dy.create(RO, E1, M, N);
        Gadgetron::clear(dy);

        T* pTarget = const_cast<T*>(target.begin());
        T* pSource = const_cast<T*>(source.begin());

        // perform moco

        hoNDArray<T> TT, S;
        TT.create(RO, E1, M, N);
        S.create(RO, E1, M, N);

        size_t n, m;
        for (n = 0; n < N; n++)
        {
            for (m = 0; m < M; m++)
            {
                memcpy(&TT(0, 0, m, n), pTarget + n*RO*E1, sizeof(T)*RO*E1);
                memcpy(&S(0, 0, m, n), pSource + m*RO*E1, sizeof(T)*RO*E1);
            }
        }

        //std::string debug_folder = "D:/gtuser/mrprogs/gtprep/ut/result/";

        //ImageIOAnalyze gt_exporter_;
        //gt_exporter_.export_array(TT, debug_folder + "CmrSashaMOCO_TT");
        //gt_exporter_.export_array(S, debug_folder + "CmrSashaMOCO_SS");

        bool warp_input = false;

        hoNDArray<T> T_moco(RO, E1, M*N, TT.begin());
        hoNDArray<T> S_moco(RO, E1, M*N, S.begin());
        Gadgetron::perform_moco_pair_wise_frame_2DT(T_moco, S_moco, warp_input, reg);

        hoNDArray<double> dx_reg, dy_reg;
        reg.deformation_field_[0].to_NDArray(0, dx_reg);
        reg.deformation_field_[1].to_NDArray(0, dy_reg);

        //gt_exporter_.export_array(dx_reg, debug_folder + "CmrSashaMOCO_dx");
        //gt_exporter_.export_array(dy_reg, debug_folder + "CmrSashaMOCO_dy");

        std::vector<double> mean_deform, max_deform, mean_log_jac, max_log_jac;
        compute_deformation_jacobian(dx_reg, dy_reg, mean_deform, max_deform, mean_log_jac, max_log_jac);

        memcpy(dx.begin(), dx_reg.begin(), dx_reg.get_number_of_bytes());
        memcpy(dy.begin(), dy_reg.begin(), dy_reg.get_number_of_bytes());

        for (n = 0; n < N; n++)
        {
            T total_deform = 0;
            for (m = 0; m < M; m++)
            {
                total_deform += mean_deform[m + n*M];
            }

            moco_quality[n].first = total_deform;
            moco_quality[n].second = n;
        }

        std::sort(moco_quality.begin(), moco_quality.end(), cmr_moco_ave_compObj());

        key_frame = moco_quality[0].second;
    }
    catch (...)
    {
        GADGET_THROW("Exceptions happened in find_key_frame_use_deformation_cross_series(...) ... ");
    }
}

template EXPORTCMR void find_key_frame_use_deformation_cross_series(const Gadgetron::hoNDArray<float>& target, const Gadgetron::hoNDArray<float>& source, Gadgetron::hoImageRegContainer2DRegistration<Gadgetron::hoNDImage<float, 2>, Gadgetron::hoNDImage<float, 2>, double>& reg, size_t& key_frame, std::vector< std::pair<double, size_t> >& moco_quality, hoNDArray<double>& dx, hoNDArray<double>& dy);

template EXPORTCMR void find_key_frame_use_deformation_cross_series(const Gadgetron::hoNDArray<double>& target, const Gadgetron::hoNDArray<double>& source, Gadgetron::hoImageRegContainer2DRegistration<Gadgetron::hoNDImage<double, 2>, Gadgetron::hoNDImage<double, 2>, double>& reg, size_t& key_frame, std::vector< std::pair<double, size_t> >& moco_quality, hoNDArray<double>& dx, hoNDArray<double>& dy);

// ------------------------------------------------------------------------
}
