/** \file   cmr_parametric_mapping.cpp
    \brief  Implement CMR parametric mapping for 2D acquisition
    \author Hui Xue
*/

#include "cmr_parametric_mapping.h"
#include "log.h"

#include "hoNDArray_reductions.h"
#include "hoNDArray_elemwise.h"
#include "hoNDArray_math.h"
#include "ho2DArray.h"
#include "ho3DArray.h"

#include "morphology.h"
#include "BSplineFFD2D.h"
#include "hoNDBSpline.h"

#include "hoNDArray_linalg.h"

#include <boost/math/special_functions/sign.hpp>

namespace Gadgetron { 

template <typename T>
void perform_hole_filling(hoNDArray<T>& map, T hole, size_t max_size_of_holes, bool is_8_connected)
{
    try
    {
        hoNDArray<unsigned int> label;

        size_t RO = map.get_size(0);
        size_t E1 = map.get_size(1);

        size_t num = map.get_number_of_elements() / (RO*E1);

        std::vector<size_t> dim(2);
        dim[0] = RO;
        dim[1] = E1;

        size_t t;
        for (t = 0; t < num; t++)
        {
            hoNDImage<T, 2> curr_map;
            curr_map.create(dim, map.begin() + t*RO*E1);

            T v = Gadgetron::nrm2(curr_map);
            if (v <= 1e-3) continue; // empty map

            T maxV;
            Gadgetron::maxValue(curr_map, maxV);

            T minV;
            Gadgetron::minValue(map, minV);

            label.clear();
            GADGET_CATCH_THROW(Gadgetron::bwlabel_2d(curr_map, hole, label, is_8_connected));

            std::vector<unsigned int> labels, areas;
            GADGET_CATCH_THROW(Gadgetron::bwlabel_area_2d(label, labels, areas));

            bool needHoleFilling = false;
            size_t numHoles = labels.size();

            size_t n;
            for (n = 0; n < numHoles; n++)
            {
                if (areas[n] <= max_size_of_holes)
                {
                    needHoleFilling = true;
                    break;
                }
            }

            if (needHoleFilling)
            {
                size_t gridSize[2];
                gridSize[0] = 3;
                gridSize[1] = 3;

                size_t numOfRefinement = 0;
                size_t numOfControlPts = gridSize[0];
                while (numOfControlPts < RO)
                {
                    numOfRefinement++;
                    numOfControlPts = 2 * numOfControlPts - 1;
                }

                typedef typename Gadgetron::BSplineFFD2D<T, double, 1>::MaskArrayType MaskArrayType;

                MaskArrayType mask;
                mask.create(RO, E1);
                Gadgetron::fill(mask, float(1));

                for (n = 0; n < RO*E1; n++)
                {
                    unsigned int l = label(n);

                    if (l > 0)
                    {
                        size_t p;
                        for (p = 0; p < labels.size(); p++)
                        {
                            if (l == labels[p]) break;
                        }

                        if (areas[p] <= max_size_of_holes)
                        {
                            mask(n) = 0;
                        }
                    }
                }

                hoNDImage<double, 2> curr_map_used;
                curr_map_used.copyFrom(curr_map);

                double totalResidual(0);
                Gadgetron::BSplineFFD2D<double, double, 1> ffd(curr_map_used, gridSize[0], gridSize[1]);

                ffd.ffdApproxImage(&curr_map_used, mask, totalResidual, numOfRefinement);

                hoNDImage<double, 2> flowDual(curr_map_used);
                ffd.evaluateFFDOnImage(flowDual);

                for (n = 0; n < RO*E1; n++)
                {
                    if (mask(n) == 0)
                    {
                        std::vector<size_t> pt(3);
                        curr_map_used.calculate_index(n, pt);

                        T v(0);
                        v = flowDual(n);
                        if (v >= 0.9*minV && v <= 1.1*maxV) curr_map(n) = (T)(v);
                    }
                }
            }
        }
    }
    catch (...)
    {
        GADGET_THROW("Error happened in perform_hole_filling(...) ... ");
    }
}

template EXPORTCMR void perform_hole_filling(hoNDArray<float>& map, float hole, size_t max_size_of_holes, bool is_8_connected);
template EXPORTCMR void perform_hole_filling(hoNDArray<double>& map, double hole, size_t max_size_of_holes, bool is_8_connected);

// ---------------------------------------------------------------------

template <typename T>
void perform_hole_filling(Gadgetron::hoNDImageContainer2D< hoMRImage<T, 2> >& maps, T hole, size_t max_size_of_holes, bool is_8_connected) 
{
    try
    {
        size_t row = maps.rows();
        std::vector<size_t> cols = maps.cols();

        size_t r, c;
        for (r = 0; r < row; r++)
        {
            for (c = 0; c < cols[r]; c++)
            {
                GADGET_CATCH_THROW( Gadgetron::perform_hole_filling(maps(r, c), hole, max_size_of_holes, is_8_connected) );
            }
        }
    }
    catch (...)
    {
        GADGET_THROW("Error happened in perform_hole_filling(hoNDImageContainer2D) ... ");
    }
}

template EXPORTCMR void perform_hole_filling(Gadgetron::hoNDImageContainer2D< hoMRImage<float, 2> >& maps, float hole, size_t max_size_of_holes, bool is_8_connected);
template EXPORTCMR void perform_hole_filling(Gadgetron::hoNDImageContainer2D< hoMRImage<double, 2> >& maps, double hole, size_t max_size_of_holes, bool is_8_connected);

// ---------------------------------------------------------------------

template <typename T> 
CmrParametricMapping<T>::CmrParametricMapping()
{
    fill_holes_in_maps_ = true;
    max_size_of_holes_ = 20;
    hole_marking_value_ = 0;

    compute_SD_maps_ = false;;

    max_iter_ = 50;
    max_fun_eval_ = 100;
    thres_fun_ = 1e-5;

    max_map_value_ = -1;
    min_map_value_ = 0;

    verbose_ = false;
    perform_timing_ = false;

    gt_timer_.set_timing_in_destruction(false);
    gt_timer_local_.set_timing_in_destruction(false);
}

template <typename T> 
CmrParametricMapping<T>::~CmrParametricMapping()
{
}

template <typename T>
void CmrParametricMapping<T>::perform_parametric_mapping()
{
    try
    {
        size_t n, s, slc;

        size_t RO = data_.get_size(0);
        size_t E1 = data_.get_size(1);
        size_t N = data_.get_size(2);
        size_t S = data_.get_size(3);
        size_t SLC = data_.get_size(4);

        GADGET_CHECK_THROW(N >= ti_.size());

        size_t num_ti = ti_.size();
        size_t NUM = this->get_num_of_paras();

        map_.create(RO, E1, S, SLC);
        Gadgetron::clear(map_);

        para_.create(RO, E1, NUM, S, SLC);
        Gadgetron::clear(para_);

        if (this->compute_SD_maps_)
        {
            sd_map_.create(RO, E1, S, SLC);
            Gadgetron::clear(sd_map_);

            sd_para_.create(RO, E1, NUM, S, SLC);
            Gadgetron::clear(sd_para_);
        }

        if (!debug_folder_.empty()) gt_exporter_.export_array(this->data_, debug_folder_ + "CmrParametricMapping_data");

        if (this->verbose_)
        {
            this->data_.print(std::cout);

            GDEBUG_STREAM("Time points of mapping : " << ti_.size());
            for (n = 0; n < ti_.size(); n++)
            {
                GDEBUG_STREAM("Time point " << n << " - " << ti_[n]);
            }
        }

        T* pMask = NULL;
        if (mask_for_mapping_.get_size(0) == RO && mask_for_mapping_.get_size(1) == E1 && mask_for_mapping_.get_size(2) == SLC)
        {
            pMask = mask_for_mapping_.get_data_ptr();

            if (!debug_folder_.empty()) gt_exporter_.export_array(this->mask_for_mapping_, debug_folder_ + "CmrParametricMapping_mask_for_mapping");
        }

        if (this->perform_timing_) { gt_timer_.start("perform pixel-wise mapping ... "); }

        long long ro, e1;

        for (slc = 0; slc < SLC; slc++)
        {
            for (s = 0; s < S; s++)
            {
                T* pData = &data_(0, 0, 0, s, slc);

                T* pMap = &map_(0, 0, s, slc);
                T* pPara = &para_(0, 0, 0, s, slc);

                T* pMapSD = NULL;
                T* pParaSD = NULL;
                if (this->compute_SD_maps_)
                {
                    pMapSD = &sd_map_(0, 0, s, slc);
                    pParaSD = &sd_para_(0, 0, 0, s, slc);
                }

                T* pMaskCurr = NULL;
                if (pMask != NULL)
                {
                    pMaskCurr = pMask + s*RO*E1 + slc*S*RO*E1;
                }

#pragma omp parallel private(e1, ro, n) shared(RO, E1, pMask, pMaskCurr, pData, pMap, pMapSD, pPara, pParaSD, num_ti, NUM)
                {
                    std::vector<T> yi(num_ti, 0);
                    std::vector<T> guess(NUM + 1, 0);
                    std::vector<T> bi(NUM + 1, 0);
                    std::vector<T> sd(NUM + 1, 0);

                    T map_v(0), map_sd(0);

#pragma omp for 
                    for (e1 = 0; e1 < E1; e1++)
                    {
                        for (ro = 0; ro < RO; ro++)
                        {
                            long long offset = ro + e1*RO;

                            if (pMask != NULL)
                            {
                                if (pMaskCurr[offset] <= 0)
                                    continue;
                            }

                            // get data vector
                            for (n = 0; n < num_ti; n++)
                            {
                                yi[n] = pData[offset + n*RO*E1];
                            }

                            // estimate initial para
                            this->get_initial_guess(ti_, yi, guess);

                            // perform mapping
                            this->compute_map(ti_, yi, guess, bi, map_v);

                            pMap[offset] = map_v;
                            for (n = 0; n < NUM; n++)
                            {
                                pPara[offset + n*RO*E1] = bi[n];
                            }

                            // compute SD if needed
                            if (this->compute_SD_maps_)
                            {
                                try
                                {
                                    this->compute_sd(ti_, yi, bi, sd, map_sd);
                                }
                                catch(...)
                                {
                                    for (n = 0; n < NUM; n++)
                                    {
                                        sd[n] = 0;
                                    }

                                    map_sd = 0;
                                }

                                pMapSD[offset] = map_sd;
                                for (n = 0; n < NUM; n++)
                                {
                                    pParaSD[offset + n*RO*E1] = sd[n];
                                }
                            }
                        }
                    }
                } // openmp
            }
        }

        if (this->perform_timing_) { gt_timer_.stop(); }

        if (!debug_folder_.empty()) gt_exporter_.export_array(this->map_, debug_folder_ + "CmrParametricMapping_map");
        if (!debug_folder_.empty()) gt_exporter_.export_array(this->para_, debug_folder_ + "CmrParametricMapping_para");

        if (this->compute_SD_maps_)
        {
            if (!debug_folder_.empty()) gt_exporter_.export_array(this->sd_map_, debug_folder_ + "CmrParametricMapping_sd_map");
            if (!debug_folder_.empty()) gt_exporter_.export_array(this->sd_para_, debug_folder_ + "CmrParametricMapping_sd_para");
        }

        if (this->fill_holes_in_maps_)
        {
            if (this->perform_timing_) { gt_timer_.start("perform hole filling for pixel-wise mapping ... "); }

            bool is_8_connected = true;

            GADGET_CATCH_THROW( Gadgetron::perform_hole_filling(this->map_, (T)(this->hole_marking_value_), this->max_size_of_holes_, is_8_connected) );
            GADGET_CATCH_THROW( Gadgetron::perform_hole_filling(this->para_, (T)(this->hole_marking_value_), this->max_size_of_holes_, is_8_connected) );

            if (this->compute_SD_maps_)
            {
                GADGET_CATCH_THROW( Gadgetron::perform_hole_filling(this->sd_map_, (T)(this->hole_marking_value_), this->max_size_of_holes_, is_8_connected) );
                GADGET_CATCH_THROW( Gadgetron::perform_hole_filling(this->sd_para_, (T)(this->hole_marking_value_), this->max_size_of_holes_, is_8_connected) );
            }

            if (this->perform_timing_) { gt_timer_.stop(); }

            if (!debug_folder_.empty()) gt_exporter_.export_array(this->map_, debug_folder_ + "CmrParametricMapping_map_after_hole_filling");
            if (!debug_folder_.empty()) gt_exporter_.export_array(this->para_, debug_folder_ + "CmrParametricMapping_para_after_hole_filling");

            if (this->compute_SD_maps_)
            {
                if (!debug_folder_.empty()) gt_exporter_.export_array(this->sd_map_, debug_folder_ + "CmrParametricMapping_sd_map_after_hole_filling");
                if (!debug_folder_.empty()) gt_exporter_.export_array(this->sd_para_, debug_folder_ + "CmrParametricMapping_sd_para_after_hole_filling");
            }
        }

        size_t N_pts = this->map_.get_number_of_elements();

        size_t ii;
        for (ii = 0; ii<N_pts; ii++)
        {
            T v = this->map_(ii);
            if (v >= this->max_map_value_) v = this->max_map_value_;
            if (v <= this->min_map_value_) v = this->min_map_value_;

            this->map_(ii) = v;
        }

        if (this->compute_SD_maps_)
        {
            for (ii = 0; ii<N_pts; ii++)
            {
                T v = this->sd_map_(ii);
                if (v <= 0) v = 0;

                this->sd_map_(ii) = v;
            }
        }
    }
    catch (...)
    {
        GADGET_THROW("Error happened in CmrParametricMapping<T>::perform_parametric_mapping ... ");
    }
}

template <typename T>
void CmrParametricMapping<T>::get_initial_guess(const std::vector<T>& ti, const std::vector<T>& yi, std::vector<T>& guess)
{
    guess.clear();
    guess.resize(this->get_num_of_paras(), 0);
}

template <typename T>
void CmrParametricMapping<T>::compute_map(const std::vector<T>& ti, const std::vector<T>& yi, const std::vector<T>& guess, std::vector<T>& bi, T& map_v)
{
    bi = guess;
    map_v = 0;
}

template <typename T>
void CmrParametricMapping<T>::compute_sd(const std::vector<T>& ti, const std::vector<T>& yi, const std::vector<T>& bi, std::vector<T>& sd, T& map_sd)
{
    sd.clear();
    sd.resize(bi.size(), 0);

    map_sd = 0;
}

template <typename T>
void CmrParametricMapping<T>::compute_sd_impl(const VectorType& ti, const VectorType& yi, const VectorType& bi, const VectorType& res, const hoNDArray<T>& grad, VectorType& sd)
{
    try
    {
        sd.clear();
        sd.resize(bi.size(), 0);

        size_t num = ti.size();
        size_t N = this->get_num_of_paras();

        size_t n;

        // rank
        size_t rank = num - (N - 1);

        VectorType res_sorted(res);
        std::sort(res_sorted.begin(), res_sorted.end());

        //hoNDArray<T> res_array(res.size(), const_cast<T*>(&res[0]));
        //hoNDArray<T> res_array
        //Gadgetron::sort(res_array, res_array, true);

        VectorType res_trunc(rank);
        memcpy(&res_trunc[0], &res_sorted[0] + N - 1, sizeof(T)*res_trunc.size());

        T std;
        if (res_trunc.size() % 2 == 0)
        {
            std = (res_trunc[res_trunc.size() / 2] + res_trunc[res_trunc.size() / 2 - 1]) / 2 / 0.6745;
        }
        else
        {
            std = res_trunc[std::floor(res_trunc.size() / 2)] / 0.6745;
        }

        if (std::abs(std) < FLT_EPSILON)
        {
            return; // deviation is too small to compute sd
        }

        hoNDArray<T> hessian(N, N);
        Gadgetron::clear(hessian);

        size_t r, c;
        for (n = 0; n<num; n++)
        {
            for (c = 0; c < N; c++)
            {
                for (r = 0; r < N; r++)
                {
                    hessian(r, c) += grad(r, n) * grad(c, n);
                }
            }
        }

        Gadgetron::scal((T)(1.0 / (std*std)), hessian);

        // take inversion
        Gadgetron::invert(hessian);

        for (n = 0; n < N; n++)
        {
            sd[n] = std::sqrt(hessian(n, n));
        }
    }
    catch (...)
    {
        GADGET_THROW("Exceptions happened in CmrParametricMapping<T>::compute_sd_impl(...) ... ");
    }
}

template <typename T>
size_t CmrParametricMapping<T>::get_num_of_paras() const
{
    return 1;
}

// ------------------------------------------------------------
// Instantiation
// ------------------------------------------------------------

template class EXPORTCMR CmrParametricMapping< float >;

}
