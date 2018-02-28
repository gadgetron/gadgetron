/** \file   cmr_time_stamp.h
    \brief  Implement functionalities to handle cardiac time stamps
    \author Hui Xue
*/

#include "cmr_time_stamp.h"
#include "hoNDArray_reductions.h"
#include "hoNDArray_elemwise.h"
#include <boost/math/special_functions/sign.hpp>

namespace Gadgetron { 

    template <typename T> 
    void simple_line_fit(const std::vector<T>& x, const std::vector<T>& y, T& a, T& b)
    {
        try
        {
            size_t num = x.size();

            if(num<2)
            {
                a = 0;
                b = 0;
                return;
            }

            T sx(0), sy(0);

            size_t n;
            for (n=0; n<num; n++)
            {
                sx += x[n];
                sy += y[n];
            }

            T mx = sx / (T)(num);
            T syy = 0;
            b = 0;
            for (n=0; n<num; n++)
            {
                T v = (x[n] - mx);
                syy += v*v;
                b += v*y[n];
            }

            syy = (std::abs(syy) > 0 ? syy : boost::math::sign(syy)*FLT_EPSILON);
            b /= syy;
            a = (sy - sx*b) / (T)(num);
        }
        catch(...)
        {
            GADGET_THROW("Exceptions happened in simple_line_fit ... ");
        }
    }

    void correct_time_stamp_with_fitting(hoNDArray<float>& time_stamp, size_t startE1, size_t endE1)
    {
        try
        {
            size_t E1 = time_stamp.get_size(0);
            size_t N = time_stamp.get_size(1);
            size_t rE1 = endE1 - startE1 + 1;

            size_t e1, n;

            size_t num_acq_read_outs = 0;
            for ( n=0; n<N; n++ )
            {
                for ( e1=0; e1<E1; e1++ )
                {
                    if ( time_stamp(e1, n) > 0 )
                    {
                        num_acq_read_outs++;
                    }
                }
            }

            GDEBUG_STREAM(" Number of acquired lines : " << num_acq_read_outs);

            float a, b; // y = a + b*x
            {
                std::vector<float> x(num_acq_read_outs), y(num_acq_read_outs);

                size_t ind = 0;
                for ( n=0; n<N; n++ )
                {
                    for ( e1=startE1; e1<=endE1; e1++ )
                    {
                        float acq_time = time_stamp(e1, n);
                        if ( acq_time > 0 )
                        {
                            x[ind] = (float)(e1-startE1 + n*rE1);
                            y[ind] = acq_time;
                            ind++;
                        }
                    }
                }

                Gadgetron::simple_line_fit(x, y, a, b);
            }

            for ( n=0; n<N; n++ )
            {
                for ( e1=startE1; e1<=endE1; e1++ )
                {
                    float x_v = (float)(e1-startE1 + n*rE1);
                    time_stamp(e1, n) = a + b*x_v;
                }
            }
        }
        catch(...)
        {
            GADGET_THROW("Exceptions happened in correct_time_stamp_with_fitting(...) ... ");
        }
    }

    void detect_heart_beat_with_time_stamp(hoNDArray<float>& cpt_time_stamp, hoNDArray<int>& ind_hb, 
                                        std::vector<size_t>& start_e1_hb, std::vector<size_t>& end_e1_hb, 
                                        std::vector<size_t>& start_n_hb, std::vector<size_t>& end_n_hb )
    {
        try
        {
            size_t E1 = cpt_time_stamp.get_size(0);
            size_t N = cpt_time_stamp.get_size(1);

            size_t e1, n, ind, ii;

            size_t num_acq_read_outs = 0;
            for ( n=0; n<N; n++ )
            {
                for ( e1=0; e1<E1; e1++ )
                {
                    if ( cpt_time_stamp(e1, n) >= 0 )
                    {
                        num_acq_read_outs++;
                    }
                }
            }

            ind_hb.create(E1, N);
            Gadgetron::clear(ind_hb);

            // --------------------------------------------------------
            // cpt time stamps
            // --------------------------------------------------------

            std::vector<float> acquired_cpt(num_acq_read_outs);
            std::vector<size_t> ind_acquired_cpt(num_acq_read_outs);

            ind = 0;
            for ( n=0; n<N; n++ )
            {
                for ( e1=0; e1<E1; e1++ )
                {
                    if ( cpt_time_stamp(e1, n) > -1 )
                    {
                        acquired_cpt[ind] = cpt_time_stamp(e1, n);
                        ind_acquired_cpt[ind] = e1 + n*E1;
                        ind++;
                    }
                }
            }

            // --------------------------------------------------------
            // find the number of heart beats
            // --------------------------------------------------------
            size_t numOfHB = 0;

            // store the line indexes for every heart beat
            std::vector<size_t> ind_HB_start, ind_HB_end;
            ind_HB_start.push_back(0);

            for ( ind=1; ind<num_acq_read_outs; ind++ )
            {
                if ( acquired_cpt[ind] < acquired_cpt[ind-1] )
                {
                    // find a new heart beat
                    numOfHB++;

                    size_t end_ind_prev_HB = ind_acquired_cpt[ind-1];
                    size_t start_ind_curr_HB = ind_acquired_cpt[ind];

                    // if there is a gap between end and start ind, fill the gap
                    if ( end_ind_prev_HB+1 != start_ind_curr_HB )
                    {
                        long long gap = start_ind_curr_HB - end_ind_prev_HB - 1;
                        if ( gap % 2 == 0 )
                        {
                            end_ind_prev_HB += gap;
                        }
                        else
                        {
                            end_ind_prev_HB += gap;
                        }

                        if ( end_ind_prev_HB+1 != start_ind_curr_HB )
                        {
                            GWARN_STREAM("end_ind_prev_HB+1 ~= start_ind_curr_HB : " << end_ind_prev_HB << " " << start_ind_curr_HB);
                        }
                    }

                    ind_HB_end.push_back( end_ind_prev_HB );
                    ind_HB_start.push_back( start_ind_curr_HB );
                }
            }

            ind_HB_end.push_back( E1*N-1 );
            numOfHB = ind_HB_end.size();

            // --------------------------------------------------------
            // fill the start and end indexes
            // --------------------------------------------------------
            start_e1_hb.resize(numOfHB, 0);
            end_e1_hb.resize(numOfHB, 0);

            start_n_hb.resize(numOfHB, 0);
            end_n_hb.resize(numOfHB, 0);

            std::vector<size_t> start, end;
            for ( ii=0; ii<numOfHB; ii++ )
            {
                start_n_hb[ii] = ind_HB_start[ii] / E1;
                start_e1_hb[ii] = ind_HB_start[ii] - start_n_hb[ii] * E1;

                end_n_hb[ii] = ind_HB_end[ii] / E1;
                end_e1_hb[ii] = ind_HB_end[ii] - end_n_hb[ii]*E1;

                for (ind=ind_HB_start[ii]; ind<=ind_HB_end[ii]; ind++)
                {
                    ind_hb(ind) = ii;
                }
            }
        }
        catch(...)
        {
            GADGET_THROW("Exceptions happened in detect_heart_beat_with_time_stamp(...) ... ");
        }
    }

    void correct_heart_beat_time_stamp_with_fitting(hoNDArray<float>& cpt_time_stamp, hoNDArray<int>& ind_hb, size_t startE1, size_t endE1, 
                                                const std::vector<size_t>& start_e1_hb, const std::vector<size_t>& end_e1_hb, 
                                                const std::vector<size_t>& start_n_hb, const std::vector<size_t>& end_n_hb )
    {
        try
        {
            size_t E1 = cpt_time_stamp.get_size(0);
            size_t N = cpt_time_stamp.get_size(1);
            size_t rE1 = endE1-startE1+1;

            size_t e1, n, ind, ii;

            size_t num_acq_read_outs = 0;
            for ( n=0; n<N; n++ )
            {
                for ( e1=0; e1<E1; e1++ )
                {
                    if ( cpt_time_stamp(e1, n) >= 0 )
                    {
                        num_acq_read_outs++;
                    }
                }
            }

            size_t numOfHB = start_e1_hb.size();

            std::vector<size_t> ind_HB_start(numOfHB);
            std::vector<size_t> ind_HB_end(numOfHB);

            for ( ind=0; ind<numOfHB; ind++ )
            {
                ind_HB_start[ind] = start_e1_hb[ind] + start_n_hb[ind] * E1;
                ind_HB_end[ind] = end_e1_hb[ind] + end_n_hb[ind] * E1;
            }

            // --------------------------------------------------------
            // fit a line to every heart beat
            // --------------------------------------------------------
            float a, b;
            std::vector<float> A(numOfHB, 0.0f), B(numOfHB, 0.0f);
            for ( ind=0; ind<numOfHB; ind++ )
            {
                std::vector<float> x, y;

                size_t cpt;
                for ( cpt=ind_HB_start[ind]; cpt<=ind_HB_end[ind]; cpt++ )
                {
                    size_t n = cpt / E1;
                    size_t e1 = cpt - n*E1;

                    if(e1>=startE1 && e1<=endE1)
                    {
                        if ( cpt_time_stamp[cpt] > -1 )
                        {
                            size_t x_ind = (e1-startE1) + n*rE1;
                            x.push_back( (float)x_ind );
                            y.push_back(cpt_time_stamp[cpt]);
                        }
                    }
                }

                if ( !x.empty() )
                {
                    Gadgetron::simple_line_fit(x, y, a, b);
                    A[ind] = a;
                    B[ind] = b;
                }
            }

            // --------------------------------------------------------
            // compute cpt time stamp for every line
            // --------------------------------------------------------
            size_t num = cpt_time_stamp.get_number_of_elements();
            for ( ind=0; ind<num; ind++ )
            {
                n = ind / E1;
                e1 = ind - n*E1;

                if(e1>=startE1 && e1<=endE1)
                {
                    // find to which heart beat this line belongs
                    bool foundHB = false;
                    for ( ii=0; ii<numOfHB; ii++ )
                    {
                        size_t startHB = ind_HB_start[ii];
                        size_t endHB = ind_HB_end[ii];

                        if ( ii==0 && ind<=startHB )
                        {
                            foundHB = true;
                            break;
                        }

                        if ( ii==numOfHB-1 && ind>=endHB )
                        {
                            foundHB = true;
                            break;
                        }

                        if ( ind>=startHB && ind<=endHB )
                        {
                            foundHB = true;
                            break;
                        }
                    }

                    // if cannot find a heart beat, this kspace line will not be used
                    if ( foundHB && (std::abs(B[ii])>0) )
                    {
                        ind_hb(e1, n) = ii;

                        size_t x_ind = (e1-startE1) + n*rE1;
                        cpt_time_stamp(e1, n) = (float)(A[ii] + B[ii]*x_ind);
                    }
                    else
                    {
                        ind_hb(e1, n) = -1;
                    }
                }
                else
                {
                    cpt_time_stamp(e1, n) = -1;
                    ind_hb(e1, n) = -1;
                }
            }
        }
        catch(...)
        {
            GADGET_THROW("Exceptions happened in correct_heart_beat_time_stamp_with_fitting(...) ... ");
        }
    }

    void compute_phase_time_stamp(const hoNDArray<float>& time_stamp, const hoNDArray<float>& cpt_time_stamp, size_t startE1, size_t endE1, 
        hoNDArray<float>& phs_time_stamp, hoNDArray<float>& phs_cpt_time_stamp)
    {
        try
        {
            size_t E1 = time_stamp.get_size(0);
            size_t N = time_stamp.get_size(1);
            size_t rE1 = endE1 - startE1 + 1;

            size_t e1, n;

            for ( n=0; n<N; n++ )
            {
                // phase time stamp as the mean of all aquired lines
                size_t num = 0;
                float tt = 0.0f;
                for ( e1=startE1; e1<=endE1; e1++ )
                {
                    if(time_stamp(e1, n)>0)
                    {
                        tt += time_stamp(e1, n);
                        num++;
                    }
                }
                phs_time_stamp(n, 0) = tt/((num>0) ? num : 1);

                //// phase cpt time as the median of all acquired lines
                //std::vector<float> cpt_buf(rE1, 0);
                //for ( e1=startE1; e1<=endE1; e1++ )
                //{
                //    if(cpt_time_stamp(e1, n)>=0)
                //        cpt_buf[e1-startE1] = cpt_time_stamp(e1, n);
                //}

                //std::sort(cpt_buf.begin(), cpt_buf.end());
                //phs_cpt_time_stamp(n, 0) = cpt_buf[E1/2-startE1];

                // phase cpt time as the cpt time of center line
                phs_cpt_time_stamp(n, 0) = cpt_time_stamp(E1/2, n);
            }
        }
        catch(...)
        {
            GADGET_THROW("Exceptions happened in compute_phase_time_stamp(...) ... ");
        }
    }

    template <typename T> 
    void resample_cardiac_phase_cmr_array(const hoNDArray<T>& data, size_t output_N, hoNDArray<T>& res, size_t spline_degree)
    {
        try
        {
            std::vector<size_t> dims;
            data.get_dimensions(dims);

            size_t num_of_dims = dims.size();

            GADGET_CHECK_THROW(num_of_dims >=2);

            size_t N = dims[num_of_dims -1];

            GADGET_CHECK_THROW(N>2);

            std::vector<size_t> dims_res(dims);
            dims_res[num_of_dims - 1] = output_N;

            res.create(dims_res);
            Gadgetron::clear(res);

            long long num = data.get_number_of_elements()/N;

            long long n;

#pragma omp parallel default(none) shared(N, output_N, num, data, res, spline_degree, num_of_dims) private (n)
            {
                hoNDArray<T> data_in(N);
                hoNDArray<T> data_res(output_N);
                hoNDArray<T> coeff(N);
                hoNDBSpline< T, 1 > interp;

                float dx = (float(N)-1)/(float(output_N)-1);

                std::vector<size_t> ind;

                #pragma omp for 
                for (n=0; n<num; n++)
                {
                    ind = data.calculate_index(n);

                    size_t ii;
                    for(ii=0; ii<N; ii++)
                    {
                        ind[num_of_dims - 1] = ii;
                        data_in(ii) = data(ind);
                    }

                    interp.computeBSplineCoefficients(data_in, spline_degree, coeff);

                    for (ii = 0; ii<output_N; ii++)
                    {
                        ind[num_of_dims - 1] = ii;
                        data_res(ii) = interp.evaluateBSpline(coeff.begin(), N, spline_degree, 0, ii*dx);
                        res(ind) = data_res(ii);
                    }
                }
            }
        }
        catch(...)
        {
            GADGET_THROW("Exceptions happened in resample_cardiac_phase_cmr_array(...) ... ");
        }
    }

    template EXPORTCMR void resample_cardiac_phase_cmr_array(const hoNDArray<float>& data, size_t output_N, hoNDArray<float>& res, size_t spline_degree);
    template EXPORTCMR void resample_cardiac_phase_cmr_array(const hoNDArray<double>& data, size_t output_N, hoNDArray<double>& res, size_t spline_degree);
    template EXPORTCMR void resample_cardiac_phase_cmr_array(const hoNDArray< std::complex<float> >& data, size_t output_N, hoNDArray< std::complex<float> >& res, size_t spline_degree);
    template EXPORTCMR void resample_cardiac_phase_cmr_array(const hoNDArray< std::complex<double> >& data, size_t output_N, hoNDArray< std::complex<double> >& res, size_t spline_degree);
}
