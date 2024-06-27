
#include "hoWaveletOperator.h"

namespace Gadgetron
{
    template <typename T>
    hoWaveletOperator<T>::hoWaveletOperator(const std::vector<size_t>& dims) : input_in_kspace_(false), no_null_space_(true), num_of_wav_levels_(1), with_approx_coeff_(false), proximity_across_cha_(false), BaseClass(dims)
    {
        p_active_wav_ = &harr_wav_;

        //gt_timer1_.set_timing_in_destruction(false);
        //gt_timer2_.set_timing_in_destruction(false);
        //gt_timer3_.set_timing_in_destruction(false);
    }

    template <typename T>
    hoWaveletOperator<T>::~hoWaveletOperator()
    {
    }

    template <typename T>
    void hoWaveletOperator<T>::select_wavelet(const std::string& wav_name)
    {
        if (wav_name == "db2" || wav_name == "db3" || wav_name == "db4" || wav_name == "db5")
        {
            redundant_wav_.compute_wavelet_filter(wav_name);
            p_active_wav_ = &redundant_wav_;
        }
        else
        {
            p_active_wav_ = &harr_wav_;
        }
    }

    template <typename T>
    void hoWaveletOperator<T>::restore_acquired_kspace(const ARRAY_TYPE& acquired, ARRAY_TYPE& y)
    {
        try
        {
            GADGET_CHECK_THROW(acquired.get_number_of_elements() == y.get_number_of_elements());

            size_t N = acquired.get_number_of_elements();

            const T* pA = acquired.get_data_ptr();
            T* pY = y.get_data_ptr();

            int n(0);
#pragma omp parallel for default(none) private(n) shared(N, pA, pY)
            for (n = 0; n<(int)N; n++)
            {
                if (std::abs(pA[n]) > 0)
                {
                    pY[n] = pA[n];
                }
            }
        }
        catch (...)
        {
            GADGET_THROW("Errors happened in hoWaveletOperator<T>::restore_acquired_kspace(...) ... ");
        }
    }

    template <typename T>
    void hoWaveletOperator<T>::restore_acquired_kspace(ARRAY_TYPE& y)
    {
        this->restore_acquired_kspace(acquired_points_, y);
    }

    template <typename T>
    void hoWaveletOperator<T>::set_acquired_points(ARRAY_TYPE& kspace)
    {
        try
        {
            std::vector<size_t> dim;
            kspace.get_dimensions(dim);
            acquired_points_.create(dim, kspace.begin());

            acquired_points_indicator_.create(kspace.dimensions());
            Gadgetron::clear(acquired_points_indicator_);

            unacquired_points_indicator_.create(kspace.dimensions());
            Gadgetron::clear(unacquired_points_indicator_);

            size_t N = kspace.get_number_of_elements();

            long long ii(0);

#pragma omp parallel for default(shared) private(ii) shared(N, kspace)
            for (ii = 0; ii<(long long)N; ii++)
            {
                if (std::abs(kspace(ii)) < DBL_EPSILON)
                {
                    unacquired_points_indicator_(ii) = T(1.0);
                }
                else
                {
                    acquired_points_indicator_(ii) = T(1.0);
                }
            }

            // allocate the helper memory
            kspace_.create(kspace.dimensions());
            complexIm_.create(kspace.dimensions());
        }
        catch (...)
        {
            GADGET_THROW("Errors in hoWaveletOperator<T>::set_acquired_points(...) ... ");
        }
    }

    // ------------------------------------------------------------
    // Instantiation
    // ------------------------------------------------------------

    template class EXPORTCPUOPERATOR hoWaveletOperator< float >;
    template class EXPORTCPUOPERATOR hoWaveletOperator< double >;

    template class EXPORTCPUOPERATOR hoWaveletOperator< std::complex<float> >;
    template class EXPORTCPUOPERATOR hoWaveletOperator< std::complex<double> >;

}
