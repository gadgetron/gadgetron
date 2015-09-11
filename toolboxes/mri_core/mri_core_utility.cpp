
/** \file   mri_core_utility.cpp
    \brief  Implementation useful utility functionalities for 2D and 3D MRI parallel imaging
    \author Hui Xue
*/

#include "mri_core_utility.h"
#include "hoNDArray_elemwise.h"
#include "mri_core_kspace_filter.h"

namespace Gadgetron
{
    template <typename T>
    void detect_readout_sampling_status(const hoNDArray<T>& data, hoNDArray<float>& sampled)
    {
        try
        {
            typedef typename realType<T>::Type value_type;

            size_t RO = data.get_size(0);
            size_t E1 = data.get_size(1);
            size_t E2 = data.get_size(2);
            size_t CHA = data.get_size(3);
            size_t N = data.get_size(4);
            size_t S = data.get_size(5);
            size_t SLC = data.get_size(6);

            sampled.create(E1, E2, N, S, SLC);
            Gadgetron::clear(sampled);

            size_t slc, s, n, e2, e1;

            for (slc = 0; slc < SLC; slc++)
            {
                for (s = 0; s < S; s++)
                {
                    for (n = 0; n < N; n++)
                    {
                        for (e2 = 0; e2 < E2; e2++)
                        {
                            for (e1 = 0; e1 < E1; e1++)
                            {
                                value_type v = 0;
                                const T* pData = &(data(0, e1, e2, 0, n, s, slc));

                                size_t ro;
                                for (ro = 0; ro < RO; ro++)
                                {
                                    v += std::abs(pData[ro]);
                                }

                                if (v > 0)
                                {
                                    sampled(e1, e2, n, s, slc) = 1;
                                }
                            }
                        }
                    }
                }
            }
        }
        catch (...)
        {
            GADGET_THROW("Errors in detect_readout_sampling_status(...) ... ");
        }
    }

    template EXPORTMRICORE void detect_readout_sampling_status(const hoNDArray< float >& data, hoNDArray<float>& sampled);
    template EXPORTMRICORE void detect_readout_sampling_status(const hoNDArray< double >& data, hoNDArray<float>& sampled);
    template EXPORTMRICORE void detect_readout_sampling_status(const hoNDArray< std::complex<float> >& data, hoNDArray<float>& sampled);
    template EXPORTMRICORE void detect_readout_sampling_status(const hoNDArray< std::complex<double> >& data, hoNDArray<float>& sampled);

    // ------------------------------------------------------------------------

    template <typename T>
    void detect_sampled_region_E1(const hoNDArray<T>& data, size_t& start_E1, size_t& end_E1)
    {
        try
        {
            hoNDArray<float> sampled;
            Gadgetron::detect_readout_sampling_status(data, sampled);

            size_t RO = data.get_size(0);
            size_t E1 = data.get_size(1);
            size_t E2 = data.get_size(2);
            size_t CHA = data.get_size(3);
            size_t N = data.get_size(4);
            size_t S = data.get_size(5);
            size_t SLC = data.get_size(6);

            start_E1 = E1;
            end_E1 = 0;

            size_t e1, e2, n, s, slc;

            for (slc = 0; slc < SLC; slc++)
            {
                for (s = 0; s < S; s++)
                {
                    for (n = 0; n < N; n++)
                    {
                        for (e2 = 0; e2 < E2; e2++)
                        {
                            for (e1 = 0; e1 < E1; e1++)
                            {
                                if (sampled(e1, e2, n, s, slc)>0)
                                {
                                    if (e1 < start_E1) start_E1 = e1;
                                    if (e1 > end_E1) end_E1 = e1;
                                }
                            }
                        }
                    }
                }
            }

        }
        catch (...)
        {
            GADGET_THROW("Errors in detect_sampled_region_E1(...) ... ");
        }
    }

    template EXPORTMRICORE void detect_sampled_region_E1(const hoNDArray<float>& data, size_t& start_E1, size_t& end_E1);
    template EXPORTMRICORE void detect_sampled_region_E1(const hoNDArray<double>& data, size_t& start_E1, size_t& end_E1);
    template EXPORTMRICORE void detect_sampled_region_E1(const hoNDArray< std::complex<float> >& data, size_t& start_E1, size_t& end_E1);
    template EXPORTMRICORE void detect_sampled_region_E1(const hoNDArray< std::complex<double> >& data, size_t& start_E1, size_t& end_E1);

    // ------------------------------------------------------------------------

    template <typename T>
    void detect_sampled_region_E2(const hoNDArray<T>& data, size_t& start_E2, size_t& end_E2)
    {
        try
        {
            hoNDArray<float> sampled;
            Gadgetron::detect_readout_sampling_status(data, sampled);

            size_t RO = data.get_size(0);
            size_t E1 = data.get_size(1);
            size_t E2 = data.get_size(2);
            size_t CHA = data.get_size(3);
            size_t N = data.get_size(4);
            size_t S = data.get_size(5);
            size_t SLC = data.get_size(6);

            start_E2 = E2;
            end_E2 = 0;

            size_t e1, e2, n, s, slc;

            for (slc = 0; slc < SLC; slc++)
            {
                for (s = 0; s < S; s++)
                {
                    for (n = 0; n < N; n++)
                    {
                        for (e2 = 0; e2 < E2; e2++)
                        {
                            for (e1 = 0; e1 < E1; e1++)
                            {
                                if (sampled(e1, e2, n, s, slc)>0)
                                {
                                    if (e2 < start_E2) start_E2 = e2;
                                    if (e2 > end_E2) end_E2 = e2;
                                }
                            }
                        }
                    }
                }
            }

        }
        catch (...)
        {
            GADGET_THROW("Errors in detect_sampled_region_E2(...) ... ");
        }
    }

    template EXPORTMRICORE void detect_sampled_region_E2(const hoNDArray<float>& data, size_t& start_E2, size_t& end_E2);
    template EXPORTMRICORE void detect_sampled_region_E2(const hoNDArray<double>& data, size_t& start_E2, size_t& end_E2);
    template EXPORTMRICORE void detect_sampled_region_E2(const hoNDArray< std::complex<float> >& data, size_t& start_E2, size_t& end_E2);
    template EXPORTMRICORE void detect_sampled_region_E2(const hoNDArray< std::complex<double> >& data, size_t& start_E2, size_t& end_E2);

    // ------------------------------------------------------------------------

    template <typename T>
    void generate_ref_filter_for_coil_map(const hoNDArray<T>& ref, const SamplingLimit& lim_RO, const SamplingLimit& lim_E1, const SamplingLimit& lim_E2, hoNDArray<T>& filter_RO, hoNDArray<T>& filter_E1, hoNDArray<T>& filter_E2)
    {
        try
        {
            size_t RO = ref.get_size(0);
            size_t E1 = ref.get_size(1);
            size_t E2 = ref.get_size(2);

            Gadgetron::generate_symmetric_filter_ref(RO, lim_RO.min_, lim_RO.max_, filter_RO);
            Gadgetron::generate_symmetric_filter_ref(E1, lim_E1.min_, lim_E1.max_, filter_E1);

            if (E2 > 1)
            {
                Gadgetron::generate_symmetric_filter_ref(E2, lim_E2.min_, lim_E2.max_, filter_E2);
            }
        }
        catch (...)
        {
            GADGET_THROW("Errors in generate_ref_filter_for_coil_map(...) ... ");
        }
    }

    template EXPORTMRICORE void generate_ref_filter_for_coil_map(const hoNDArray< std::complex<float> >& ref, const SamplingLimit& lim_RO, const SamplingLimit& lim_E1, const SamplingLimit& lim_E2, hoNDArray< std::complex<float> >& filter_RO, hoNDArray< std::complex<float> >& filter_E1, hoNDArray< std::complex<float> >& filter_E2);
    template EXPORTMRICORE void generate_ref_filter_for_coil_map(const hoNDArray< std::complex<double> >& ref, const SamplingLimit& lim_RO, const SamplingLimit& lim_E1, const SamplingLimit& lim_E2, hoNDArray< std::complex<double> >& filter_RO, hoNDArray< std::complex<double> >& filter_E1, hoNDArray< std::complex<double> >& filter_E2);

    // ------------------------------------------------------------------------
}
