
/** \file   mri_core_utility.cpp
    \brief  Implementation useful utility functionalities for 2D and 3D MRI parallel imaging
    \author Hui Xue
*/

#include "mri_core_utility.h"
#include "hoNDArray_elemwise.h"
#include "mri_core_kspace_filter.h"
#include <ctime>

namespace Gadgetron
{
    template <typename T>
    std::tuple<hoNDArray<bool> > detect_readout_sampling_status(const hoNDArray<T>& data)
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

            hoNDArray<bool> sampled;
            sampled.create(E1, E2, N, S, SLC);

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
                                    sampled(e1, e2, n, s, slc) = true;
                                }
                                else
                                {
                                    sampled(e1, e2, n, s, slc) =false;
                                }
                            }
                        }
                    }
                }
            }

            return std::make_tuple(std::move(sampled));
        }
        catch (...)
        {
            GADGET_THROW("Errors in detect_readout_sampling_status(...) ... ");
        }
    }

    template EXPORTMRICORE std::tuple<hoNDArray<bool> > detect_readout_sampling_status(const hoNDArray< float >& data);
    template EXPORTMRICORE std::tuple<hoNDArray<bool> > detect_readout_sampling_status(const hoNDArray< double >& data);
    template EXPORTMRICORE std::tuple<hoNDArray<bool> > detect_readout_sampling_status(const hoNDArray< std::complex<float> >& data);
    template EXPORTMRICORE std::tuple<hoNDArray<bool> > detect_readout_sampling_status(const hoNDArray< std::complex<double> >& data);

    // ------------------------------------------------------------------------

    template <typename T>
    std::tuple<size_t, size_t> detect_sampled_region_E1(const hoNDArray<T>& data)
    {
        try
        {
            hoNDArray<bool> sampled;
            sampled = std::get<0>(Gadgetron::detect_readout_sampling_status(data));

            size_t RO = data.get_size(0);
            size_t E1 = data.get_size(1);
            size_t E2 = data.get_size(2);
            size_t CHA = data.get_size(3);
            size_t N = data.get_size(4);
            size_t S = data.get_size(5);
            size_t SLC = data.get_size(6);

            size_t start_E1, end_E1;
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
                                if (sampled(e1, e2, n, s, slc)==true)
                                {
                                    if (e1 < start_E1) start_E1 = e1;
                                    if (e1 > end_E1) end_E1 = e1;
                                }
                            }
                        }
                    }
                }
            }

            return std::make_tuple(start_E1, end_E1);
        }
        catch (...)
        {
            GADGET_THROW("Errors in detect_sampled_region_E1(...) ... ");
        }
    }

    template EXPORTMRICORE std::tuple<size_t, size_t> detect_sampled_region_E1(const hoNDArray<float>& data);
    template EXPORTMRICORE std::tuple<size_t, size_t> detect_sampled_region_E1(const hoNDArray<double>& data);
    template EXPORTMRICORE std::tuple<size_t, size_t> detect_sampled_region_E1(const hoNDArray< std::complex<float> >& data);
    template EXPORTMRICORE std::tuple<size_t, size_t> detect_sampled_region_E1(const hoNDArray< std::complex<double> >& data);

    // ------------------------------------------------------------------------

    template <typename T>
    std::tuple<size_t, size_t> detect_sampled_region_E2(const hoNDArray<T>& data)
    {
        try
        {
            hoNDArray<bool> sampled;
            sampled = std::get<0>(Gadgetron::detect_readout_sampling_status(data));

            size_t RO = data.get_size(0);
            size_t E1 = data.get_size(1);
            size_t E2 = data.get_size(2);
            size_t CHA = data.get_size(3);
            size_t N = data.get_size(4);
            size_t S = data.get_size(5);
            size_t SLC = data.get_size(6);

            size_t start_E2, end_E2;
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
                                if (sampled(e1, e2, n, s, slc)==true)
                                {
                                    if (e2 < start_E2) start_E2 = e2;
                                    if (e2 > end_E2) end_E2 = e2;
                                }
                            }
                        }
                    }
                }
            }

            return std::make_tuple(start_E2, end_E2);
        }
        catch (...)
        {
            GADGET_THROW("Errors in detect_sampled_region_E2(...) ... ");
        }
    }

    template EXPORTMRICORE std::tuple<size_t, size_t> detect_sampled_region_E2(const hoNDArray<float>& data);
    template EXPORTMRICORE std::tuple<size_t, size_t> detect_sampled_region_E2(const hoNDArray<double>& data);
    template EXPORTMRICORE std::tuple<size_t, size_t> detect_sampled_region_E2(const hoNDArray< std::complex<float> >& data);
    template EXPORTMRICORE std::tuple<size_t, size_t> detect_sampled_region_E2(const hoNDArray< std::complex<double> >& data);

    // ------------------------------------------------------------------------
}
