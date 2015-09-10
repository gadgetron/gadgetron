
/** \file   mri_core_utility.cpp
    \brief  Implementation useful utility functionalities for 2D and 3D MRI parallel imaging
    \author Hui Xue
*/

#include "mri_core_utility.h"
#include "hoNDArray_elemwise.h"
#include "hoNDFFT.h"
#include "hoNDArray_utils.h"

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

                                if ( v > 0 )
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
    void average_kspace_across_N(const hoNDArray<T>& data, hoNDArray<T>& averaged)
    {
        typedef typename realType<T>::Type value_type;

        try
        {
            size_t RO = data.get_size(0);
            size_t E1 = data.get_size(1);
            size_t E2 = data.get_size(2);
            size_t CHA = data.get_size(3);
            size_t N = data.get_size(4);
            size_t S = data.get_size(5);
            size_t SLC = data.get_size(6);

            hoNDArray<float> sampled;
            Gadgetron::detect_readout_sampling_status(data, sampled);

            Gadgetron::sum_over_dimension(data, averaged, 4);

            hoNDArray<float> sampled_across_N;
            Gadgetron::sum_over_dimension(sampled, sampled_across_N, 2);

            size_t ro, e1, e2, cha, s, slc;

            for (slc = 0; slc < SLC; slc++)
            {
                for (s = 0; s < S; s++)
                {
                    for (cha = 0; cha < CHA; cha++)
                    {
                        for (e2 = 0; e2 < E2; e2++)
                        {
                            for (e1 = 0; e1 < E1; e1++)
                            {
                                value_type freq = sampled_across_N(e1, e2, 0, s, slc);

                                if (freq > 1)
                                {
                                    value_type freq_reciprocal = (value_type)(1.0 / freq);

                                    T* pAve = &(averaged(0, e1, e2, cha, 0, s, slc));
                                    for (ro = 0; ro < RO; ro++)
                                    {
                                        pAve[ro] *= freq_reciprocal;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        catch (...)
        {
            GADGET_THROW("Errors in average_kspace_across_N(...) ... ");
        }
    }

    template EXPORTMRICORE void average_kspace_across_N(const hoNDArray< float >& data, hoNDArray< float >& averaged);
    template EXPORTMRICORE void average_kspace_across_N(const hoNDArray< double >& data, hoNDArray< double >& averaged);
    template EXPORTMRICORE void average_kspace_across_N(const hoNDArray< std::complex<float> >& data, hoNDArray< std::complex<float> >& averaged);
    template EXPORTMRICORE void average_kspace_across_N(const hoNDArray< std::complex<double> >& data, hoNDArray< std::complex<double> >& averaged);

    // ------------------------------------------------------------------------

    template <typename T> 
    void average_kspace_across_S(const hoNDArray<T>& data, hoNDArray<T>& averaged)
    {
        typedef typename realType<T>::Type value_type;

        try
        {
            size_t RO = data.get_size(0);
            size_t E1 = data.get_size(1);
            size_t E2 = data.get_size(2);
            size_t CHA = data.get_size(3);
            size_t N = data.get_size(4);
            size_t S = data.get_size(5);
            size_t SLC = data.get_size(6);

            hoNDArray<float> sampled;
            Gadgetron::detect_readout_sampling_status(data, sampled);

            Gadgetron::sum_over_dimension(data, averaged, 5);

            hoNDArray<float> sampled_across_S;
            Gadgetron::sum_over_dimension(sampled, sampled_across_S, 3);

            size_t ro, e1, e2, cha, n, slc;

            for (slc = 0; slc < SLC; slc++)
            {
                for (n = 0; n < N; n++)
                {
                    for (cha = 0; cha < CHA; cha++)
                    {
                        for (e2 = 0; e2 < E2; e2++)
                        {
                            for (e1 = 0; e1 < E1; e1++)
                            {
                                value_type freq = sampled_across_S(e1, e2, n, 0, slc);

                                if (freq > 1)
                                {
                                    value_type freq_reciprocal = (value_type)(1.0 / freq);

                                    T* pAve = &(averaged(0, e1, e2, cha, n, 0, slc));
                                    for (ro = 0; ro < RO; ro++)
                                    {
                                        pAve[ro] *= freq_reciprocal;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        catch (...)
        {
            GADGET_THROW("Errors in average_kspace_across_S(...) ... ");
        }
    }

    template EXPORTMRICORE void average_kspace_across_S(const hoNDArray< float >& data, hoNDArray< float >& averaged);
    template EXPORTMRICORE void average_kspace_across_S(const hoNDArray< double >& data, hoNDArray< double >& averaged);
    template EXPORTMRICORE void average_kspace_across_S(const hoNDArray< std::complex<float> >& data, hoNDArray< std::complex<float> >& averaged);
    template EXPORTMRICORE void average_kspace_across_S(const hoNDArray< std::complex<double> >& data, hoNDArray< std::complex<double> >& averaged);

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

    template <typename T> 
    void zero_pad_resize(const hoNDArray<T>& complexIm, size_t sizeRO, size_t sizeE1, size_t sizeE2, hoNDArray<T>& complexImResized)
    {
        try
        {
            std::vector<size_t> dim;
            complexIm.get_dimensions(dim);

            size_t RO = complexIm.get_size(0);
            size_t E1 = complexIm.get_size(1);
            size_t E2 = complexIm.get_size(2);

            GADGET_CHECK_THROW(sizeRO >= RO);
            GADGET_CHECK_THROW(sizeE1 >= E1);

            if (sizeE2 > 1)
            {
                if (RO == sizeRO && E1 == sizeE1 && E2 == sizeE2)
                {
                    complexImResized = complexIm;
                    return;
                }

                if (complexImResized.get_size(0) != sizeRO || complexImResized.get_size(1) != sizeE1 || complexImResized.get_size(2) != sizeE2)
                {
                    dim[0] = sizeRO;
                    dim[1] = sizeE1;
                    dim[2] = sizeE2;
                    complexImResized.create(dim);
                }

                Gadgetron::clear(&complexImResized);

                hoNDArray<T> kspace(complexIm);
                Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->fft3c(complexIm, kspace);

                Gadgetron::pad(sizeRO, sizeE1, sizeE2, &kspace, &complexImResized);

                Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft3c(complexImResized);

                typename realType<T>::Type scaling = (typename realType<T>::Type)(std::sqrt((double)sizeRO*sizeE1*sizeE2) / std::sqrt((double)RO*E1));
                Gadgetron::scal(scaling, complexImResized);
            }
            else
            {
                if (RO == sizeRO && E1 == sizeE1)
                {
                    complexImResized = complexIm;
                    return;
                }

                if (complexImResized.get_size(0) != sizeRO || complexImResized.get_size(1) != sizeE1)
                {
                    dim[0] = sizeRO;
                    dim[1] = sizeE1;
                    complexImResized.create(dim);
                }

                Gadgetron::clear(&complexImResized);

                hoNDArray<T> kspace(complexIm);
                Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->fft2c(complexIm, kspace);

                Gadgetron::pad(sizeRO, sizeE1, &kspace, &complexImResized);
                Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft2c(complexImResized);

                typename realType<T>::Type scaling = (typename realType<T>::Type)(std::sqrt((double)sizeRO*sizeE1) / std::sqrt((double)RO*E1));
                Gadgetron::scal(scaling, complexImResized);
            }
        }
        catch (...)
        {
            GADGET_THROW("Errors in zero_pad_resize(...) ... ");
        }
    }

    template EXPORTMRICORE void zero_pad_resize(const hoNDArray< std::complex<float> >& complexIm, size_t sizeRO, size_t sizeE1, size_t sizeE2, hoNDArray< std::complex<float> >& complexImResized);
    template EXPORTMRICORE void zero_pad_resize(const hoNDArray< std::complex<double> >& complexIm, size_t sizeRO, size_t sizeE1, size_t sizeE2, hoNDArray< std::complex<double> >& complexImResized);

    // ------------------------------------------------------------------------

    void get_debug_folder_path(const std::string& debugFolder, std::string& debugFolderPath)
    {
        char* v = std::getenv("GADGETRON_DEBUG_FOLDER");
        if (v == NULL)
        {
#ifdef _WIN32
            debugFolderPath = "c:/temp/gadgetron";
#else
            debugFolderPath = "/tmp/gadgetron";
#endif // _WIN32
        }
        else
        {
            debugFolderPath = std::string(v);
        }

        debugFolderPath.append("/");
        debugFolderPath.append(debugFolder);
        debugFolderPath.append("/");
    }
}
