
/** \file   mri_core_utility.cpp
    \brief  Implementation useful utility functionalities for 2D and 3D MRI parallel imaging
    \author Hui Xue
*/

#include "mri_core_utility.h"
#include "hoNDArray_elemwise.h"
#include "hoNDArray_utils.h"
#include "hoNDFFT.h"
#include "mri_core_kspace_filter.h"
#include <ctime>

namespace Gadgetron
{
    template <typename T>
    hoNDArray<bool> detect_readout_sampling_status(const hoNDArray<T>& data)
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

            return sampled;
        }
        catch (...)
        {
            GADGET_THROW("Errors in detect_readout_sampling_status(...) ... ");
        }
    }

    template EXPORTMRICORE hoNDArray<bool> detect_readout_sampling_status(const hoNDArray< float >& data);
    template EXPORTMRICORE hoNDArray<bool> detect_readout_sampling_status(const hoNDArray< double >& data);
    template EXPORTMRICORE hoNDArray<bool> detect_readout_sampling_status(const hoNDArray< std::complex<float> >& data);
    template EXPORTMRICORE hoNDArray<bool> detect_readout_sampling_status(const hoNDArray< std::complex<double> >& data);

    // ------------------------------------------------------------------------

    template <typename T>
    std::tuple<size_t, size_t> detect_sampled_region_E1(const hoNDArray<T>& data)
    {
        try
        {
            hoNDArray<bool> sampled = Gadgetron::detect_readout_sampling_status(data);

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
            hoNDArray<bool> sampled = Gadgetron::detect_readout_sampling_status(data);

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

    template <typename T> 
    void compute_averaged_data_N_S(const hoNDArray<T>& data, bool average_N, bool average_S, bool count_sampling_freq, hoNDArray<T>& res)
    {
        try
        {
            size_t RO = data.get_size(0);
            size_t E1 = data.get_size(1);
            size_t E2 = data.get_size(2);
            size_t CHA = data.get_size(3);
            size_t N = data.get_size(4);
            size_t S = data.get_size(5);
            size_t SLC = data.get_size(6);

            typedef typename realType<T>::Type value_type;

            if (average_N)
            {
                if (N > 1)
                {
                    if (count_sampling_freq)
                    {
                        hoNDArray<bool> sampled = Gadgetron::detect_readout_sampling_status(data);
                        Gadgetron::sum_over_dimension(data, res, 4);

                        // for every E1/E2 location, count how many times it is sampled for all N
                        size_t ro, e1, e2, cha, n, s, slc;
                        for (slc = 0; slc < SLC; slc++)
                        {
                            for (cha = 0; cha < CHA; cha++)
                            {
                                for (e2 = 0; e2 < E2; e2++)
                                {
                                    for (e1 = 0; e1 < E1; e1++)
                                    {
                                        float freq = 0;

                                        for (s = 0; s < S; s++)
                                        {
                                            for (n = 0; n < N; n++)
                                            {
                                                if (sampled(e1, e2, n, s, slc)) freq += 1;
                                            }

                                            if (freq > 1)
                                            {
                                                value_type freq_reciprocal = (value_type)(1.0 / freq);
                                                T* pAve = &(res(0, e1, e2, cha, 0, s, slc));
                                                for (ro = 0; ro < RO; ro++) pAve[ro] *= freq_reciprocal;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                    else
                    {
                        Gadgetron::sum_over_dimension(data, res, (size_t)4);
                        Gadgetron::scal((value_type)(1.0 / N), res);
                    }
                }
                else
                {
                    res = data;
                }
            }
            else
            {
                res = data;
            }

            if (average_S)
            {
                if (S > 1)
                {
                    hoNDArray<T> dataTmp;
                    Gadgetron::sum_over_dimension(res, dataTmp, 5);
                    Gadgetron::scal((value_type)(1.0 / S), dataTmp);
                    res = dataTmp;
                }
            }
        }
        catch (...)
        {
            GADGET_THROW("Errors in compute_averaged_data_N_S(...) ... ");
        }
    }

    template EXPORTMRICORE void compute_averaged_data_N_S(const hoNDArray< std::complex<float> >& data, bool average_N, bool average_S, bool count_sampling_freq, hoNDArray< std::complex<float> >& res);
    template EXPORTMRICORE void compute_averaged_data_N_S(const hoNDArray< std::complex<double> >& data, bool average_N, bool average_S, bool count_sampling_freq, hoNDArray< std::complex<double> >& res);

    // ------------------------------------------------------------------------

    template <typename T> 
    void compute_eigen_channel_coefficients(const hoNDArray<T>& data, bool average_N, bool average_S, bool count_sampling_freq, size_t N, size_t S, double coil_compression_thres, size_t compression_num_modesKept, std::vector< std::vector< std::vector< hoNDKLT<T> > > >& KLT)
    {
        try
        {
            size_t RO = data.get_size(0);
            size_t E1 = data.get_size(1);
            size_t E2 = data.get_size(2);
            size_t CHA = data.get_size(3);
            size_t dataN = data.get_size(4);
            size_t dataS = data.get_size(5);
            size_t SLC = data.get_size(6);

            typedef typename realType<T>::Type value_type;

            size_t n, s, slc;

            hoNDArray<T> dataAve;
            GADGET_CATCH_THROW( Gadgetron::compute_averaged_data_N_S(data, average_N, average_S, count_sampling_freq, dataAve) );

            size_t dataAveN = dataAve.get_size(4);
            size_t dataAveS = dataAve.get_size(5);

            if(KLT.size()!=SLC) KLT.resize(SLC);
            for (slc = 0; slc < SLC; slc++)
            {
                if (KLT[slc].size() != S) KLT[slc].resize(S);
                for (s = 0; s < S; s++)
                {
                    if (KLT[slc][s].size() != N) KLT[slc][s].resize(N);
                }
            }

            for (slc = 0; slc < SLC; slc++)
            {
                for (s = 0; s < S; s++)
                {
                    size_t s_used = s;
                    if (s_used >= dataAveS) s_used = dataAveS - 1;

                    for (n = 0; n < N; n++)
                    {
                        size_t n_used = n;
                        if (n_used >= dataAveN) n_used = dataAveN - 1;

                        T* pDataAve = &(dataAve(0, 0, 0, 0, n_used, s_used, slc));
                        hoNDArray<T> dataUsed(RO, E1, E2, CHA, pDataAve);

                        if (slc == 0 && n == 0 && s == 0)
                        {
                            if (compression_num_modesKept > 0)
                            {
                                KLT[slc][s][n].prepare(dataUsed, 3, compression_num_modesKept);
                            }
                            else if (coil_compression_thres > 0)
                            {
                                KLT[slc][s][n].prepare(dataUsed, 3, (value_type)(coil_compression_thres));
                            }
                            else
                            {
                                KLT[slc][s][n].prepare(dataUsed, 3, (size_t)(0));
                            }
                        }
                        else
                        {
                            KLT[slc][s][n].prepare(dataUsed, 3, KLT[0][0][0].output_length());
                        }
                    }
                }
            }
        }
        catch (...)
        {
            GADGET_THROW("Errors in compute_eigen_channel_coefficients(...) ... ");
        }
    }

    template EXPORTMRICORE void compute_eigen_channel_coefficients(const hoNDArray< std::complex<float> >& data, bool average_N, bool average_S, bool count_sampling_freq, size_t N, size_t S, double coil_compression_thres, size_t compression_num_modesKept, std::vector< std::vector< std::vector< hoNDKLT< std::complex<float> > > > >& KLT);
    template EXPORTMRICORE void compute_eigen_channel_coefficients(const hoNDArray< std::complex<double> >& data, bool average_N, bool average_S, bool count_sampling_freq, size_t N, size_t S, double coil_compression_thres, size_t compression_num_modesKept, std::vector< std::vector< std::vector< hoNDKLT< std::complex<double> > > > >& KLT);

    // ------------------------------------------------------------------------

    template <typename T> 
    void apply_eigen_channel_coefficients(const std::vector< std::vector< std::vector< hoNDKLT<T> > > >& KLT, hoNDArray<T>& data)
    {
        try
        {
            size_t RO = data.get_size(0);
            size_t E1 = data.get_size(1);
            size_t E2 = data.get_size(2);
            size_t CHA = data.get_size(3);
            size_t N = data.get_size(4);
            size_t S = data.get_size(5);
            size_t SLC = data.get_size(6);

            GADGET_CHECK_THROW(KLT.size() == SLC);

            size_t dstCHA = KLT[0][0][0].output_length();

            hoNDArray<T> dstData;
            dstData.create(RO, E1, E2, dstCHA, N, S, SLC);

            size_t n, s, slc;
            for (slc = 0; slc < SLC; slc++)
            {
                for (s = 0; s < S; s++)
                {
                    size_t s_KLT = s;
                    if (s_KLT >= KLT[slc].size()) s_KLT = KLT[slc].size()-1;

                    for (n = 0; n < N; n++)
                    {
                        size_t n_KLT = n;
                        if (n_KLT >= KLT[slc][s_KLT].size()) n_KLT = KLT[slc][s_KLT].size()-1;

                        T* pData = &(data(0, 0, 0, 0, n, s, slc));
                        hoNDArray<T> data_in(RO, E1, E2, CHA, pData);

                        T* pDstData = &(dstData(0, 0, 0, 0, n, s, slc));
                        hoNDArray<T> data_out(RO, E1, E2, dstCHA, pDstData);

                        KLT[slc][s_KLT][n_KLT].transform(data_in, data_out, 3);
                    }
                }
            }

            data = dstData;
        }
        catch (...)
        {
            GADGET_THROW("Errors in apply_eigen_channel_coefficients(...) ... ");
        }
    }

    template EXPORTMRICORE void apply_eigen_channel_coefficients(const std::vector< std::vector< std::vector< hoNDKLT< std::complex<float> > > > >& KLT, hoNDArray< std::complex<float> >& data);
    template EXPORTMRICORE void apply_eigen_channel_coefficients(const std::vector< std::vector< std::vector< hoNDKLT< std::complex<double> > > > >& KLT, hoNDArray< std::complex<double> >& data);

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
