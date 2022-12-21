
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

            hoNDArray<bool> sampled(E1, E2, N, S, SLC);



            for (size_t slc = 0; slc < SLC; slc++)
            {
                for (size_t s = 0; s < S; s++)
                {
                    for (size_t n = 0; n < N; n++)
                    {
                        for (size_t e2 = 0; e2 < E2; e2++)
                        {
                            for (size_t e1 = 0; e1 < E1; e1++)
                            {
                                value_type v = 0;
                                const T* pData = &(data(0, e1, e2, 0, n, s, slc));


                                for (size_t ro = 0; ro < RO; ro++)
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



        for (size_t slc = 0; slc < SLC; slc++)
        {
            for (size_t s = 0; s < S; s++)
            {
                for (size_t n = 0; n < N; n++)
                {
                    for (size_t e2 = 0; e2 < E2; e2++)
                    {
                        for (size_t e1 = 0; e1 < E1; e1++)
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

                Gadgetron::pad(sizeRO, sizeE1, sizeE2, kspace, complexImResized);

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

                Gadgetron::pad(sizeRO, sizeE1, kspace, complexImResized);
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

    template EXPORTMRICORE void compute_averaged_data_N_S(const hoNDArray< std::complex<float> >& data, bool average_N, bool average_S, bool count_sampling_freq, hoNDArray< std::complex<float> >& res);
    template EXPORTMRICORE void compute_averaged_data_N_S(const hoNDArray< std::complex<double> >& data, bool average_N, bool average_S, bool count_sampling_freq, hoNDArray< std::complex<double> >& res);

    // ------------------------------------------------------------------------

    template <typename T> 
    void select_data_N_S(const hoNDArray<T>& data, bool select_N, size_t n, bool select_S, size_t s, hoNDArray<T>& res)
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

            std::vector<size_t> dim;
            data.get_dimensions(dim);

            std::vector<size_t> dimR(dim);

            size_t slc, ss, nn;

            if(!select_N && !select_S)
            {
                res = data;
            }
            else if(select_N && !select_S)
            {
                if(n<N)
                {
                    dimR[4] = 1;
                    res.create(dimR);

                    for (slc = 0; slc < SLC; slc++)
                    {
                        for (ss = 0; ss < S; ss++)
                        {
                            const T* pData = &(data(0, 0, 0, 0, n, ss, slc));
                            T* pRes = &res(0, 0, 0, 0, 0, ss, slc);
                            memcpy(pRes, pData, sizeof(T)*RO*E1*E2*CHA);
                        }
                    }
                }
                else
                {
                    GADGET_THROW("select_data_N_S, n>=N");
                }
            }
            else if (!select_N && select_S)
            {
                if (s<S)
                {
                    dimR[5] = 1;
                    res.create(dimR);

                    for (slc = 0; slc < SLC; slc++)
                    {
                        for (nn = 0; nn < N; nn++)
                        {
                            const T* pData = &(data(0, 0, 0, 0, nn, s, slc));
                            T* pRes = &res(0, 0, 0, 0, nn, 0, slc);
                            memcpy(pRes, pData, sizeof(T)*RO*E1*E2*CHA);
                        }
                    }
                }
                else
                {
                    GADGET_THROW("select_data_N_S, s>=S");
                }
            }
            else
            {
                GADGET_CHECK_THROW(n < N);
                GADGET_CHECK_THROW(s < S);

                dimR[4] = 1;
                dimR[5] = 1;
                res.create(dimR);

                for (slc = 0; slc < SLC; slc++)
                {
                    const T* pData = &(data(0, 0, 0, 0, n, s, slc));
                    T* pRes = &res(0, 0, 0, 0, 0, 0, slc);
                    memcpy(pRes, pData, sizeof(T)*RO*E1*E2*CHA);
                }
            }
        }
        catch (...)
        {
            GADGET_THROW("Errors in select_data_N_S(...) ... ");
        }
    }

    template EXPORTMRICORE void select_data_N_S(const hoNDArray< std::complex<float> >& data, bool select_N, size_t n, bool select_S, size_t s, hoNDArray< std::complex<float> >& res);
    template EXPORTMRICORE void select_data_N_S(const hoNDArray< std::complex<double> >& data, bool select_N, size_t n, bool select_S, size_t s, hoNDArray< std::complex<double> >& res);

    // ------------------------------------------------------------------------

    template <typename T> 
    void compute_eigen_channel_coefficients(const hoNDArray<T>& data, bool average_N, bool average_S, bool count_sampling_freq, size_t N, size_t S, double coil_compression_thres, size_t compression_num_modesKept, std::vector< std::vector< std::vector< hoNDKLT<T> > > >& KLT)
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
                        if(n>=dataAveN && s>=dataAveS)
                        {
                            KLT[slc][s][n] = KLT[slc][dataAveS - 1][dataAveN-1];
                        }
                        else
                        {
                            KLT[slc][s][n].prepare(dataUsed, 3, KLT[0][0][0].output_length());
                        }
                    }
                }
            }
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

    // ------------------------------------------------------------------------

    void find_calib_mode(ISMRMRD::IsmrmrdHeader& h, Gadgetron::ismrmrdCALIBMODE& CalibMode, Gadgetron::IsmrmrdDIM& InterleaveDim, double& acceFactorE1, double& acceFactorE2, bool verbose)
    {
        try
        {
            if (!h.encoding[0].parallelImaging)
            {
                GADGET_THROW("Parallel Imaging section not found in header");
            }

            ISMRMRD::ParallelImaging p_imaging = *h.encoding[0].parallelImaging;

            acceFactorE1 = (double)(p_imaging.accelerationFactor.kspace_encoding_step_1);
            acceFactorE2 = (double)(p_imaging.accelerationFactor.kspace_encoding_step_2);

            GDEBUG_CONDITION_STREAM(verbose, "acceFactorE1 is " << acceFactorE1);
            GDEBUG_CONDITION_STREAM(verbose, "acceFactorE2 is " << acceFactorE2);

            if (!p_imaging.calibrationMode.is_present())
            {
                GADGET_THROW("Parallel calibration mode not found in header");
            }

            std::string calib = *p_imaging.calibrationMode;
            if (calib.compare("interleaved") == 0)
            {
                CalibMode = Gadgetron::ISMRMRD_interleaved;
                GDEBUG_CONDITION_STREAM(verbose, "Calibration mode is interleaved");

                if (p_imaging.interleavingDimension)
                {
                    if (p_imaging.interleavingDimension->compare("phase") == 0)
                    {
                        InterleaveDim = Gadgetron::DIM_Phase;
                    }
                    else if (p_imaging.interleavingDimension->compare("repetition") == 0)
                    {
                        InterleaveDim = Gadgetron::DIM_Repetition;
                    }
                    else if (p_imaging.interleavingDimension->compare("average") == 0)
                    {
                        InterleaveDim = Gadgetron::DIM_Average;
                    }
                    else if (p_imaging.interleavingDimension->compare("contrast") == 0)
                    {
                        InterleaveDim = Gadgetron::DIM_Contrast;
                    }
                    else if (p_imaging.interleavingDimension->compare("other") == 0)
                    {
                        InterleaveDim = Gadgetron::DIM_other1;
                    }
                    else
                    {
                        GADGET_THROW("Unknown interleaving dimension. Bailing out");
                    }
                }
            }
            else if (calib.compare("embedded") == 0)
            {
                CalibMode = Gadgetron::ISMRMRD_embedded;
                GDEBUG_CONDITION_STREAM(verbose, "Calibration mode is embedded");
            }
            else if (calib.compare("separate") == 0)
            {
                CalibMode = Gadgetron::ISMRMRD_separate;
                GDEBUG_CONDITION_STREAM(verbose, "Calibration mode is separate");
            }
            else if (calib.compare("external") == 0)
            {
                CalibMode = Gadgetron::ISMRMRD_external;
            }
            else if ((calib.compare("other") == 0) && acceFactorE1 == 1 && acceFactorE2 == 1)
            {
                CalibMode = Gadgetron::ISMRMRD_noacceleration;
                acceFactorE1 = 1;
            }
            else if ((calib.compare("other") == 0) && (acceFactorE1>1 || acceFactorE2>1))
            {
                CalibMode = Gadgetron::ISMRMRD_interleaved;
                acceFactorE1 = 2;
                InterleaveDim = Gadgetron::DIM_Phase;
            }
            else
            {
                GADGET_THROW("Failed to process parallel imaging calibration mode");
            }
        }
        catch (...)
        {
            GADGET_THROW("Error happened in findCalibMode(...) ... ");
        }
    }

    // ------------------------------------------------------------------------

    void find_encoding_limits(ISMRMRD::IsmrmrdHeader& h, ISMRMRD::EncodingCounters& meas_max_idx, bool verbose)
    {
        try
        {
            ISMRMRD::EncodingSpace e_space = h.encoding[0].encodedSpace;
            ISMRMRD::EncodingSpace r_space = h.encoding[0].reconSpace;
            ISMRMRD::EncodingLimits e_limits = h.encoding[0].encodingLimits;

            meas_max_idx.kspace_encode_step_1 = (uint16_t)e_space.matrixSize.y - 1;

            meas_max_idx.set = (e_limits.set && (e_limits.set->maximum>0)) ? e_limits.set->maximum : 0;
            meas_max_idx.phase = (e_limits.phase && (e_limits.phase->maximum>0)) ? e_limits.phase->maximum : 0;

            meas_max_idx.kspace_encode_step_2 = (uint16_t)e_space.matrixSize.z - 1;

            meas_max_idx.contrast = (e_limits.contrast && (e_limits.contrast->maximum > 0)) ? e_limits.contrast->maximum : 0;

            meas_max_idx.slice = (e_limits.slice && (e_limits.slice->maximum > 0)) ? e_limits.slice->maximum : 0;

            meas_max_idx.repetition = e_limits.repetition ? e_limits.repetition->maximum : 0;

            meas_max_idx.average = e_limits.average ? e_limits.average->maximum : 0;

            // always combine the SEG
            meas_max_idx.segment = 0;
        }
        catch (...)
        {
            GADGET_THROW("Error happened in findEncodingLimits(...) ... ");
        }
    }

    // ------------------------------------------------------------------------

    void find_matrix_size_encoding(ISMRMRD::IsmrmrdHeader& h, size_t matrix_size_encoding[3])
    {
        matrix_size_encoding[0] = h.encoding[0].encodedSpace.matrixSize.x;
        matrix_size_encoding[1] = h.encoding[0].encodedSpace.matrixSize.y;
        matrix_size_encoding[2] = h.encoding[0].encodedSpace.matrixSize.z;
    }

    // ------------------------------------------------------------------------

    void find_FOV_encoding(ISMRMRD::IsmrmrdHeader& h, float field_of_view_encoding[3])
    {
        field_of_view_encoding[0] = h.encoding[0].encodedSpace.fieldOfView_mm.x;
        field_of_view_encoding[1] = h.encoding[0].encodedSpace.fieldOfView_mm.y;
        field_of_view_encoding[2] = h.encoding[0].encodedSpace.fieldOfView_mm.z;
    }

    // ------------------------------------------------------------------------

    void find_matrix_size_recon(ISMRMRD::IsmrmrdHeader& h, size_t matrix_size_recon[3])
    {
        matrix_size_recon[0] = h.encoding[0].reconSpace.matrixSize.x;
        matrix_size_recon[1] = h.encoding[0].reconSpace.matrixSize.y;
        matrix_size_recon[2] = h.encoding[0].reconSpace.matrixSize.z;
    }

    // ------------------------------------------------------------------------

    void find_FOV_recon(ISMRMRD::IsmrmrdHeader& h, float field_of_view_recon[3])
    {
        field_of_view_recon[0] = h.encoding[0].reconSpace.fieldOfView_mm.x;
        field_of_view_recon[1] = h.encoding[0].reconSpace.fieldOfView_mm.y;
        field_of_view_recon[2] = h.encoding[0].reconSpace.fieldOfView_mm.z;
    }

    // ------------------------------------------------------------------------

    void get_current_moment(std::string& procTime)
    {
        char timestamp[100];
        time_t mytime;
        struct tm *mytm;
        mytime = time(NULL);
        mytm = localtime(&mytime);
        strftime(timestamp, sizeof(timestamp), "%a, %b %d %Y, %H:%M:%S", mytm);
        procTime = timestamp;
    }

    // ------------------------------------------------------------------------

    void get_ismrmrd_meta_values(const ISMRMRD::MetaContainer& attrib, const std::string& name, std::vector<long>& v)
    {
        try
        {
            size_t num = attrib.length(name.c_str());
            if (num == 0)
            {
                v.clear();
                GWARN_STREAM("get_ismrmrd_meta_values, can not find field : " << name);
                return;
            }

            v.resize(num);

            size_t ii;
            for (ii = 0; ii<num; ii++)
            {
                v[ii] = attrib.as_long(name.c_str(), ii);
            }
        }
        catch (...)
        {
            GADGET_THROW("Error happened in get_ismrmrd_meta_values(const ISMRMRD::MetaContainer& attrib, const std::string& name, std::vector<long>& v) ... ");
        }
    }

    void get_ismrmrd_meta_values(const ISMRMRD::MetaContainer& attrib, const std::string& name, std::vector<double>& v)
    {
        try
        {
            size_t num = attrib.length(name.c_str());
            if (num == 0)
            {
                v.clear();
                GWARN_STREAM("get_ismrmrd_meta_values, can not find field : " << name);
                return;
            }

            v.resize(num);

            size_t ii;
            for (ii = 0; ii<num; ii++)
            {
                v[ii] = attrib.as_double(name.c_str(), ii);
            }
        }
        catch (...)
        {
            GADGET_THROW("Error happened in get_ismrmrd_meta_values(const ISMRMRD::MetaContainer& attrib, const std::string& name, std::vector<double>& v) ... ");
        }
    }

    void get_ismrmrd_meta_values(const ISMRMRD::MetaContainer& attrib, const std::string& name, std::vector<std::string>& v)
    {
        try
        {
            size_t num = attrib.length(name.c_str());
            if (num == 0)
            {
                v.clear();
                GWARN_STREAM("get_ismrmrd_meta_values, can not find field : " << name);
                return;
            }

            v.resize(num);

            size_t ii;
            for (ii = 0; ii<num; ii++)
            {
                v[ii] = std::string(attrib.as_str(name.c_str(), ii));
            }
        }
        catch (...)
        {
            GADGET_THROW("Error happened in get_ismrmrd_meta_values(const ISMRMRD::MetaContainer& attrib, const std::string& name, std::vector<std::string>& v) ... ");
        }
    }

    // ------------------------------------------------------------------------

    template <typename T>
    void set_ismrmrd_meta_values(ISMRMRD::MetaContainer& attrib, const std::string& name, const std::vector<T>& v)
    {
        try
        {
            size_t num = v.size();
            if (num == 0)
            {
                GWARN_STREAM("setISMRMRMetaValues, input vector is empty ... " << name);
                return;
            }

            attrib.set(name.c_str(), v[0]);

            size_t ii;
            for (ii = 1; ii<v.size(); ii++)
            {
                attrib.append(name.c_str(), v[ii]);
            }
        }
        catch (...)
        {
            GADGET_THROW("Error happened in setISMRMRMetaValues(ISMRMRD::MetaContainer& attrib, const std::string& name, const std::vector<T>& v) ... ");
        }
    }

    template EXPORTMRICORE void set_ismrmrd_meta_values(ISMRMRD::MetaContainer& attrib, const std::string& name, const std::vector<long>& v);
    template EXPORTMRICORE void set_ismrmrd_meta_values(ISMRMRD::MetaContainer& attrib, const std::string& name, const std::vector<double>& v);

    void set_ismrmrd_meta_values(ISMRMRD::MetaContainer& attrib, const std::string& name, const std::vector<std::string>& v)
    {
        try
        {
            size_t num = v.size();
            if (num == 0)
            {
                GWARN_STREAM("setISMRMRMetaValues, input vector is empty ... " << name);
                return;
            }

            attrib.set(name.c_str(), v[0].c_str());

            size_t ii;
            for (ii = 1; ii<v.size(); ii++)
            {
                attrib.append(name.c_str(), v[ii].c_str());
            }
        }
        catch (...)
        {
            GADGET_THROW("Error happened in setISMRMRMetaValues(ISMRMRD::MetaContainer& attrib, const std::string& name, const std::vector<std::string>& v) ... ");
        }
    }

    // ------------------------------------------------------------------------

    template <typename T>
    void append_ismrmrd_meta_values(ISMRMRD::MetaContainer& attrib, const std::string& name, const std::vector<T>& v)
    {
        try
        {
            size_t num = v.size();
            if (num == 0)
            {
                GWARN_STREAM("append_ismrmrd_meta_values, input vector is empty ... " << name);
                return;
            }

            attrib.append(name.c_str(), v[0]);

            size_t ii;
            for (ii = 1; ii<v.size(); ii++)
            {
                attrib.append(name.c_str(), v[ii]);
            }
        }
        catch (...)
        {
            GADGET_THROW("Error happened in append_ismrmrd_meta_values(ISMRMRD::MetaContainer& attrib, const std::string& name, const std::vector<T>& v) ... ");
        }
    }

    template EXPORTMRICORE void append_ismrmrd_meta_values(ISMRMRD::MetaContainer& attrib, const std::string& name, const std::vector<long>& v);
    template EXPORTMRICORE void append_ismrmrd_meta_values(ISMRMRD::MetaContainer& attrib, const std::string& name, const std::vector<double>& v);

    void append_ismrmrd_meta_values(ISMRMRD::MetaContainer& attrib, const std::string& name, const std::vector<std::string>& v)
    {
        try
        {
            size_t num = v.size();
            if (num == 0)
            {
                GWARN_STREAM("append_ismrmrd_meta_values, input vector is empty ... " << name);
                return;
            }

            attrib.append(name.c_str(), v[0].c_str());

            size_t ii;
            for (ii = 1; ii<v.size(); ii++)
            {
                attrib.append(name.c_str(), v[ii].c_str());
            }
        }
        catch (...)
        {
            GADGET_THROW("Error happened in append_ismrmrd_meta_values(ISMRMRD::MetaContainer& attrib, const std::string& name, const std::vector<std::string>& v) ... ");
        }
    }

    // ------------------------------------------------------------------------

    void PatientCoordinateSystem_to_DeviceCoordinateSystem(double& x, double& y, double& z, const std::string& position)
    {
        // this is following dicom tag (0020, 0037)

        if (position == "HFS") // Head-first supine (HFS)
        {
            y = -y;
            z = -z;
        }
        else if (position == "HFP") // Head-first prone (HFP)
        {
            x = -x;
            z = -z;
        }
        else if (position == "HFDR") // Head-first decubitus-right
        {
            double v = x;
            x = y;
            y = v;
            z = -z;
        }
        else if (position == "HFDL") // Head-first decubitus-left (HFDL)
        {
            double v = x;
            x = y;
            y = v;

            x = -x;
            y = -y;
            z = -z;
        }
        else if (position == "FFDR") // Feet-first decubitus-right (FFDR)
        {
            double v = x;
            x = y;
            y = v;

            x = -x;
        }
        else if (position == "FFDL") // Feet-first decubitus-left (FFDL)
        {
            double v = x;
            x = y;
            y = v;

            y = -y;
        }
        else if (position == "FFP") // Feet-first prone (FFP)
        {
        }
        else if (position == "FFS") // Feet-first supine (FFS)
        {
            x = -x;
            y = -y;
        }
        else
        {
            GADGET_THROW("Unknown position string ... ");
        }
    }

    void DeviceCoordinateSystem_to_PatientCoordinateSystem(double& x, double& y, double& z, const std::string& position)
    {
        if (position == "HFS") // Head-first supine (HFS)
        {
            y = -y;
            z = -z;
        }
        else if (position == "HFP") // Head-first prone (HFP)
        {
            x = -x;
            z = -z;
        }
        else if (position == "HFDR") // Head-first decubitus-right
        {
            double v = x;
            x = y;
            y = v;
            z = -z;
        }
        else if (position == "HFDL") // Head-first decubitus-left (HFDL)
        {
            double v = x;
            x = y;
            y = v;

            x = -x;
            y = -y;
            z = -z;
        }
        else if (position == "FFDR") // Feet-first decubitus-right (FFDR)
        {
            double v = x;
            x = y;
            y = v;

            y = -y;
        }
        else if (position == "FFDL") // Feet-first decubitus-left (FFDL)
        {
            double v = x;
            x = y;
            y = v;

            x = -x;
        }
        else if (position == "FFP") // Feet-first prone (FFP)
        {
        }
        else if (position == "FFS") // Feet-first supine (FFS)
        {
            x = -x;
            y = -y;
        }
        else
        {
            GADGET_THROW("Unknown position string ... ");
        }
    }

    // ------------------------------------------------------------------------

    bool check_idential_slice_prescription(ISMRMRD::ISMRMRD_ImageHeader a, ISMRMRD::ISMRMRD_ImageHeader b)
    {
        GADGET_CHECK_RETURN_FALSE(a.matrix_size[0] == b.matrix_size[0]);
        GADGET_CHECK_RETURN_FALSE(a.matrix_size[1] == b.matrix_size[1]);
        GADGET_CHECK_RETURN_FALSE(a.matrix_size[2] == b.matrix_size[2]);

        GADGET_CHECK_RETURN_FALSE(std::abs(a.field_of_view[0] - b.field_of_view[0])<0.1);
        GADGET_CHECK_RETURN_FALSE(std::abs(a.field_of_view[1] - b.field_of_view[1])<0.1);
        GADGET_CHECK_RETURN_FALSE(std::abs(a.field_of_view[2] - b.field_of_view[2])<0.1);

        GADGET_CHECK_RETURN_FALSE(std::abs(a.patient_table_position[0] - b.patient_table_position[0])<0.001);
        GADGET_CHECK_RETURN_FALSE(std::abs(a.patient_table_position[1] - b.patient_table_position[1])<0.001);
        GADGET_CHECK_RETURN_FALSE(std::abs(a.patient_table_position[2] - b.patient_table_position[2])<0.001);

        GADGET_CHECK_RETURN_FALSE(std::abs(a.position[0] - b.position[0])<0.001);
        GADGET_CHECK_RETURN_FALSE(std::abs(a.position[1] - b.position[1])<0.001);
        GADGET_CHECK_RETURN_FALSE(std::abs(a.position[2] - b.position[2])<0.001);

        GADGET_CHECK_RETURN_FALSE(std::abs(a.read_dir[0] - b.read_dir[0])<0.001);
        GADGET_CHECK_RETURN_FALSE(std::abs(a.read_dir[1] - b.read_dir[1])<0.001);
        GADGET_CHECK_RETURN_FALSE(std::abs(a.read_dir[2] - b.read_dir[2])<0.001);

        GADGET_CHECK_RETURN_FALSE(std::abs(a.phase_dir[0] - b.phase_dir[0])<0.001);
        GADGET_CHECK_RETURN_FALSE(std::abs(a.phase_dir[1] - b.phase_dir[1])<0.001);
        GADGET_CHECK_RETURN_FALSE(std::abs(a.phase_dir[2] - b.phase_dir[2])<0.001);

        return true;
    }

    // ------------------------------------------------------------------------

    std::string get_ismrmrd_dim_name(const IsmrmrdDIM& dim)
    {
        std::ostringstream os;
        switch (dim)
        {
        case DIM_ReadOut:
            os << "DIM_ReadOut";
            break;

        case DIM_Encoding1:
            os << "DIM_Encoding1";
            break;

        case DIM_Channel:
            os << "DIM_Channel";
            break;

        case DIM_Slice:
            os << "DIM_Slice";
            break;

        case DIM_Encoding2:
            os << "DIM_Encoding2";
            break;

        case DIM_Contrast:
            os << "DIM_Contrast";
            break;

        case DIM_Phase:
            os << "DIM_Phase";
            break;

        case DIM_Repetition:
            os << "DIM_Repetition";
            break;

        case DIM_Set:
            os << "DIM_Set";
            break;

        case DIM_Segment:
            os << "DIM_Segment";
            break;

        case DIM_Average:
            os << "DIM_Average";
            break;

        case DIM_other1:
            os << "DIM_other1";
            break;

        case DIM_other2:
            os << "DIM_other2";
            break;

        case DIM_other3:
            os << "DIM_other3";
            break;

        default:
            os << "DIM_NONE";
        }

        std::string dimStr(os.str());
        return dimStr;
    }

    // ------------------------------------------------------------------------

    IsmrmrdDIM get_ismrmrd_dim_from_name(const std::string& name)
    {
        if (name == "DIM_ReadOut") return DIM_ReadOut;
        if (name == "DIM_Encoding1") return DIM_Encoding1;
        if (name == "DIM_Channel") return DIM_Channel;
        if (name == "DIM_Slice") return DIM_Slice;
        if (name == "DIM_Encoding2") return DIM_Encoding2;
        if (name == "DIM_Contrast") return DIM_Contrast;
        if (name == "DIM_Phase") return DIM_Phase;
        if (name == "DIM_Repetition") return DIM_Repetition;
        if (name == "DIM_Set") return DIM_Set;
        if (name == "DIM_Segment") return DIM_Segment;
        if (name == "DIM_Average") return DIM_Average;
        if (name == "DIM_other1") return DIM_other1;
        if (name == "DIM_other2") return DIM_other2;
        if (name == "DIM_other3") return DIM_other3;

        return DIM_NONE;
    }
    namespace {
        template<class T> std::map<std::string, decltype(T::value)> to_map_internal(const std::vector<T>& userparameters){
            std::map<std::string, decltype(T::value)> output_map;
            for (auto element : userparameters){
                output_map.insert({element.name, element.value});
            }
            return output_map;
        }
    }

    std::map<std::string, long> to_map(const std::vector<ISMRMRD::UserParameterLong> & userparameters) {
        return to_map_internal(userparameters);
    }

    std::map<std::string, double> to_map(const std::vector<ISMRMRD::UserParameterDouble> & userparameters) {
        return to_map_internal(userparameters);
    }

    std::map<std::string, std::string> to_map(const std::vector<ISMRMRD::UserParameterString> & userparameters) {
        return to_map_internal(userparameters);
    }

    ISMRMRD::ImageHeader image_header_from_acquisition(
        const ISMRMRD::AcquisitionHeader& acq_header, const ISMRMRD::IsmrmrdHeader& header, const hoNDArray<std::complex<float>>& data) {

        ISMRMRD::ImageHeader im_header;

        im_header.version         = acq_header.version;
        im_header.data_type       = ISMRMRD::ISMRMRD_CXFLOAT;
        im_header.measurement_uid = acq_header.measurement_uid;

        im_header.matrix_size[0] = (uint16_t)data.get_size(0);
        im_header.matrix_size[1] = (uint16_t)data.get_size(1);
        im_header.matrix_size[2] = (uint16_t)data.get_size(2);

        im_header.field_of_view[0] = header.encoding[0].reconSpace.fieldOfView_mm.x;
        im_header.field_of_view[1] = header.encoding[0].reconSpace.fieldOfView_mm.y;
        im_header.field_of_view[2] = header.encoding[0].reconSpace.fieldOfView_mm.z;

        im_header.channels = (uint16_t)data.get_size(3);

        im_header.position[0] = acq_header.position[0];
        im_header.position[1] = acq_header.position[1];
        im_header.position[2] = acq_header.position[2];

        im_header.read_dir[0] = acq_header.read_dir[0];
        im_header.read_dir[1] = acq_header.read_dir[1];
        im_header.read_dir[2] = acq_header.read_dir[2];

        im_header.phase_dir[0] = acq_header.phase_dir[0];
        im_header.phase_dir[1] = acq_header.phase_dir[1];
        im_header.phase_dir[2] = acq_header.phase_dir[2];

        im_header.slice_dir[0] = acq_header.slice_dir[0];
        im_header.slice_dir[1] = acq_header.slice_dir[1];
        im_header.slice_dir[2] = acq_header.slice_dir[2];

        im_header.patient_table_position[0] = acq_header.patient_table_position[0];
        im_header.patient_table_position[1] = acq_header.patient_table_position[1];
        im_header.patient_table_position[2] = acq_header.patient_table_position[2];

        im_header.average    = acq_header.idx.average;
        im_header.slice      = acq_header.idx.slice;
        im_header.contrast   = acq_header.idx.contrast;
        im_header.phase      = acq_header.idx.phase;
        im_header.repetition = acq_header.idx.repetition;
        im_header.set        = acq_header.idx.set;

        im_header.acquisition_time_stamp = acq_header.acquisition_time_stamp;

        im_header.physiology_time_stamp[0] = acq_header.physiology_time_stamp[0];
        im_header.physiology_time_stamp[1] = acq_header.physiology_time_stamp[1];
        im_header.physiology_time_stamp[2] = acq_header.physiology_time_stamp[2];

        im_header.image_type         = ISMRMRD::ISMRMRD_IMTYPE_COMPLEX;
        im_header.image_series_index = 0;

        memcpy(im_header.user_int, acq_header.user_int, sizeof(int32_t) * ISMRMRD::ISMRMRD_USER_INTS);
        memcpy(im_header.user_float, acq_header.user_float, sizeof(float) * ISMRMRD::ISMRMRD_USER_FLOATS);

        im_header.attribute_string_len = 0;
        return im_header;
    }
}
