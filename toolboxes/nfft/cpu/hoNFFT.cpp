/** hoNFFT.cpp */

#include "hoNFFT.h"

#include "hoNDFFT.h"
#include "hoNDArray_elemwise.h"
#include "hoNDArray_reductions.h"
#include "hoNDArray_utils.h"

#include "vector_td_utilities.h"
#include "vector_td_operators.h"
#include "vector_td_io.h"

#include <algorithm>
#include <vector>
#include <cmath>
#include <stdexcept>

#include "hoNFFT_sparseMatrix.h"
#include "KaiserBessel_kernel.h"
#include <boost/range/algorithm/transform.hpp>

#include "GadgetronTimer.h"

#include "NFFT.hpp"
#include "NDArray_utils.h"

using namespace std;

namespace Gadgetron {


    namespace {

        template<typename T, unsigned int D>
        struct FFTD {
        };

        template<typename T>
        struct FFTD<T, 1> {
            using REAL = typename realType<T>::Type;

            static void fft(hoNDArray<T> &array, NFFT_fft_mode mode, bool do_scale) {
                if (mode == NFFT_fft_mode::FORWARDS) {
                    hoNDFFT<REAL>::instance()->fft1c(array);
                } else {
                    hoNDFFT<REAL>::instance()->ifft1c(array);
                }
                if (!do_scale) array *= std::sqrt(REAL(array.get_size(0)));
            }
        };

        template<typename T>
        struct FFTD<T, 2> {
            using REAL = typename realType<T>::Type;

            static void fft(hoNDArray<T> &array, NFFT_fft_mode mode, bool do_scale ) {
                if (mode == NFFT_fft_mode::FORWARDS) {
                    hoNDFFT<REAL>::instance()->fft2c(array);
                } else {
                    hoNDFFT<REAL>::instance()->ifft2c(array);
                }

                if (!do_scale) array *= std::sqrt(REAL(array.get_size(0)*array.get_size(1)));

            }
        };

        template<typename T>
        struct FFTD<T, 3> {
            using REAL = typename realType<T>::Type;

            static void fft(hoNDArray<T> &array, NFFT_fft_mode mode, bool do_scale) {
                if (mode == NFFT_fft_mode::FORWARDS) {
                    hoNDFFT<REAL>::instance()->fft3c(array);
                } else {
                    hoNDFFT<REAL>::instance()->ifft3c(array);
                }

                if (!do_scale) array *= std::sqrt(REAL(array.get_size(0)*array.get_size(1)*array.get_size(2)));
            }
        };

        template<class REAL> hoNDArray<std::complex<REAL>>
        compute_deapodization_filter(const vector_td<size_t,1>& image_dims, const vector_td<REAL,1>& beta, REAL W){

            hoNDArray<std::complex<REAL>> deapodization(to_std_vector(image_dims));
            vector_td<REAL,1> image_dims_real(image_dims);
            for (int x = 0; x < image_dims[0]; x++){
                auto offset = x - image_dims_real[0]/2;
                deapodization(x) = std::abs(offset) < W/2 ? KaiserBessel(offset,image_dims_real[0],REAL(1)/W,beta[0]) : REAL(0);
            }
            return deapodization;
        }
        template<class REAL> hoNDArray<std::complex<REAL>>
        compute_deapodization_filter(const vector_td<size_t,2>& image_dims, const vector_td<REAL,2>& beta, REAL W){

            hoNDArray<std::complex<REAL>> deapodization(to_std_vector(image_dims));
            vector_td<REAL,2> image_dims_real(image_dims);
            for (int y = 0; y < image_dims[1]; y++) {
                auto offset_y = y - image_dims_real[1]/2;
                auto weight_y = std::abs(offset_y) < W/2 ? KaiserBessel(offset_y,image_dims_real[1],REAL(1)/W,beta[1]) : REAL(0);

                for (int x = 0; x < image_dims[0]; x++) {
                    auto offset_x = x - image_dims_real[0]/2;
                    auto weight_x = std::abs(offset_x) < W/2 ? KaiserBessel(offset_x,image_dims_real[0],REAL(1)/W,beta[0]) : REAL(0);

                    deapodization(x,y) = weight_x*weight_y;
                }
            }
            return deapodization;
        }

         template<class REAL> hoNDArray<std::complex<REAL>>
        compute_deapodization_filter(const vector_td<size_t,3>& image_dims, const vector_td<REAL,3>& beta, REAL W){

            hoNDArray<std::complex<REAL>> deapodization(to_std_vector(image_dims));
            vector_td<REAL,3> image_dims_real(image_dims);
            for (int z = 0; z < image_dims[2]; z++) {
                auto offset_z = z - image_dims_real[2]/2;
                auto weight_z = std::abs(offset_z) < W/2 ? KaiserBessel(offset_z,image_dims_real[2],REAL(1)/W,beta[2]) : REAL(0);

                for (int y = 0; y < image_dims[1]; y++) {
                    auto offset_y = y - image_dims_real[1]/2;
                    auto weight_y = std::abs(offset_y) < W/2 ? KaiserBessel(offset_y,image_dims_real[1],REAL(1)/W,beta[1]) : REAL(0);

                    for (int x = 0; x < image_dims[0]; x++) {
                        auto offset_x = x - image_dims_real[0]/2;
                        auto weight_x = std::abs(offset_x) < W/2 ? KaiserBessel(offset_x,image_dims_real[0],REAL(1)/W,beta[0]) : REAL(0);

                        deapodization(x,y,z) = weight_x*weight_y*weight_z;
                    }
                }
            }
            return deapodization;
        }
    }



    template<class REAL, unsigned int D>
    hoNFFT_plan<REAL, D>::hoNFFT_plan(
            const vector_td<size_t, D> &matrix_size,
            const vector_td<size_t, D> &matrix_size_os,
            REAL W
    ) : NFFT_plan<hoNDArray,REAL,D>(matrix_size,matrix_size_os,W) {

        this->beta = compute_beta(W,matrix_size,matrix_size_os);
        this->deapodization_filter_IFFT = compute_deapodization_filter(this->matrix_size_os,this->beta, this->W);
        this->deapodization_filter_FFT = deapodization_filter_IFFT;
        FFTD<std::complex<REAL>,D>::fft(deapodization_filter_IFFT,NFFT_fft_mode::BACKWARDS,true);
        FFTD<std::complex<REAL>,D>::fft(deapodization_filter_FFT,NFFT_fft_mode::FORWARDS,true);

        boost::transform(deapodization_filter_IFFT,deapodization_filter_IFFT.begin(),[](auto val){return REAL(1)/val;});
        boost::transform(deapodization_filter_FFT,deapodization_filter_FFT.begin(),[](auto val){return REAL(1)/val;});
    }

    template<class REAL, unsigned int D>
    hoNFFT_plan<REAL, D>::hoNFFT_plan(const vector_td<size_t, D> &matrix_size, REAL oversampling_factor, REAL W)
            : NFFT_plan<hoNDArray, REAL, D>(matrix_size, oversampling_factor, W) {


        this->beta = compute_beta(W, matrix_size, this->matrix_size_os);
        this->deapodization_filter_IFFT = compute_deapodization_filter(this->matrix_size_os, this->beta, this->W);
        this->deapodization_filter_FFT = deapodization_filter_IFFT;

        FFTD<std::complex<REAL>, D>::fft(deapodization_filter_IFFT, NFFT_fft_mode::BACKWARDS, false);
        FFTD<std::complex<REAL>, D>::fft(deapodization_filter_FFT, NFFT_fft_mode::FORWARDS, false);
        
        boost::transform(deapodization_filter_IFFT, deapodization_filter_IFFT.begin(),
                         [](auto val) { return REAL(1) / val; });
        boost::transform(deapodization_filter_FFT, deapodization_filter_FFT.begin(),
                         [](auto val) { return REAL(1) / val; });
    }


    template<class REAL, unsigned int D>
    void hoNFFT_plan<REAL, D>::preprocess(
            const hoNDArray<vector_td<REAL, D>> &trajectories, NFFT_prep_mode mode) {

        NFFT_plan<hoNDArray,REAL,D>::preprocess(trajectories,mode);
        auto trajectories_scaled = trajectories;
        auto matrix_size_os_real = vector_td<REAL,D>(this->matrix_size_os);
        std::transform(trajectories_scaled.begin(),trajectories_scaled.end(),trajectories_scaled.begin(),[matrix_size_os_real](auto point){
           return (point+REAL(0.5))*matrix_size_os_real;
        });

        convolution_matrix.reserve(this->number_of_frames);
        convolution_matrix_T.reserve(this->number_of_frames);

        for (auto traj : NDArrayViewRange<hoNDArray<vector_td<REAL,D>>>(trajectories_scaled,0)){
            convolution_matrix.push_back(NFFT_internal::make_NFFT_matrix(traj, this->matrix_size_os, this->W, beta));
            if (mode == NFFT_prep_mode::ALL || mode == NFFT_prep_mode::NC2C) {
                convolution_matrix_T.push_back(NFFT_internal::transpose(convolution_matrix.back()));
            }
        }


    }

    template<class REAL, unsigned int D>
    void hoNFFT_plan<REAL, D>::compute(
            const hoNDArray<ComplexType> &d,
            hoNDArray<ComplexType> &m,
            const hoNDArray<REAL> *dcw,
            NFFT_comp_mode mode
    ) {
        const auto *pd = reinterpret_cast<const hoNDArray<complext<REAL>> *>(&d);
        auto *pm = reinterpret_cast<hoNDArray<complext<REAL>> *>(&m);

        this->compute(*pd, *pm, dcw, mode);
    }

    template<class REAL, unsigned int D>
    void hoNFFT_plan<REAL, D>::compute(
            const hoNDArray<complext<REAL>> &d,
            hoNDArray<complext<REAL>> &m,
            const hoNDArray<REAL> *dcw,
            NFFT_comp_mode mode
    ) {
       NFFT_plan<hoNDArray,REAL,D>::compute(d,m,dcw,mode);
    }



    template<class REAL, unsigned int D>
    void hoNFFT_plan<REAL, D>::mult_MH_M(
            const hoNDArray<complext<REAL>> &in,
            hoNDArray<complext<REAL>> &out,
            const hoNDArray<REAL>* dcw
    ) {
        const hoNDArray<ComplexType> *pin = reinterpret_cast<const hoNDArray<ComplexType> *>(&in);
        hoNDArray<ComplexType> *pout = reinterpret_cast<hoNDArray<ComplexType> *>(&out);

        this->mult_MH_M(*pin, *pout, dcw);
    }

    template<class REAL, unsigned int D>
    void hoNFFT_plan<REAL, D>::mult_MH_M(
            const hoNDArray<ComplexType> &in,
            hoNDArray<ComplexType> &out,
            const hoNDArray<REAL>* dcw
    ) {
        std::vector<size_t> dims = {this->number_of_samples,this->number_of_frames};
        auto batches = in.get_number_of_elements()/(prod(this->matrix_size)*this->number_of_frames);
        dims.push_back(batches);

        hoNDArray<ComplexType> tmp(dims);
        compute(in, tmp, dcw, NFFT_comp_mode::FORWARDS_C2NC);
        compute(tmp, out,dcw, NFFT_comp_mode::BACKWARDS_NC2C);
    }

    template<class REAL, unsigned int D>
    void hoNFFT_plan<REAL, D>::convolve(
            const hoNDArray<ComplexType> &d,
            hoNDArray<ComplexType> &m,
            NFFT_conv_mode mode,
            bool accumulate
    ) {
        if (mode == NFFT_conv_mode::NC2C)
            convolve_NFFT_NC2C(d, m,accumulate);
        else
            convolve_NFFT_C2NC(d, m, accumulate);
    }
        template<class REAL, unsigned int D>
    void hoNFFT_plan<REAL, D>::convolve(
            const hoNDArray<complext<REAL>> &in,
            hoNDArray<complext<REAL>> &out,
            NFFT_conv_mode mode,
            bool accumulate
    ) {

        const hoNDArray<ComplexType> *pin = reinterpret_cast<const hoNDArray<ComplexType> *>(&in);
        hoNDArray<ComplexType> *pout = reinterpret_cast<hoNDArray<ComplexType> *>(&out);
        this->convolve(*pin,*pout,mode,accumulate);

    }

    template<class REAL, unsigned int D>
    void hoNFFT_plan<REAL, D>::fft(
            hoNDArray<ComplexType> &d,
            NFFT_fft_mode mode,
            bool do_scale
    ) {
        FFTD<std::complex<REAL>, D>::fft(d, mode,do_scale);
    }
    template<class REAL, unsigned int D>
    void hoNFFT_plan<REAL, D>::fft(
            hoNDArray<complext<REAL>> &d,
            NFFT_fft_mode mode,
            bool do_scale
    ) {
        hoNDArray<ComplexType> *pd = reinterpret_cast<hoNDArray<ComplexType> *>(&d);
        this->fft(*pd,mode,do_scale);
    }

    template<class REAL, unsigned int D>
    void hoNFFT_plan<REAL, D>::deapodize(
            hoNDArray<ComplexType> &d,
            bool fourierDomain
    ) {
        if (fourierDomain){
            d *= deapodization_filter_FFT;
        } else {
            d *= deapodization_filter_IFFT;
        }
    }
        template<class REAL, unsigned int D>
    void hoNFFT_plan<REAL, D>::deapodize(
            hoNDArray<complext<REAL>> &d,
            bool fourierDomain
    ) {
        hoNDArray<ComplexType> *pd = reinterpret_cast<hoNDArray<ComplexType> *>(&d);
        this->deapodize(*pd,fourierDomain);
    }



    namespace {
        template<class REAL> void
        matrix_vector_multiply(const Gadgetron::NFFT_internal::NFFT_Matrix<REAL>& matrix, const complext<REAL>* vector, complext<REAL>* result) {

            for (size_t i = 0; i < matrix.n_cols; i++) {
                auto &row_indices = matrix.indices[i];
                auto &weights = matrix.weights[i];

#ifndef WIN32
    #pragma omp simd
#endif // WIN32
                for (size_t n = 0; n < row_indices.size(); n++) {
                    result[i] += vector[row_indices[n]] * weights[n];
                }
            }
        }

    }
    template<class REAL, unsigned int D>
    void hoNFFT_plan<REAL, D>::convolve_NFFT_C2NC(
            const hoNDArray<ComplexType> &cartesian,
            hoNDArray<ComplexType> &non_cartesian, bool accumulate
    ) {

        size_t nbatches = cartesian.get_number_of_elements()/convolution_matrix.front().n_rows;
        assert(nbatches == non_cartesian.get_number_of_elements()/convolution_matrix.front().n_cols);

        if (!accumulate) clear(&non_cartesian);

#pragma omp parallel for
        for (int b = 0; b < (int)nbatches; b++) {

            const ComplexType* cartesian_view = cartesian.get_data_ptr()+b*convolution_matrix.front().n_rows;
            ComplexType* non_cartesian_view = non_cartesian.get_data_ptr()+b*convolution_matrix.front().n_cols;
            size_t matrix_index = b%convolution_matrix.size();
            matrix_vector_multiply(convolution_matrix[matrix_index],(complext<REAL>*)cartesian_view,(complext<REAL>*)non_cartesian_view);
        }

    }

    template<class REAL, unsigned int D>
    void hoNFFT_plan<REAL, D>::convolve_NFFT_NC2C(
            const hoNDArray<ComplexType> &non_cartesian,
            hoNDArray<ComplexType> &cartesian, bool accumulate
    ) {
                size_t nbatches = cartesian.get_number_of_elements()/convolution_matrix.front().n_rows;
        assert(nbatches == non_cartesian.get_number_of_elements()/convolution_matrix.front().n_cols);
        if (!accumulate) clear(&cartesian);
#pragma omp parallel for
        for (int b = 0; b < (int)nbatches; b++) {

            ComplexType *cartesian_view = cartesian.get_data_ptr() + b * convolution_matrix.front().n_rows;
            const ComplexType *non_cartesian_view = non_cartesian.get_data_ptr() + b * convolution_matrix.front().n_cols;
            size_t matrix_index = b%convolution_matrix.size();
            matrix_vector_multiply(convolution_matrix_T[matrix_index], (complext<REAL>*)non_cartesian_view, (complext<REAL>*)cartesian_view);

        }
    }


    template<class REAL, unsigned int D>
    vector_td<REAL, D> hoNFFT_plan<REAL, D>::compute_beta(REAL W, const Gadgetron::vector_td<size_t, D> &matrix_size,
                                                          const Gadgetron::vector_td<size_t, D> &matrix_size_os) {
        // Compute Kaiser-Bessel beta paramter according to the formula provided in
        // Beatty et. al. IEEE TMI 2005;24(6):799-808.
        using boost::math::constants::pi;
        vector_td<REAL, D> beta;

        auto alpha = matrix_size_os / matrix_size;

        for (int d = 0; d < D; d++) {
            beta[d] = (pi<REAL>() * std::sqrt(
                    (W * W) / (alpha[d] * alpha[d]) * (alpha[d] - REAL(0.5)) * (alpha[d] - REAL(0.5)) -
                    REAL(0.8)));
        }

        return beta;
    }

    template<class REAL, unsigned int D>
    boost::shared_ptr<hoNFFT_plan<REAL,D>> NFFT<hoNDArray,REAL,D>::make_plan(const Gadgetron::vector_td<size_t, D> &matrix_size,
                                        const Gadgetron::vector_td<size_t, D> &matrix_size_os, REAL W) {
        return boost::make_shared<hoNFFT_plan<REAL,D>>(matrix_size,matrix_size_os,W);
    }



}

template
class Gadgetron::hoNFFT_plan<float, 1>;

template
class Gadgetron::hoNFFT_plan<float, 2>;

template
class Gadgetron::hoNFFT_plan<float, 3>;

template
class Gadgetron::hoNFFT_plan<double, 1>;

template
class Gadgetron::hoNFFT_plan<double, 2>;

template
class Gadgetron::hoNFFT_plan<double, 3>;

template class Gadgetron::NFFT<Gadgetron::hoNDArray,float,1>;
template class Gadgetron::NFFT<Gadgetron::hoNDArray,float,2>;
template class Gadgetron::NFFT<Gadgetron::hoNDArray,float,3>;



template class Gadgetron::NFFT<Gadgetron::hoNDArray,double,1>;
template class Gadgetron::NFFT<Gadgetron::hoNDArray,double,2>;
template class Gadgetron::NFFT<Gadgetron::hoNDArray,double,3>;
