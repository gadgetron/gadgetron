//
// Created by dchansen on 9/19/18.
//

#include <cpu/hoNDArray_utils.h>
#include "hoNDFFT.h"
#include "mri_core_girf_correction.h"
#include "hoNDArray_utils.h"
#include "hoNDArray_math.h"
#include <boost/filesystem.hpp>
#include <cpu/hoNDArray_fileio.h>
#include <boost/math/constants/constants.hpp>

#include <boost/config/warning_disable.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_stl.hpp>
#include <boost/fusion/include/adapt_adt.hpp>

using namespace std::complex_literals;

BOOST_FUSION_ADAPT_ADT(
            std::complex<float>,
    (float, float, obj.real(),obj.real(val))
    (float,float,obj.imag(),obj.imag(val))
)
namespace Gadgetron {
    namespace {
        constexpr float PI = boost::math::constants::pi<float>();
        template<unsigned int D>
        hoNDArray<vector_td<float,D>> girf_correct_D(const hoNDArray<vector_td<float,D>> &gradients,
                                        const hoNDArray<std::complex<float>> &girf_kernel,
                                        const arma::fmat33 &rotation_matrix, float gradient_sampling_time,
                                        float girf_sampling_time, float TE) {

            size_t nlines = gradients.get_size(0);
//            size_t nlines_new = size_t(nlines * girf_sampling_time / gradient_sampling_time);
            size_t nlines_new = girf_kernel.get_size(0);

            size_t nbatches = gradients.get_size(1);

            hoNDArray<vector_td<float,D>> result(gradients.get_dimensions());

            for (int batch = 0; batch < nbatches; batch++) {

                hoNDArray<std::complex<float>> rotated(nlines_new, 3);
                std::fill(rotated.begin(), rotated.end(), 0);


                for (size_t i = 0; i < nlines; i++) {
                    arma::fvec3 gradient_vector(arma::fill::zeros);
                    for (size_t k = 0; k < D; k++)
                        gradient_vector[k] = gradients(i,batch)[k];

                    gradient_vector = rotation_matrix * gradient_vector;
                    for (size_t k = 0; k < 3; k++)
                        rotated(i, k) = gradient_vector[k];
                }


                hoNDFFT<float>::instance()->fft1(rotated);
                hoNDFFT<float>::instance()->fftshift1D(rotated);
//            auto padded_kernel = std::move(*pad<std::complex<float>,1>(uint64d1(nlines_new), &girf_kernel));

//            rotated *= padded_kernel;
            rotated *= girf_kernel;
            float ADCShift = -TE-0.5*gradient_sampling_time;
                float frequency_resolution = 1.0 / (girf_sampling_time * nlines_new);
            for (size_t k = 0; k < 3; k++) {
                for (size_t i = 0; i < nlines_new; i++){
                    rotated(i,k) *= std::exp(-2if*ADCShift*PI*frequency_resolution*float(i));
                }
            }

                size_t full_domain = size_t(std::round(1.0 / (2 * frequency_resolution * gradient_sampling_time))) * 2;

//                rotated = *pad<std::complex<float>, 1>(uint64d1(full_domain), &rotated);

                hoNDFFT<float>::instance()->ifftshift1D(rotated);
                hoNDFFT<float>::instance()->ifft1(rotated);


//            hoNDFFT<float>::instance()->ifft(&rotated,0);

                arma::fmat33 inv_rotation = arma::inv(rotation_matrix);

                for (size_t i = 0; i < nlines; i++) {
                    // arma::cx_vec3 gradient_vector = {rotated(i, 0), rotated(i, 1), rotated(i, 2)};
                    arma::cx_vec3 gradient_vector;
                    gradient_vector[0] = rotated(i, 0);
                    gradient_vector[1] = rotated(i, 1);
                    gradient_vector[2] = rotated(i, 2);
                    gradient_vector = inv_rotation * gradient_vector;
                    for (size_t k = 0; k < D; k++) {
                        result(i,batch)[k] = std::copysign(std::abs(gradient_vector[k]), std::real(gradient_vector[k]));
                    }
                };
            }

            return result;
        }
    }


    namespace {
        template<typename Iterator>
        bool parse_girf_text(Iterator first, Iterator last, std::vector<std::vector<std::complex<float>>> &v) {
            using namespace boost::spirit::qi;
            namespace qi = boost::spirit::qi;
            using boost::spirit::ascii::space;

            rule<Iterator, std::complex<float>(), ascii::space_type> base = float_ >> ((float_ >> (lit('i') | lit('j'))) | attr(0.0f));
            decltype(base) complex = (lit('(') >> base >> lit(')')) | base;
            rule<Iterator, std::vector<std::complex<float>>(), ascii::space_type> row = (lit('[') >> complex % ','
                                                                                                  >> lit(']'));
            bool r = phrase_parse(first, last,
                                  lit('[') >> row >> qi::repeat(2)[lit(',') >> row] >> lit(']'),
                                  space, v);
            if (!r || first != last) // fail if we did not get a full match
                return false;
            return r;
        }
    }



    hoNDArray<std::complex<float>> GIRF::load_girf_kernel(const std::string& girf_string) {
        std::vector<std::vector<std::complex<float>>> temp_array;

        bool success = parse_girf_text(girf_string.begin(),girf_string.end(),temp_array);
        if (!success) throw std::runtime_error("Failed to parse girf string");

        hoNDArray<std::complex<float>> result(temp_array[0].size(),temp_array.size());
        auto data_ptr = result.get_data_ptr();
        for (auto& v: temp_array)
            for (auto c: v)
                *(data_ptr++) = c;


        return result;

    }

    hoNDArray<floatd3>
    GIRF::girf_correct(const hoNDArray<floatd3> &gradients, const hoNDArray<std::complex<float>> &girf_kernel,
                       const arma::fmat33 &rotation_matrix, float gradient_sampling_time, float girf_sampling_time,
                       float TE) {
        return girf_correct_D(gradients,girf_kernel,rotation_matrix,girf_sampling_time,girf_sampling_time,TE);
    }

    hoNDArray<floatd2>
    GIRF::girf_correct(const hoNDArray<floatd2> &gradients, const hoNDArray<std::complex<float>> &girf_kernel,
                       const arma::fmat33 &rotation_matrix, float gradient_sampling_time, float girf_sampling_time,
                       float TE) {
        return girf_correct_D(gradients,girf_kernel,rotation_matrix,girf_sampling_time,girf_sampling_time,TE);
    }


}
