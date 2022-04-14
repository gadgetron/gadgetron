//
// Created by dchansen on 9/19/18.
// Edited by ajaved on  5/20/2020 - Fixed bugs and matched the performance with results achieved using ACW GIRF in MATLAB
//
#include "mri_core_girf_correction.h"

#include <boost/filesystem.hpp>
#include <boost/config/warning_disable.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_stl.hpp>
#include <boost/fusion/include/adapt_adt.hpp>
#include "mri_core_partial_fourier.h"
#include "mri_core_kspace_filter.h"
using namespace std::complex_literals;
using namespace Gadgetron;

BOOST_FUSION_ADAPT_ADT(
    std::complex<float>,
    (float, float, obj.real(), obj.real(val))(float, float, obj.imag(), obj.imag(val)))
namespace Gadgetron
{
    namespace GIRF
    {
        constexpr float PI = boost::math::constants::pi<float>();
        template <unsigned int D>
        hoNDArray<vector_td<float, D>> girf_correct_D(const hoNDArray<vector_td<float, D>> &gradients,
                                                      const hoNDArray<std::complex<float>> &girf_kernel,
                                                      const arma::fmat33 &rotation_matrix, float gradient_sampling_time,
                                                      float girf_sampling_time, float TE) // Times are in us
        {

            using namespace std;

            size_t nlines = gradients.get_size(0);
            size_t nlines_new = girf_kernel.get_size(0) * std::round(girf_sampling_time / gradient_sampling_time); //std::max(girf_kernel.get_size(0) * girf_sampling_time / gradient_sampling_time,nlines * std::round(girf_sampling_time / gradient_sampling_time)); // This logic here may need to be checked
            size_t nbatches = gradients.get_size(1);

            //auto padded_kernel=GIRF::zeropadding(girf_kernel,int(girf_sampling_time/gradient_sampling_time));

            auto padded_kernel = pad<std::complex<float>, 1>(uint64d1(nlines_new), girf_kernel, 0);
            //auto padded_kernel = GIRF::zeropadding(girf_kernel, 5);
            auto dims = gradients.get_dimensions();
            auto result = hoNDArray<vector_td<float, D>>(gradients.get_dimensions());

            hoNDArray<float> filter;
            using namespace Gadgetron::Indexing;

            // Generate Hanning window filter (using ACW parameters) and normalize the max to 1
            Gadgetron::generate_symmetric_filter(400, filter, ISMRMRD_FILTER_HANNING, 1.5, 400);
            auto fmax = max(filter);
            filter /= fmax;

            hoNDArray<std::complex<float>> rotated(nlines_new, 3);

            for (int batch = 0; batch < nbatches; batch++)
            {
                std::fill(rotated.begin(), rotated.end(), 0);
                for (size_t i = 0; i < nlines; i++)
                {
                    arma::fvec3 gradient_vector(arma::fill::zeros);
                    for (size_t k = 0; k < D; k++)
                        gradient_vector[k] = gradients(i, batch)[k];

                    gradient_vector = rotation_matrix * gradient_vector;
                    for (size_t k = 0; k < 3; k++)
                        rotated(i, k) = gradient_vector[k];
                }
                float maxTx;
                float minTx;
                auto temp = real(permute(rotated, {0, 1}));
                Gadgetron::maxValue(hoNDArray<float>(temp(slice, 0)), maxTx);
                Gadgetron::minValue(hoNDArray<float>(temp(slice, 0)), minTx);

                if (maxTx > 100 || minTx < -100)
                    GERROR("What the heck");
                // Add filter to smooth out the transition for the gradients after zero-padding
                // This reduces ringing !
                for (size_t k = 0; k < 3; k++)
                {
                    for (size_t i = nlines; i < nlines + filter.size() / 2; i++)
                    {
                        rotated(i, k) = filter[filter.size() / 2 + i - nlines] * (rotated(nlines - 1, k)); // last term for fft scaling
                        if (real(rotated(i, k)) > 100 || real(rotated(i, k)) < -100)
                        {
                            GERROR("What the heck");
                        }
                    }
                    std::rotate(rotated.begin() + k * rotated.get_size(0), rotated.begin() + k * rotated.get_size(0) + (nlines / 2 + nlines_new / 2),
                                rotated.begin() + (k + 1) * rotated.get_size(0));
                }

                temp = real(permute(rotated, {0, 1}));
                Gadgetron::maxValue(hoNDArray<float>(temp(slice, 0)), maxTx);
                Gadgetron::minValue(hoNDArray<float>(temp(slice, 0)), minTx);

                if (maxTx > 100 || minTx < -100)
                    GERROR("What the heck");

                hoNDFFT<float>::instance()->fftshift1D(rotated);
                hoNDFFT<float>::instance()->fft1(rotated);
                hoNDFFT<float>::instance()->fftshift1D(rotated);

                rotated *= padded_kernel;

                float ADCShift = 1.0*(TE + 0.5 * gradient_sampling_time); //Added -ve sign inspired by NAM
                float frequency_resolution = 1.0 / (gradient_sampling_time * nlines_new);
                for (size_t k = 0; k < 3; k++)
                {
                    for (size_t i = 0; i < nlines_new; i++)
                    {
                        // Removed -ve sign needs more testing but was trying to match MATLAB implementation
                        rotated(i, k) = rotated(i, k) * std::exp(2if * ADCShift * PI * frequency_resolution * float(i));
                    }
                }

                size_t full_domain = size_t(std::round(1.0 / (2 * frequency_resolution * gradient_sampling_time))) * 2;

                // rotated = *pad<std::complex<float>, 1>(uint64d1(full_domain), &rotated);

                hoNDFFT<float>::instance()->ifftshift1D(rotated);
                hoNDFFT<float>::instance()->ifft1(rotated);
                hoNDFFT<float>::instance()->ifftshift1D(rotated);

                temp = real(permute(rotated, {0, 1}));
                Gadgetron::maxValue(hoNDArray<float>(temp(slice, 0)), maxTx);
                Gadgetron::minValue(hoNDArray<float>(temp(slice, 0)), minTx);

                if (maxTx > 100 || minTx < -100)
                    GERROR("What the heck");

                for (size_t k = 0; k < 3; k++)
                {

                    using namespace std;
                    std::rotate(rotated.begin() + k * rotated.get_size(0), rotated.begin() + k * rotated.get_size(0) + nlines_new - (nlines / 2 + nlines_new / 2) + 1,
                                rotated.begin() + (k + 1) * rotated.get_size(0));
                }

                arma::fmat33 inv_rotation = arma::inv(rotation_matrix);

                for (size_t i = 0; i < nlines; i++)
                {
                    // arma::cx_vec3 gradient_vector = {rotated(i, 0), rotated(i, 1), rotated(i, 2)};
                    arma::cx_vec3 gradient_vector;
                    gradient_vector[0] = rotated(i, 0);
                    gradient_vector[1] = rotated(i, 1);
                    gradient_vector[2] = rotated(i, 2);
                    gradient_vector = inv_rotation * gradient_vector;
                    for (size_t k = 0; k < D; k++)
                    {
                        result(i, batch)[k] = std::copysign(std::abs(gradient_vector[k]), std::real(gradient_vector[k]));
                    }
                };
            }

            return result;
        }

        template <typename Iterator>
        bool parse_girf_text(Iterator first, Iterator last, std::vector<std::vector<std::complex<float>>> &v)
        {
            using namespace boost::spirit::qi;
            namespace qi = boost::spirit::qi;
            using boost::spirit::ascii::space;

            rule<Iterator, std::complex<float>(), ascii::space_type> base = float_ >> ((float_ >> (lit('i') | lit('j'))) | attr(0.0f));
            decltype(base) complex = (lit('(') >> base >> lit(')')) | base;
            rule<Iterator, std::vector<std::complex<float>>(), ascii::space_type> row = (lit('[') >> complex % ',' >> lit(']'));
            bool r = phrase_parse(first, last,
                                  lit('[') >> row >> qi::repeat(2)[lit(',') >> row] >> lit(']'),
                                  space, v);
            if (!r || first != last) // fail if we did not get a full match
                return false;
            return r;
        }

        hoNDArray<std::complex<float>> load_girf_kernel(const std::string &girf_string)
        {
            std::vector<std::vector<std::complex<float>>> temp_array;

            bool success = Gadgetron::GIRF::parse_girf_text(girf_string.begin(), girf_string.end(), temp_array);
            if (!success)
                throw std::runtime_error("Failed to parse girf string");

            hoNDArray<std::complex<float>> result(temp_array[0].size(), temp_array.size());
            auto data_ptr = result.get_data_ptr();
            for (auto &v : temp_array)
                for (auto c : v)
                    *(data_ptr++) = c;

            return result;
        }

        hoNDArray<floatd3>
        girf_correct(const hoNDArray<floatd3> &gradients, const hoNDArray<std::complex<float>> &girf_kernel,
                                          const arma::fmat33 &rotation_matrix, float gradient_sampling_time, float girf_sampling_time,
                                          float TE)
        {
            return girf_correct_D(gradients, girf_kernel, rotation_matrix, girf_sampling_time, girf_sampling_time, TE);
        }

        hoNDArray<floatd2>
        girf_correct(const hoNDArray<floatd2> &gradients, const hoNDArray<std::complex<float>> &girf_kernel,
                                          const arma::fmat33 &rotation_matrix, float gradient_sampling_time, float girf_sampling_time,
                                          float TE)
        {
            return girf_correct_D(gradients, girf_kernel, rotation_matrix, gradient_sampling_time, girf_sampling_time, TE);
        }

        hoNDArray<std::complex<float>> zeropadding(hoNDArray<std::complex<float>> input, int zpadFactor)
        {
            using namespace Gadgetron::Indexing;
            hoNDArray<std::complex<float>> output(input.get_size(0) * zpadFactor, input.get_size(1));
            std::fill(output.begin(), output.end(), 0);
            for (int ii = 0; ii < output.get_size(0); ii++)
            {
                if (ii + 1 > output.get_size(0) / 2 - input.get_size(0) / 2 && ii < output.get_size(0) / 2 + (input.get_size(0) / 2))
                {
                    output(ii, slice) = input(ii - (output.get_size(0) / 2 - input.get_size(0) / 2), slice);
                    //   GDEBUG("output [%d]: value: %0.2f + i %0.2f\n", ii, real(output[ii]), imag(output[ii]));
                }
            }
            return output;
        }
        hoNDArray<std::complex<float>> readGIRFKernel(std::string folder)
        {
            using namespace std;
            using namespace boost::filesystem;
            using namespace Gadgetron::Indexing;
            path fnamex = folder + "GIRFx.txt";
            path fnamey = folder + "GIRFy.txt";
            path fnamez = folder + "GIRFz.txt";

            boost::filesystem::fstream filex, filey, filez;
            filex.open(fnamex);
            filey.open(fnamey);
            filez.open(fnamez);

            if (!filex.is_open() || !filey.is_open() || !filez.is_open())
                throw std::runtime_error("GIRF Kernel is not loaded and Waveform to Trajectories need this");
            vector<std::complex<float>> girfx;
            vector<std::complex<float>> girfy;
            vector<std::complex<float>> girfz;
            string temp_line;
            int girf_numpoint;
            float girf_sampletime;
            // Header first 4 lines
            for (int ii = 0; ii < 4; ii++)
            {
                getline(filex, temp_line);
                if (ii == 0)
                    girf_numpoint = stoi(temp_line);
                if (ii == 2)
                    girf_sampletime = stod(temp_line);

                getline(filey, temp_line); // got info from x no need for y and z but still need to skip the lines
                getline(filez, temp_line);
            }

            hoNDArray<std::complex<float>> gk(girf_numpoint, 3);
            string temp_liner;
            string temp_linei;
            int index = 0;
            while (getline(filex, temp_liner))
            {
                getline(filex, temp_linei);
                gk(index, 1) = complex<float>(stod(temp_liner), stod(temp_linei));

                getline(filey, temp_liner);
                getline(filey, temp_linei);
                gk(index, 0) = complex<float>(stod(temp_liner), stod(temp_linei));

                getline(filez, temp_liner);
                getline(filez, temp_linei);
                gk(index, 2) = complex<float>(stod(temp_liner), stod(temp_linei));
                index++;
            }
            return gk;
        }
    } // namespace GIRF
} // namespace Gadgetron