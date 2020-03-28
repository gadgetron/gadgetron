//
// Created by dchansen on 6/19/18.
//

#include <vector_td_utilities.h>
#include <GadgetronTimer.h>
#include "non_local_means.h"


namespace Gadgetron {
    namespace Denoise {

        namespace {


            template<class T, int D>
            vector_td<T, D * D> get_patch(const hoNDArray<T> &image, int x, int y,
                                          const vector_td<int, 2> &image_dims) {

                constexpr int N = D * D;
                vector_td<T, N> window;
                for (int ky = 0; ky < D; ky++) {
                    for (int kx = 0; kx < D; kx++) {
                        window[kx + ky * D] = image(((kx - D / 2) + x + image_dims[0]) % image_dims[0],
                                                    ((ky - D / 2) + y + image_dims[1]) % image_dims[1]);

                    }
                }
                return window;
            };


            template<class T>
            hoNDArray<T> non_local_means_single_image(const hoNDArray<T> &image, float noise_std, int search_radius) {

                constexpr int D = 5;

                hoNDArray<T> result(image.dimensions());
                const float noise_std2 = noise_std * noise_std;
                const vector_td<int, 2> image_dims = vector_td<int, 2>(
                        from_std_vector<size_t, 2>(image.dimensions()));

#pragma omp parallel for
                for (int ky = 0; ky < image.get_size(1); ky++) {
                    for (int kx = 0; kx < image.get_size(0); kx++) {
                        float sum_weight = 0;
                        T sum_value = 0;
                        auto window = get_patch<T, D>(image, kx, ky, image_dims);

                        for (int dy = -int(search_radius); dy < search_radius; dy++) {
                            for (int dx = -int(search_radius); dx < search_radius; dx++) {

                                auto window2 = get_patch<T, D>(image, kx + dx, ky + dy, image_dims);
                                auto diff = window - window2;
                                auto weight = std::exp(-norm_squared(diff) / (noise_std2 * D * D));

                                sum_weight += weight;
                                sum_value += weight * window2[D / 2 + D * (D / 2)];

                            }
                        }

                        result(kx, ky) = sum_value / sum_weight;

                    }
                }

                return result;

            }

            template<class T>
            hoNDArray<T> non_local_means_T(const hoNDArray<T> &image, float noise_std, unsigned int search_radius) {


                GadgetronTimer timer("Non local means");
                             size_t n_images = image.get_number_of_elements() / (image.get_size(0) * image.get_size(1));

                std::vector<size_t> image_dims = {image.get_size(0), image.get_size(1)};
                size_t image_elements = image_dims[0] * image_dims[1];

                auto result = hoNDArray<T>(image.dimensions());

                #pragma omp parallel for
                for (int i = 0; i < n_images; i++) {

                    auto image_view = hoNDArray<T>(image_dims, const_cast<T*>(image.get_data_ptr() + i * image_elements));
                    auto result_view = non_local_means_single_image(image_view, noise_std,search_radius) ;

                    memcpy(result.begin() + i * image_elements, result_view.begin(), result_view.get_number_of_bytes());
                    GDEBUG_STREAM("I finished one");
                }
                return result;

            }
        }



        hoNDArray<float> non_local_means(const hoNDArray<float> &image, float noise_std, unsigned int search_radius) {
            return non_local_means_T(image, noise_std, search_radius);
        }

        hoNDArray<std::complex<float>>
        non_local_means(const hoNDArray<std::complex<float>> &image, float noise_std, unsigned int search_radius) {
            return non_local_means_T(image, noise_std, search_radius);
        }


    }
}
