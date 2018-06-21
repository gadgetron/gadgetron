//
// Created by dchansen on 6/20/18.
//

#include "non_local_bayes.h"
#include "hoNDArray.h"
#include "vector_td_utilities.h"
#include <Eigen/Dense>
#include <GadgetronTimer.h>


namespace Gadgetron {
    namespace Denoise {


        namespace {

            template<class T, int D>
            Eigen::Matrix<T, D * D, 1> get_patch(const hoNDArray<T> &image, int x, int y,
                                                 const vector_td<int, 2> &image_dims) {

                constexpr int N = D * D;
                Eigen::Matrix<T, D * D, 1> window;
                for (int ky = 0; ky < D; ky++) {
                    for (int kx = 0; kx < D; kx++) {
                        window[kx + ky * D] = image(((kx - D / 2) + x + image_dims[0]) % image_dims[0],
                                                    ((ky - D / 2) + y + image_dims[1]) % image_dims[1]);

                    }
                }
                return window;
            };


            template<class T, int D>
            hoNDArray<vector_td<T, D * D>> create_patches(const hoNDArray<T> &image) {
                const vector_td<int, 2> image_dims = vector_td<int, 2>(
                        from_std_vector<size_t, 2>(*image.get_dimensions()));

                hoNDArray<vector_td<T, D * D>> result(image.get_dimensions());

                for (int ky = 0; ky < image.get_size(0); ky++) {
                    for (int kx = 0; kx < image.get_size(0); kx++) {
                        result(kx, ky) = get_patch(image, kx, ky, image_dims);
                    }
                }

                return result;

            };

//    template<class T, unsigned int D>
//    std::vector<float>
//    calculate_distances(const hoNDArray<T> &image, float noise_std, int search_radius,
//                        const vector_td<int, 2> &image_dims,
//                        int ky, int kx) {
//        std::vector<float> distances(search_radius * search_radius);
//
//        auto current_patch = get_patch(image, kx, ky, image_dims);
//        for (int dy = -search_radius; dy < search_radius; dy++) {
//            for (int dx = -search_radius; dx < search_radius; dx++) {
//
//                distances[dx + dy * search_radius] = distance(current_patch,
//                                                              get_patch(image, kx + dx, ky + dy, image_dims),
//                                                              noise_std);
//
//            }
//        }
//        return distances;
//    }

            template<class T, unsigned int D>
            Eigen::Matrix<T, D * D, 1>
            get_mean_patch(const hoNDArray<T> &image, int search_radius, const vector_td<int, 2> &image_dims,
                           int kx,
                           int ky) {


                Eigen::Matrix<T, D * D, 1> mean_patch = Eigen::Matrix<T, D * D, 1>::Constant(T(0));

                for (int dy = -search_radius / 2; dy < search_radius / 2; dy++) {
                    for (int dx = -search_radius / 2; dx < search_radius / 2; dx++) {
                        auto other_patch = get_patch<T, D>(image, kx + dx, ky + dy, image_dims);
                        mean_patch += other_patch;
                    }
                }
                mean_patch /= search_radius * search_radius;
                return mean_patch;
            }

            template<class T, unsigned int D>
            Eigen::Matrix<T, D * D, D * D>
            get_covariance_matrix(const hoNDArray<T> &image, int search_radius, Eigen::Matrix<T, D * D, 1> mean_patch,
                                  const vector_td<int, 2> &image_dims, int kx,
                                  int ky) {

                using CovMat = Eigen::Matrix<T, D * D, D * D>;
                CovMat covariance_matrix = CovMat::Constant(T(0));

                for (int dy = -search_radius / 2; dy < search_radius / 2; dy++) {
                    for (int dx = -search_radius / 2; dx < search_radius / 2; dx++) {
                        auto other_patch = get_patch<T, D>(image, kx + dx, ky + dy, image_dims);
                        covariance_matrix += (other_patch - mean_patch) * (other_patch - mean_patch).transpose();
                    }
                }
                covariance_matrix /= search_radius * search_radius - 1;
                return covariance_matrix;
            }


            template<class T>
            hoNDArray<T> non_local_bayes_single_image(const hoNDArray<T> &image, float noise_std, int search_radius) {

                constexpr int D = 5;

                hoNDArray<T> result(image.get_dimensions());
                result.fill(0);

                const float noise_std2 = noise_std * noise_std;
                const vector_td<int, 2> image_dims = vector_td<int, 2>(
                        from_std_vector<size_t, 2>(*image.get_dimensions()));

//                auto patches = create_patches<T,D>(image);


//#pragma omp parallel for collapse(2)
                for (int ky = 0; ky < image.get_size(1); ky++) {
                    for (int kx = 0; kx < image.get_size(0); kx++) {

                        auto mean_patch = get_mean_patch<T, D>(image, search_radius, image_dims, kx, ky);
                        auto covariance_matrix = get_covariance_matrix<T, D>(image, search_radius, mean_patch,
                                                                             image_dims,
                                                                             kx,
                                                                             ky);


                        auto window = get_patch<T, D>(image, kx, ky, image_dims);
//                        window = mean_patch;
                        window =
                                 (covariance_matrix - noise_std2 * decltype(covariance_matrix)::Ones()) *
                                 covariance_matrix.inverse() * (window - mean_patch);

                        result(kx, ky) = window[D / 2 + (D / 2) * D];

                    }
                }

                return result;

            }


            template<class T>
            hoNDArray<T> non_local_bayes_T(const hoNDArray<T> &image, float noise_std, unsigned int search_radius) {

                GadgetronTimer("Non local bayes");
                size_t n_images = image.get_number_of_elements() / (image.get_size(0) * image.get_size(1));

                std::vector<size_t> image_dims = {image.get_size(0), image.get_size(1)};
                size_t image_elements = image_dims[0] * image_dims[1];

                auto result = hoNDArray<T>(image.get_dimensions());
                for (int i = 0; i < n_images; i++) {

                    auto image_view = hoNDArray<T>(image_dims, image.get_data_ptr() + i * image_elements);
                    auto result_view = hoNDArray<T>(image_dims, result.get_data_ptr() + i * image_elements);

                    result_view = non_local_bayes_single_image(image_view, noise_std, search_radius);


                }
                return result;
            }
        }


        hoNDArray<float> non_local_bayes(const hoNDArray<float> &image, float noise_std, unsigned int search_radius) {
            return non_local_bayes_T(image, noise_std, search_radius);
        }

        hoNDArray<std::complex<float>>
        non_local_bayes(const hoNDArray<std::complex<float>> &image, float noise_std, unsigned int search_radius) {
            return non_local_bayes_T(image, noise_std, search_radius);
        }
    }

}