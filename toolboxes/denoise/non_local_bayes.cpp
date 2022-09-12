#include "non_local_bayes.h"
#include "hoNDArray.h"
#include "vector_td_utilities.h"
#include <GadgetronTimer.h>
#include "hoArmadillo.h"
#include <numeric>

namespace Gadgetron {
    namespace Denoise {

        namespace {

            template<class T>
            arma::Col<T> get_patch(const hoNDArray<T> &image, int x, int y, int patch_size,
                                   const vector_td<int, 2> &image_dims) {

                const int N = patch_size * patch_size;
                arma::Col<T> window = arma::Col<T>(N);
                for (int ky = 0; ky < patch_size; ky++) {
                    for (int kx = 0; kx < patch_size; kx++) {
                        window[kx + ky * patch_size] = image(
                                ((kx - patch_size / 2) + x + image_dims[0]) % image_dims[0],
                                ((ky - patch_size / 2) + y + image_dims[1]) % image_dims[1]);

                    }
                }
                return window;
            };



            template<class T>
            struct ImagePatch {
                arma::Col<T> patch;
                int center_x, center_y;
            };


            template<class T>
            std::vector<ImagePatch<T>>
            create_patches(const hoNDArray<T> &image, int kx, int ky, int patch_size, int search_window,
                           const vector_td<int, 2> &image_dims) {

                std::vector<ImagePatch<T>> result;
                result.reserve(search_window);

                for (int dy = std::max(ky - search_window / 2, 0);
                     dy < std::min(search_window / 2 + ky, image_dims[1]); dy++) {
                    for (int dx = std::max(kx - search_window / 2, 0);
                         dx < std::min(search_window / 2 + kx, image_dims[0]); dx++) {

                        result.push_back(ImagePatch<T>{get_patch(image, dx, dy, patch_size, image_dims), dx, dy});
                    }
                }

                return result;

            };


            template<class T>
            void add_patch(ImagePatch<T> &patch, hoNDArray<T> &image, hoNDArray<int> &count, int patch_size,
                           const vector_td<int, 2> &image_dims) {


                for (int ky = 0; ky < patch_size; ky++) {
                    auto output_ky = (patch.center_y + ky - patch_size / 2 + image_dims[1]) % image_dims[1];
                    for (int kx = 0; kx < patch_size; kx++) {
                        auto output_kx = (patch.center_x + kx - patch_size / 2 + image_dims[0]) % image_dims[0];
                        image(output_kx, output_ky) += patch.patch[kx + ky * patch_size];
                        count(output_kx, output_ky)++;
                    }
                }

            };


            template<class T>
            float distance(const arma::Col<T> &patch1, const arma::Col<T> &patch2) {
                arma::Col<T> diff = patch1 - patch2;

                float result = 0;
                for (auto d : diff) result += std::norm(d);

                result /= patch1.size() * patch1.size();
                return result;
            };

            template<class T>
            arma::Col<T> get_mean_patch(const std::vector<ImagePatch<T>> &patches) {

                auto mean_patch = arma::Col<T>(patches.front().patch.size(), arma::fill::zeros);

                for (auto &patch : patches) {
                    mean_patch += patch.patch;
                }

                mean_patch /= patches.size();

                return mean_patch;
            }

            template<class T>
            arma::Mat<T>
            get_covariance_matrix(const std::vector<ImagePatch<T>> &patches, const arma::Col<T> &mean_patch) {

                auto covariance_matrix = arma::Mat<T>(mean_patch.size(), mean_patch.size(), arma::fill::zeros);

                for (auto &patch : patches) {
                    covariance_matrix += (patch.patch - mean_patch) * (patch.patch - mean_patch).t();
                }
                covariance_matrix /= patches.size() - 1;

                return covariance_matrix;
            }


            template<class T>
            void
            filter_patches(std::vector<ImagePatch<T>> &patches, int max_n_patches, const arma::Col<T> &reference_patch) {

                auto distances = std::vector<float>(patches.size());
                transform(patches.begin(), patches.end(), distances.begin(),
                          [&](auto patch) { return distance(patch.patch, reference_patch); });

                std::vector<size_t> patch_indices(patches.size());
                std::iota(patch_indices.begin(), patch_indices.end(), 0);

                sort(patch_indices.begin(), patch_indices.end(),
                     [&](auto v1, auto v2) { return distances[v1] < distances[v2]; });

                int n_patches = std::min<int>(patches.size(),max_n_patches);
                std::vector<ImagePatch<T>> best_patches(n_patches);

                for (int i = 0; i < n_patches; i++) {
                    best_patches[i] = std::move(patches[patch_indices[i]]);
                }

                patches = std::move(best_patches);
            }

            template<class T>
            bool is_homogenous_area(std::vector<ImagePatch<T>> &patches, float noise_std) {

                float std2 = std::accumulate(patches.begin(), patches.end(), 0.0f,
                                             [](auto cur, auto patch) {
                                                 float std = arma::stddev(patch.patch);
                                                 return cur + std * std;
                                             }
                ) * patches.size() / float(patches.size() - 1);

                return std2 < noise_std * noise_std * 1.1;
            }


            template<class T>
            void denoise_patches(std::vector<ImagePatch<T>> &patches, float noise_std) {
                auto mean_patch = get_mean_patch(patches);

                if (is_homogenous_area(patches, noise_std)) {
                    auto mean_value = arma::mean(mean_patch);
                    for (auto &patch: patches) patch.patch.fill(mean_value);
                    return;
                }

                auto covariance_matrix = get_covariance_matrix(patches, mean_patch);
                arma::Mat<T> noise_covariance = covariance_matrix + noise_std * noise_std * arma::eye<arma::Mat<T>>(
                        arma::size(covariance_matrix));
                auto inv_cov = arma::Mat<T>(arma::size(covariance_matrix));
                if (inv(inv_cov, noise_covariance)) {
                    for (auto &patch : patches) {
                        patch.patch = mean_patch + inv_cov * covariance_matrix * (patch.patch - mean_patch);
                    }
                }

            }

            template<class T>
            hoNDArray<T> non_local_bayes_single_image(const hoNDArray<T> &image, float noise_std, int search_window) {


                if (image.get_number_of_dimensions() != 2)
                    throw std::invalid_argument("non_local_bayes: image must be 2 dimensional");

                constexpr int patch_size = 5;
                constexpr int n_patches = 50;

                hoNDArray<T> result(image.dimensions());
                result.fill(0);

                hoNDArray<bool> mask(image.dimensions());
                mask.fill(true);

                hoNDArray<int> count(image.dimensions());
                count.fill(0);

                const vector_td<int, 2> image_dims = vector_td<int, 2>(
                        from_std_vector<size_t, 2>(image.dimensions())
                );

#pragma omp parallel for num_threads(4)
                for (int ky = 0; ky < image.get_size(1); ky++) {
                    for (int kx = 0; kx < image.get_size(0); kx++) {

                        if (mask(kx, ky)) {
                            auto reference_patch = get_patch(image, kx, ky, patch_size, image_dims);
                            auto patches = create_patches(image, kx, ky, patch_size, search_window, image_dims);

                            filter_patches(patches, n_patches, reference_patch);
                            denoise_patches(patches, noise_std);

                            for (auto &patch : patches) {
                                #pragma omp critical
                                add_patch(patch, result, count, patch_size, image_dims);
                                mask(patch.center_x, patch.center_y) = false;
                            }
                        }
                    }
                }

                for (size_t i = 0; i < result.get_number_of_elements(); i++) {
                    result[i] /= count[i];
                }

                return result;

            }


            template<class T>
            hoNDArray<T> non_local_bayes_T(const hoNDArray<T> &image, float noise_std, unsigned int search_window) {

                size_t n_images = image.get_number_of_elements() / (image.get_size(0) * image.get_size(1));

                std::vector<size_t> image_dims = {image.get_size(0), image.get_size(1)};
                size_t image_elements = image_dims[0] * image_dims[1];

                auto result = hoNDArray<T>(image.dimensions());

                hoNDArray<T> result_view;
                result_view.create(image_dims);
                #pragma omp parallel for
                for (int i = 0; i < n_images; i++) {

                    auto image_view = hoNDArray<T>(image_dims, const_cast<T*>(image.get_data_ptr() + i * image_elements));
                    result_view = non_local_bayes_single_image(image_view, noise_std, search_window);

                    memcpy(result.begin() + i * image_elements, result_view.begin(), result_view.get_number_of_bytes());
                }
                return result;
            }
        }


        hoNDArray<float> non_local_bayes(const hoNDArray<float> &image, float noise_std, unsigned int search_window) {
            return non_local_bayes_T(image, noise_std, search_window);
        }

        hoNDArray<std::complex<float>>
        non_local_bayes(const hoNDArray<std::complex<float>> &image, float noise_std, unsigned int search_window) {
            return non_local_bayes_T(image, noise_std, search_window);
        }
    }

}
