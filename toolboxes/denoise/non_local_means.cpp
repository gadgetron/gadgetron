//
// Created by dchansen on 6/19/18.
//

#include <vector_td_utilities.h>
#include "non_local_means.h"


namespace Gadgetron {
    namespace Denoise {

        namespace {

            template<int D> vector_td<float,D*D> get_window(const hoNDArray<float>& image, int x, int y, const vector_td<int,2>& image_dims){

                constexpr int N = D*D;
                vector_td<float,N > window;
                for (int ky = 0; ky < D; ky++){
                    for (int kx = 0; kx < D; kx++){
                        window[kx+ky*D] = image(((kx-D/2)+x+image_dims[0])%image_dims[0],((ky-D/2)+y+image_dims[1])%image_dims[1]);

                    }
                }
                return window;
            };
        }

        hoNDArray<float> non_local_means(const hoNDArray<float>& image, float noise_std, unsigned int search_radius){

            constexpr int D = 5;

            hoNDArray<float> result(image.get_dimensions());
            const float noise_std2 = noise_std*noise_std;
            const vector_td<int,2> image_dims =  vector_td<int,2>(from_std_vector<size_t ,2>(*image.get_dimensions()));


//#pragma omp parallel for collapse(2)
            for (int ky = 0; ky < image.get_size(1); ky++){
                for (int kx = 0; kx < image.get_size(0); kx++){
                    float sum_weight = 0;
                    float sum_value = 0;
                    auto window = get_window<D>(image,kx,ky,image_dims);

                    for (int dy = -search_radius; dy < search_radius; dy++){
                        for (int dx = -search_radius; dx < search_radius; dx++){

                            auto window2 = get_window<D>(image,kx+dx,ky+dy,image_dims);
                            auto diff = window-window2;
                            auto weight = std::exp(-norm_squared(diff)/noise_std2);

                            sum_weight += weight;
                            sum_value += weight*window2[D/2+D*(D/2)];

                        }
                    }

                    result(kx,ky) = sum_value/sum_weight;

                }
            }

            return result;

        }






    }
}
