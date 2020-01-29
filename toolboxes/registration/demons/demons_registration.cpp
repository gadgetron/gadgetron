//
// Created by dchansen on 1/29/20.
//

#include "demons_registration.h"
using namespace Gadgetron;

namespace {

    //    template<class T, class R, unsigned int D, class...INDEX>
    //    std::enable_if_t< sizeof...(INDEX) < D> interpolation_loop(hoNDArray<T>& output,const hoNDArray<T>& image,
    //    const hoNDArray<vector_td<R,D>>& deformation_field, const vector_td<size_t,D> image_dims, INDEX... index){
    //        constexpr size_t COUNTER= sizeof...(INDEX);
    //        for (size_t i =0; i < image_dims[COUNTER]; i++){
    //            interpolation_loop(output,image,deformation_field,image_dims,i,index...);
    //        }
    //    }
    //
    //    template<class T, class R, unsigned int D, class...INDEX>
    //    std::enable_if_t< sizeof...(INDEX) == D> interpolation_loop(hoNDArray<T>& output,const hoNDArray<T>& image,
    //    const hoNDArray<vector_td<R,D>>& deformation_field, const vector_td<size_t,D> image_dims, INDEX... index){
    //
    //        auto index_vector = vector_td<size_t,D>{index...};
    //        auto idx = dot(index_vector,image_dims);
    //        const auto coords = deformation_field[idx]+index_vector;
    //        auto deformed_index =
    //
    //
    //
    //    }

    template <class T, class R>
    T interpolate_point(const hoNDArray<T>& image, R x, R y, R z, const vector_td<size_t, 3>& dims) {
        using namespace std;
        auto x1       = min(max(int(x), 0), int(dims[0] - 1));
        auto x2       = min(max(int(x) + 1, 0), int(dims[0] - 1));
        auto x2weight = x - x1;
        auto x1weight = 1 - x2weight;

        auto y1       = min(max(int(y), 0), int(dims[1] - 1));
        auto y2       = min(max(int(y) + 1, 0), int(dims[1] - 1));
        auto y2weight = y - y1;
        auto y1weight = 1 - y2weight;

        auto y1stride = y1 * dims[0];
        auto y2stride = y2 * dims[0];

        auto z1       = min(max(int(z), 0), int(dims[2] - 1));
        auto z2       = min(max(int(z) + 1, 0), int(dims[2] - 1));
        auto z2weight = z - z1;
        auto z1weight = 1 - z2weight;

        auto z1stride = z1 * dims[1] * dims[0];
        auto z2stride = z2 * dims[1] * dims[0];

        return   image[x1 + y1stride + z1stride] * z1weight * x1weight * y1weight
               + image[x2 + y1stride + z1stride] * z1weight * x2weight * y1weight
               + image[x1 + y2stride + z1stride] * z1weight * x1weight * y2weight
               + image[x2 + y2stride + z1stride] * z1weight * x2weight * y2weight

               + image[x1 + y1stride + z2stride] * z2weight * x1weight * y1weight
               + image[x2 + y1stride + z2stride] * z2weight * x2weight * y1weight
               + image[x1 + y2stride + z2stride] * z2weight * x1weight * y2weight
               + image[x2 + y2stride + z2stride] * z2weight * x2weight * y2weight;
    }
    template <class T, class R>
    void interpolation_loop(
        hoNDArray<T>& output, const hoNDArray<T>& image, const hoNDArray<vector_td<R, 3>>& deformation_field) {
        const vector_td<size_t, 3> dims{ output.dimensions()[0], output.dimensions()[1] ,output.dimensions()[2]};
        for (size_t z = 0; z < dims[2]; z++) {
            for (size_t y = 0; y < dims[1]; y++) {
                size_t offset = y * dims[0];
                for (size_t x = 0; x < dims[0]; x++) {
                    const auto& deformation = deformation_field[x + offset];
                    output[x + offset]      = interpolate_point(image, deformation[0] + x, deformation[1] + y, deformation[2]+z, dims);
                }
            }
        }
    }

    template <class T, class R>
    T interpolate_point(const hoNDArray<T>& image, R x, R y, const vector_td<size_t, 2>& dims) {
        using namespace std;
        auto x1       = min(max(int(x), 0), int(dims[0] - 1));
        auto x2       = min(max(int(x) + 1, 0), int(dims[0] - 1));
        auto x2weight = x - x1;
        auto x1weight = 1 - x2weight;

        auto y1       = min(max(int(y), 0), int(dims[1] - 1));
        auto y2       = min(max(int(y) + 1, 0), int(dims[1] - 1));
        auto y2weight = y - y1;
        auto y1weight = 1 - y2weight;

        auto y1stride = y1 * dims[0];
        auto y2stride = y2 * dims[0];

        return image[x1 + y1stride] * x1weight * y1weight + image[x2 + y1stride] * x2weight * y1weight
               + image[x1 + y2stride] * x1weight * y2weight + image[x2 + y2stride] * x2weight * y2weight;
    }

    template <class T, class R>
    void interpolation_loop(
        hoNDArray<T>& output, const hoNDArray<T>& image, const hoNDArray<vector_td<R, 2>>& deformation_field) {
        const vector_td<size_t, 2> dims{ output.dimensions()[0], output.dimensions()[1] };
        for (size_t y = 0; y < dims[1]; y++) {
            size_t offset = y * dims[0];
            for (size_t x = 0; x < dims[0]; x++) {
                const auto& deformation = deformation_field[x + offset];
                output[x + offset]      = interpolate_point(image, deformation[0] + x, deformation[1] + y, dims);
            }
        }
    }


}

template <class T, class R, unsigned int D>
Gadgetron::hoNDArray<T> Gadgetron::Registration::deform_image(
    const hoNDArray<T>& image, const hoNDArray<vector_td<R, D>>& deformation_field) {

    assert(image.dimensions().size() == D);
    assert(deformation_field.dimensions() == image.dimensions());

    auto output = hoNDArray<T>(image.dimensions());
    interpolation_loop(output, image, deformation_field);
    return output;
}



template hoNDArray<float> Gadgetron::Registration::deform_image(const hoNDArray<float> &image, const hoNDArray<vector_td<float, 2> > &deformation_field);
template hoNDArray<float> Gadgetron::Registration::deform_image(const hoNDArray<float> &image, const hoNDArray<vector_td<float, 3> > &deformation_field);
template hoNDArray<vector_td<float,2>> Gadgetron::Registration::deform_image(const hoNDArray<vector_td<float,2>> &image, const hoNDArray<vector_td<float, 2> > &deformation_field);
template hoNDArray<vector_td<float,3>> Gadgetron::Registration::deform_image(const hoNDArray<vector_td<float,3>> &image, const hoNDArray<vector_td<float, 3> > &deformation_field);
template hoNDArray<std::complex<float>> Gadgetron::Registration::deform_image(const hoNDArray<std::complex<float>> &image, const hoNDArray<vector_td<float, 2> > &deformation_field);
template hoNDArray<std::complex<float>> Gadgetron::Registration::deform_image(const hoNDArray<std::complex<float>> &image, const hoNDArray<vector_td<float, 3> > &deformation_field);

