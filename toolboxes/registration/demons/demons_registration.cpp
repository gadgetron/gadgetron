//
// Created by dchansen on 1/29/20.
//

#include "demons_registration.h"
#include "vector_td_utilities.h"
#include <numeric>
#include "hoNDArray_elemwise.h"
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

        return image[x1 + y1stride + z1stride] * z1weight * x1weight * y1weight
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
        const vector_td<size_t, 3> dims{ output.dimensions()[0], output.dimensions()[1], output.dimensions()[2] };
        for (size_t z = 0; z < dims[2]; z++) {
            for (size_t y = 0; y < dims[1]; y++) {
                size_t offset = y * dims[0];
                for (size_t x = 0; x < dims[0]; x++) {
                    const auto& deformation = deformation_field[x + offset];
                    output[x + offset]
                        = interpolate_point(image, deformation[0] + x, deformation[1] + y, deformation[2] + z, dims);
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

template <class T,unsigned int D, class R>
Gadgetron::hoNDArray<T> Gadgetron::Registration::deform_image(
    const hoNDArray<T>& image, const hoNDArray<vector_td<R, D>>& deformation_field) {

    assert(image.dimensions().size() == D);
    assert(deformation_field.dimensions() == image.dimensions());

    auto output = hoNDArray<T>(image.dimensions());
    interpolation_loop(output, image, deformation_field);
    return output;
}

template hoNDArray<float> Gadgetron::Registration::deform_image(
    const hoNDArray<float>& image, const hoNDArray<vector_td<float, 2>>& deformation_field);
template hoNDArray<float> Gadgetron::Registration::deform_image(
    const hoNDArray<float>& image, const hoNDArray<vector_td<float, 3>>& deformation_field);
template hoNDArray<vector_td<float, 2>> Gadgetron::Registration::deform_image(
    const hoNDArray<vector_td<float, 2>>& image, const hoNDArray<vector_td<float, 2>>& deformation_field);
template hoNDArray<vector_td<float, 3>> Gadgetron::Registration::deform_image(
    const hoNDArray<vector_td<float, 3>>& image, const hoNDArray<vector_td<float, 3>>& deformation_field);
template hoNDArray<std::complex<float>> Gadgetron::Registration::deform_image(
    const hoNDArray<std::complex<float>>& image, const hoNDArray<vector_td<float, 2>>& deformation_field);
template hoNDArray<std::complex<float>> Gadgetron::Registration::deform_image(
    const hoNDArray<std::complex<float>>& image, const hoNDArray<vector_td<float, 3>>& deformation_field);

namespace {

    template <class T>
    vector_td<T, 3> demons_gradient(
        const hoNDArray<T>& image, const vector_td<size_t, 3>& index, const vector_td<size_t, 3>& dims) {
        using namespace std;
        const auto x = index[0];
        const auto y = index[1];
        const auto z = index[2];

        auto x1 = min(max(int(x) - 1, 0), int(dims[0] - 1));
        auto x2 = min(max(int(x) + 1, 0), int(dims[0] - 1));
        auto y1 = min(max(int(y) - 1, 0), int(dims[1] - 1));
        auto y2 = min(max(int(y) + 1, 0), int(dims[1] - 1));
        auto z1 = min(max(int(y) - 1, 0), int(dims[2] - 1));
        auto z2 = min(max(int(y) + 1, 0), int(dims[2] - 1));

        auto zstride = dims[0] * dims[1];
        return vector_td<T,3>{ (image[x2 + y * dims[0] + z * zstride] - image[x1 + y * dims[1] + z * zstride]) / 2,
            (image[x + y2 * dims[0] + z * zstride] - image[x + y1 * dims[0] + z * zstride]) / 2,
            (image[x + y * dims[0] + z2 * zstride] - image[x + y * dims[0] + z1 * zstride]) / 2 };
    }

    template <class T>
    vector_td<T, 2> demons_gradient(
        const hoNDArray<T>& image, const vector_td<size_t, 2>& index, const vector_td<size_t, 2>& dims) {
        using namespace std;
        const auto x = index[0];
        const auto y = index[1];
        auto x1      = min(max(int(x) - 1, 0), int(dims[0] - 1));
        auto x2      = min(max(int(x) + 1, 0), int(dims[0] - 1));
        auto y1      = min(max(int(y) - 1, 0), int(dims[1] - 1));
        auto y2      = min(max(int(y) + 1, 0), int(dims[1] - 1));

        return vector_td<T,2>{ (image[x2 + y * dims[0]] - image[x1 + y * dims[1]]) / 2,
            (image[x + y2 * dims[0]] - image[x + y1 * dims[0]]) / 2 };
    }

    template <class T, unsigned int D>
    vector_td<T, D> demons_point(const hoNDArray<T>& fixed, const hoNDArray<T>& moving, T alpha,
        T beta, const vector_td<size_t, D>& index, const vector_td<size_t, D>& dims) {

        auto fixed_grad   = demons_gradient(fixed, index, dims);
        auto moving_grad  = demons_gradient(moving, index, dims);
        auto average_grad = (fixed_grad + moving_grad) / 2;

        auto it = moving[co_to_idx(index, dims)] - fixed[co_to_idx(index, dims)];
//
        auto result = it * average_grad / (norm_squared(average_grad) + (alpha * it) * (alpha * it) + beta);
        return result;
//        return average_grad*it;
    }

    template <class T, unsigned int D>
    std::enable_if_t<D == 3, hoNDArray<vector_td<T, D>>> demons_step(
        const hoNDArray<T>& fixed, const hoNDArray<T>& moving, T alpha, T beta) {
        assert(fixed.get_number_of_dimensions() == D);
        assert(fixed.dimensions() == moving.dimensions());

        auto dims   = from_std_vector<size_t, D>(fixed.dimensions());
        auto result = hoNDArray<vector_td<T, D>>(fixed.dimensions());

        auto zstride = dims[0] * dims[1];
#pragma omp parallel for if (dims[2] > 1)
        for (size_t z = 0; z < dims[2]; z++) {
            vector_td<size_t, D> index;
            index[2] = z;
            for (size_t y = 0; y < dims[1]; y++) {
                index[1] = y;
                for (size_t x = 0; x < dims[0]; x++) {
                    index[0]                              = x;
                    result[x + y * dims[0] + z * zstride] = demons_point(fixed, moving, alpha, beta, index, dims);
                }
            }
        }
        return result;
    }

    template <class T, unsigned int D>
    std::enable_if_t<D == 2, hoNDArray<vector_td<T, D>>> demons_step(
        const hoNDArray<T>& fixed, const hoNDArray<T>& moving, T alpha, T beta) {
        assert(fixed.get_number_of_dimensions() == D);
        assert(fixed.dimensions() == moving.dimensions());

        auto dims   = from_std_vector<size_t, D>(fixed.dimensions());
        auto result = hoNDArray<vector_td<T, D>>(fixed.dimensions());
        vector_td<size_t, D> index;
        for (size_t y = 0; y < dims[1]; y++) {
            index[1] = y;
            for (size_t x = 0; x < dims[0]; x++) {
                index[0]                = x;
                result[x + y * dims[0]] = demons_point(fixed, moving, alpha, beta, index, dims);
            }
        }

        return result;
    }
}

namespace {

    template <class T, class R, int DIM>
    hoNDArray<T> convolve2D(const hoNDArray<T>& input, const std::vector<R>& kernel) {
        using namespace std;

        hoNDArray<T> output(input.dimensions());
        const auto dims    = from_std_vector<size_t, 2>(input.dimensions());
        const auto strides = vector_td<size_t, 2>(1, dims[0]);

        const int kernel_size = kernel.size();
        constexpr int dim     = DIM;
        vector_td<size_t, 2> index;

        for (size_t y = 0; y < dims[1]; y++) {
            index[1] = y;
            for (size_t x = 0; x < dims[0]; x++) {
                index[0]    = x;
                T summation = T(0);
                for (int k = 0; k < kernel_size; k++) {
                    auto offset = (max(min(k - kernel_size / 2 + int(index[dim]), int(dims[dim] - 1)), 0) - index[dim])
                                  * ((long long)(strides[dim]));
                    summation += kernel[k] * input[x + y * dims[0] + offset];
                }
                output[x + y * dims[0]] = summation;
            }
        }

        return output;
    }
    template <class T, class R, int DIM>
    hoNDArray<T> convolve3D(const hoNDArray<T>& input, const std::vector<R>& kernel) {
        using namespace std;

        hoNDArray<T> output(input.dimensions());
        const auto dims    = from_std_vector<size_t, 3>(input.dimensions());
        const auto strides = vector_td<size_t, 3>(1, dims[0], dims[0] * dims[2]);

        const int kernel_size = kernel.size();
        constexpr int dim     = DIM;
        vector_td<size_t, 3> index;

        for (size_t z = 0; z < dims[2]; z++) {
            index[0] = z;
            for (size_t y = 0; y < dims[1]; y++) {
                index[1] = y;
                for (size_t x = 0; x < dims[0]; x++) {
                    index[0]    = x;
                    T summation = T(0);
                    for (int k = 0; k < kernel_size; k++) {
                        auto offset
                            = (max(min(k - kernel_size / 2 + int(index[dim]), int(dims[dim] - 1)), 0) - index[dim])
                              * ((long long)(strides[dim]));
                        summation += kernel[k] * input[x + y * dims[0] + offset];
                    }
                    output[x + y * dims[0]] = summation;
                }
            }
        }
        return output;
    }

    std::vector<float> calculate_gauss_kernel(float sigma) {
        auto kernel = std::vector<float>(std::max(int(sigma * 3) * 2 + 1, 1));

        for (long long k = 0; k < (long long)kernel.size(); k++) {
            float x   = float(k) - float(kernel.size() / 2);
            kernel[k] = std::exp(-0.5f * (x / sigma) * (x / sigma));
        }

        auto sum = std::accumulate(kernel.begin(), kernel.end(), 0.0f);
        std::transform(kernel.begin(), kernel.end(), kernel.begin(), [sum](auto& val) { return val / sum; });
        return kernel;
    }

}
template <class T> hoNDArray<T> gaussian_filter(const hoNDArray<T>& image, float sigma) {

    auto kernel = calculate_gauss_kernel(sigma);

    switch (image.get_number_of_dimensions()) {
    case 2: {
        auto result = convolve2D<T, float, 0>(image, kernel);
        result      = convolve2D<T, float, 1>(result, kernel);
        return result;
    }
    case 3: {
        auto result = convolve3D<T, float, 0>(image, kernel);
        result      = convolve3D<T, float, 1>(result, kernel);
        result      = convolve3D<T, float, 2>(result, kernel);
        return result;
    }
    default: throw std::runtime_error("Gaussian filter only support 2 and 3 D images");
    }
}
template hoNDArray<float> gaussian_filter(const hoNDArray<float>& image, float sigma);
template hoNDArray<double> gaussian_filter(const hoNDArray<double>& image, float sigma);
template hoNDArray<vector_td<float, 2>> gaussian_filter(const hoNDArray<vector_td<float, 2>>& image, float sigma);
template hoNDArray<vector_td<float, 3>> gaussian_filter(const hoNDArray<vector_td<float, 3>>& image, float sigma);

namespace {
    using namespace Gadgetron::Registration;
    template <class T, unsigned int D>
    hoNDArray<vector_td<T, D>> compose_fields(
        const hoNDArray<vector_td<T, D>>& update_field, const hoNDArray<vector_td<T, D>>& vfield) {
        auto resulting_field = deform_image(vfield, update_field);
        resulting_field += update_field;
        return resulting_field;
    }

    template <class T, unsigned int D>
    hoNDArray<vector_td<T, D>> vector_field_exponential(const hoNDArray<vector_td<T,D>>& vector_field) {

        auto maximum_vector_length = std::accumulate(vector_field.begin(), vector_field.end(), T(0),
            [](auto current, const auto& next) { return std::max(current, norm(next)); });
        int n_iteration = std::ceil(2.0 + 0.5f * std::log2(maximum_vector_length));
        auto field_exponent = hoNDArray<vector_td<T, D>>(vector_field.dimensions());
        std::transform(vector_field.begin(), vector_field.end(), field_exponent.begin(),
            [n_iteration](const auto& val) { return val * std::pow(T(2), T(-n_iteration)); });

        for (size_t i = 0; i < n_iteration; i++) {
            field_exponent = compose_fields(field_exponent, field_exponent);
        }
        return field_exponent;
    }
}


template<class T, unsigned int D>
hoNDArray<vector_td<T,D>> Gadgetron::Registration::diffeomorphic_demons(const hoNDArray<T>& fixed, const hoNDArray<T>& moving, unsigned int iterations,float sigma){
    auto vector_field = demons_step<T,D>(fixed,moving,2.0f,1e-6f);
    vector_field = gaussian_filter(vector_field,sigma);
    vector_field = vector_field_exponential(vector_field);
    for (size_t i = 1; i < iterations; i++){
        auto current_fixed = deform_image(fixed,vector_field);
        auto update_field = demons_step<T,D>(current_fixed,moving,2.0f,1e-6f);
        update_field = vector_field_exponential(update_field);
        vector_field = compose_fields(update_field,vector_field);
        vector_field = gaussian_filter(vector_field,sigma);
    }

    return vector_field;
}

template hoNDArray<vector_td<float,2>> Gadgetron::Registration::diffeomorphic_demons(const hoNDArray<float>& , const hoNDArray<float>& ,unsigned int, float );

