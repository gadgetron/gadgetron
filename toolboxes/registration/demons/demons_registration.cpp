//
// Created by dchansen on 1/29/20.
//

#include "demons_registration.h"
#include "hoNDArray_elemwise.h"
#include "hoNDArray_utils.h"
#include "vector_td_utilities.h"
#include <numeric>

#include "hoNDInterpolator.h"
#include "hoNDArray_fileio.h"
using namespace Gadgetron;

namespace {

template <class T, class R>
inline void interpolation_loop(hoNDArray<T>& output,
                               const hoNDArray<vector_td<R, 2>>& deformation_field,
                               hoNDInterpolatorBSpline<hoNDArray<T>, 2>& interpolator) {
    const vector_td<size_t, 2> dims{output.dimensions()[0], output.dimensions()[1]};
    for (size_t y = 0; y < dims[1]; y++) {
        size_t offset = y * dims[0];
        for (size_t x = 0; x < dims[0]; x++) {
            const auto& deformation = deformation_field[x + offset];
            output[x + offset] = interpolator(deformation[0] + x, deformation[1] + y);
        }
    }
}

template <class T, class R>
void interpolation_loop(hoNDArray<T>& output, const hoNDArray<vector_td<R, 3>>& deformation_field,
                        hoNDInterpolatorBSpline<hoNDArray<T>, 3>& interpolator) {
    const vector_td<size_t, 3> dims{output.dimensions()[0], output.dimensions()[1],
                                    output.dimensions()[2]};
    for (size_t z = 0; z < dims[2]; z++) {
        for (size_t y = 0; y < dims[1]; y++) {
            size_t offset = y * dims[0];
            for (size_t x = 0; x < dims[0]; x++) {
                const auto& deformation = deformation_field[x + offset];
                output[x + offset] =
                    interpolator(deformation[0] + x, deformation[1] + y, deformation[2] + z);
            }
        }
    }
}

template <class T, class R>
T interpolate_point(const hoNDArray<T>& image, R x, R y, R z, const vector_td<size_t, 3>& dims) {
    using namespace std;
    auto x1 = min(max(int(x), 0), int(dims[0] - 1));
    auto x2 = min(max(int(x) + 1, 0), int(dims[0] - 1));
    auto x2weight = x - x1;
    auto x1weight = 1 - x2weight;

    auto y1 = min(max(int(y), 0), int(dims[1] - 1));
    auto y2 = min(max(int(y) + 1, 0), int(dims[1] - 1));
    auto y2weight = y - y1;
    auto y1weight = 1 - y2weight;

    auto y1stride = y1 * dims[0];
    auto y2stride = y2 * dims[0];

    auto z1 = min(max(int(z), 0), int(dims[2] - 1));
    auto z2 = min(max(int(z) + 1, 0), int(dims[2] - 1));
    auto z2weight = z - z1;
    auto z1weight = 1 - z2weight;

    auto z1stride = z1 * dims[1] * dims[0];
    auto z2stride = z2 * dims[1] * dims[0];

    return image[x1 + y1stride + z1stride] * z1weight * x1weight * y1weight +
           image[x2 + y1stride + z1stride] * z1weight * x2weight * y1weight +
           image[x1 + y2stride + z1stride] * z1weight * x1weight * y2weight +
           image[x2 + y2stride + z1stride] * z1weight * x2weight * y2weight

           + image[x1 + y1stride + z2stride] * z2weight * x1weight * y1weight +
           image[x2 + y1stride + z2stride] * z2weight * x2weight * y1weight +
           image[x1 + y2stride + z2stride] * z2weight * x1weight * y2weight +
           image[x2 + y2stride + z2stride] * z2weight * x2weight * y2weight;
}
template <class T, class R>
void interpolation_loop(hoNDArray<T>& output, const hoNDArray<T>& image,
                        const hoNDArray<vector_td<R, 3>>& deformation_field) {
    const vector_td<size_t, 3> dims{output.dimensions()[0], output.dimensions()[1],
                                    output.dimensions()[2]};
    for (size_t z = 0; z < dims[2]; z++) {
        for (size_t y = 0; y < dims[1]; y++) {
            size_t offset = y * dims[0];
            for (size_t x = 0; x < dims[0]; x++) {
                const auto& deformation = deformation_field[x + offset];
                output[x + offset] = interpolate_point(
                    image, deformation[0] + x, deformation[1] + y, deformation[2] + z, dims);
            }
        }
    }
}

template <class T, class R>
T interpolate_point(const hoNDArray<T>& image, R x, R y, const vector_td<size_t, 2>& dims) {
    using namespace std;
    auto x1 = min(max(int(x), 0), int(dims[0] - 1));
    auto x2 = min(max(int(x) + 1, 0), int(dims[0] - 1));
    auto x2weight = x - x1;
    auto x1weight = 1 - x2weight;

    auto y1 = min(max(int(y), 0), int(dims[1] - 1));
    auto y2 = min(max(int(y) + 1, 0), int(dims[1] - 1));
    auto y2weight = y - y1;
    auto y1weight = 1 - y2weight;

    auto y1stride = y1 * dims[0];
    auto y2stride = y2 * dims[0];

    return image[x1 + y1stride] * x1weight * y1weight + image[x2 + y1stride] * x2weight * y1weight +
           image[x1 + y2stride] * x1weight * y2weight + image[x2 + y2stride] * x2weight * y2weight;
}

template <class T, class R>
void interpolation_loop(hoNDArray<T>& output, const hoNDArray<T>& image,
                        const hoNDArray<vector_td<R, 2>>& deformation_field) {
    const vector_td<size_t, 2> dims{output.dimensions()[0], output.dimensions()[1]};
    for (size_t y = 0; y < dims[1]; y++) {
        size_t offset = y * dims[0];
        for (size_t x = 0; x < dims[0]; x++) {
            const auto& deformation = deformation_field[x + offset];
            output[x + offset] =
                interpolate_point(image, deformation[0] + x, deformation[1] + y, dims);
        }
    }
}

} // namespace

template <class T, unsigned int D, class R>
Gadgetron::hoNDArray<T>
Gadgetron::Registration::deform_image(const hoNDArray<T>& image,
                                      const hoNDArray<vector_td<R, D>>& deformation_field) {

    if (image.dimensions().size() != D)
        throw std::runtime_error("Image dimensions do not match");
    if (deformation_field.dimensions() != image.dimensions())
        throw std::runtime_error("Deformation field and image dimensions do not match");

    auto output = hoNDArray<T>(image.dimensions());
    interpolation_loop(output, image, deformation_field);
    return output;
}

template hoNDArray<float>
Gadgetron::Registration::deform_image(const hoNDArray<float>& image,
                                      const hoNDArray<vector_td<float, 2>>& deformation_field);
template hoNDArray<float>
Gadgetron::Registration::deform_image(const hoNDArray<float>& image,
                                      const hoNDArray<vector_td<float, 3>>& deformation_field);
template hoNDArray<vector_td<float, 2>>
Gadgetron::Registration::deform_image(const hoNDArray<vector_td<float, 2>>& image,
                                      const hoNDArray<vector_td<float, 2>>& deformation_field);
template hoNDArray<vector_td<float, 3>>
Gadgetron::Registration::deform_image(const hoNDArray<vector_td<float, 3>>& image,
                                      const hoNDArray<vector_td<float, 3>>& deformation_field);
template hoNDArray<std::complex<float>>
Gadgetron::Registration::deform_image(const hoNDArray<std::complex<float>>& image,
                                      const hoNDArray<vector_td<float, 2>>& deformation_field);
template hoNDArray<std::complex<float>>
Gadgetron::Registration::deform_image(const hoNDArray<std::complex<float>>& image,
                                      const hoNDArray<vector_td<float, 3>>& deformation_field);

template <class T, unsigned int D, class R>
hoNDArray<T>
Gadgetron::Registration::deform_image_bspline(const hoNDArray<T>& image,
                                              const hoNDArray<vector_td<R, D>>& deformation_field) {

    assert(image.dimensions().size() == D);
    assert(deformation_field.dimensions() == image.dimensions());

    hoNDBoundaryHandlerBorderValue<hoNDArray<T>> boundaryhandler(image);
    hoNDInterpolatorBSpline<hoNDArray<T>, D> interpolator(image, boundaryhandler);

    auto output = hoNDArray<T>(image.dimensions());
    interpolation_loop(output, deformation_field, interpolator);
    return output;
}

template hoNDArray<float> Gadgetron::Registration::deform_image_bspline(
    const hoNDArray<float>& image, const hoNDArray<vector_td<float, 2>>& deformation_field);
template hoNDArray<float> Gadgetron::Registration::deform_image_bspline(
    const hoNDArray<float>& image, const hoNDArray<vector_td<float, 3>>& deformation_field);
template hoNDArray<std::complex<float>> Gadgetron::Registration::deform_image_bspline(
    const hoNDArray<std::complex<float>>& image,
    const hoNDArray<vector_td<float, 2>>& deformation_field);
template hoNDArray<std::complex<float>> Gadgetron::Registration::deform_image_bspline(
    const hoNDArray<std::complex<float>>& image,
    const hoNDArray<vector_td<float, 3>>& deformation_field);
namespace {

template <class T>
vector_td<T, 3> demons_gradient(const hoNDArray<T>& image, const vector_td<size_t, 3>& index,
                                const vector_td<size_t, 3>& dims) {
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
    return vector_td<T, 3>{
        (image[x2 + y * dims[0] + z * zstride] - image[x1 + y * dims[1] + z * zstride]) / 2,
        (image[x + y2 * dims[0] + z * zstride] - image[x + y1 * dims[0] + z * zstride]) / 2,
        (image[x + y * dims[0] + z2 * zstride] - image[x + y * dims[0] + z1 * zstride]) / 2};
}

template <class T>
vector_td<T, 2> demons_gradient(const hoNDArray<T>& image, const vector_td<size_t, 2>& index,
                                const vector_td<size_t, 2>& dims) {
    using namespace std;
    const auto x = index[0];
    const auto y = index[1];
    auto x1 = min(max(int(x) - 1, 0), int(dims[0] - 1));
    auto x2 = min(max(int(x) + 1, 0), int(dims[0] - 1));
    auto y1 = min(max(int(y) - 1, 0), int(dims[1] - 1));
    auto y2 = min(max(int(y) + 1, 0), int(dims[1] - 1));
    auto val1 = image[x2 + y * dims[0]];
    auto val2 =  image[x1 + y * dims[0]];
    auto val3 = image[x + y2 * dims[0]];
    auto val4 = image[x + y1 * dims[0]];
    auto res=  vector_td<T, 2>{(val1-val2) / 2,
                           (val3-val4) / 2};
    return res;
}

template <class T, unsigned int D>
vector_td<T, D> demons_point(const hoNDArray<T>& fixed, const hoNDArray<T>& moving, T alpha, T beta,
                             const vector_td<size_t, D>& index, const vector_td<size_t, D>& dims,
                             T noise_sigma) {

    auto fixed_grad = demons_gradient(fixed, index, dims);
    auto moving_grad = demons_gradient(moving, index, dims);
    auto average_grad = (fixed_grad + moving_grad) / 2;

    auto it = moving[co_to_idx(index, dims)] - fixed[co_to_idx(index, dims)];
    //
    auto result =
        it * average_grad / (norm_squared(average_grad) + (alpha * it) * (alpha * it) + beta);
    if (noise_sigma > 0)
        result *= T(1) - std::exp(-(it * it) / (2 * noise_sigma));
    return result;
    //        return average_grad*it;
}

template <class T, unsigned int D> struct DemonStep {};

template <class T> struct DemonStep<T, 3> {

    const hoNDArray<T>& moving;
    T alpha;
    T beta;
    T noise_sigma;

    hoNDArray<vector_td<T, 3>> operator()(const hoNDArray<T>& fixed) {
        assert(fixed.get_number_of_dimensions() == 2);
        assert(fixed.dimensions() == moving.dimensions());

        auto dims = from_std_vector<size_t, 2>(fixed.dimensions());
        auto result = hoNDArray<vector_td<T, 2>>(fixed.dimensions());

        auto zstride = dims[0] * dims[1];
#pragma omp parallel for if (dims[2] > 1)
        for (size_t z = 0; z < dims[2]; z++) {
            vector_td<size_t, 2> index;
            index[2] = z;
            for (size_t y = 0; y < dims[1]; y++) {
                index[1] = y;
                for (size_t x = 0; x < dims[0]; x++) {
                    index[0] = x;
                    result[x + y * dims[0] + z * zstride] =
                        demons_point(fixed, moving, alpha, beta, index, dims);
                }
            }
        }
        return result;
    }
};

template <class T> struct DemonStep<T, 2> {

    const hoNDArray<T>& moving;
    T alpha;
    T beta;
    T noise_sigma;

    hoNDArray<vector_td<T, 2>> operator()(const hoNDArray<T>& fixed) {
        assert(fixed.get_number_of_dimensions() == 2);
        assert(fixed.dimensions() == moving.dimensions());

        auto dims = from_std_vector<size_t, 2>(fixed.dimensions());
        auto result = hoNDArray<vector_td<T, 2>>(fixed.dimensions());
        vector_td<size_t, 2> index;
        for (size_t y = 0; y < dims[1]; y++) {
            index[1] = y;
            for (size_t x = 0; x < dims[0]; x++) {
                index[0] = x;
                result(x, y) = demons_point(fixed, moving, alpha, beta, index, dims, noise_sigma);
            }
        }

        return result;
    }
};

template <class T, unsigned int D> struct NGFDemonsStep {};
template <class T> struct NGFDemonsStep<T, 2> {
    static constexpr unsigned int D = 2;
    NGFDemonsStep(const hoNDArray<T>& moving, T alpha, T beta, T gradient_eps)
        : alpha{alpha}, beta{beta}, gradient_eps{gradient_eps}, moving_grad{create_gradient_array(
                                                                    moving)} {

//        write_nd_array(&moving,"moving.real");
    }

    const std::vector<hoNDArray<T>> moving_grad;
    const T alpha;
    const T beta;
    const T gradient_eps;
    // Implements a Demons step with normalized gradient fields.
    // Uses the estimated channel certainty of
    // "Improving Registration Using Multi-Channel Diffeomorphic Demons Combined
    // withCertainty Maps", Fosberg et. al 2011
    hoNDArray<vector_td<T, 2>> operator()(const hoNDArray<T>& fixed) {
        auto fixed_grad = create_gradient_array(fixed);

//        write_nd_array(&fixed_grad[0],"fixed_grad1.real");
//        write_nd_array(&fixed_grad[1],"fixed_grad2.real");
        assert(fixed.get_number_of_dimensions() == D);

        const auto dims = from_std_vector<size_t, D>(fixed.dimensions());
        auto result = hoNDArray<vector_td<T, D>>(fixed.dimensions());
        vector_td<size_t, 2> index;
        for (size_t y = 0; y < dims[1]; y++) {
            index[1] = y;
            for (size_t x = 0; x < dims[0]; x++) {
                index[0] = x;

                vector_td<vector_td<T, D>,D> steps;
                vector_td<T, D> certainties;
                for (int c = 0; c < D; c++) {
                    using namespace Gadgetron::Indexing;
                    steps[c] =
                        demons_point<T,D>(fixed_grad[c], moving_grad[c],
                                     alpha, beta, index, dims, 0);
                    auto loc_grad = (demons_gradient<T>(fixed_grad[c], index, dims) +
                                     demons_gradient<T>(moving_grad[c],index, dims)) /
                                    2;
                    certainties[c] = outer_norm(loc_grad);
                }

//                vector_td<T,D> step = sum(steps * certainties);
                vector_td<T,D> step = sum(steps)/2;

//                if (sum(certainties) <= 0.0f) step *= 0;
//                else step /= sum(certainties);



                if (std::isnan(step[0]) || std::isnan(step[1])) throw std::runtime_error("BATMAN!");
                if (std::isinf(step[0]) || std::isinf(step[1])) throw std::runtime_error("AND BEYOND!");
                result(x, y) = step;
            }
        }
        return result;
    }

  private:
    T outer_norm(const vector_td<T, D>& v) {
        T result = 0;
        for (int n = 0; n < D; n++) {
            for (int m = 0; m < D; m++) {
                auto vv = v[n] * v[m];
                result += vv * vv;
            }
        }
        result = std::sqrt(result);
//        if (std::isnan(result)) throw std::runtime_error("BATMAN");
        return result;
    }

    auto create_gradient_array(const hoNDArray<T>& image) {

        auto grad_dim = image.dimensions();

        auto dims = from_std_vector<size_t, D>(image.dimensions());
        auto result = std::vector<hoNDArray<T>>(D,hoNDArray<T>(grad_dim));
        vector_td<size_t, D> index;
        for (size_t y = 0; y < dims[1]; y++) {
            index[1] = y;
            for (size_t x = 0; x < dims[0]; x++) {
                index[0] = x;
                auto dg = demons_gradient(image, index, dims);
                dg /= std::sqrt(dot(dg, dg) + gradient_eps);
                for (int c = 0; c < D; c++) {
                    result[c](x, y) = dg[c];
                }
            }
        }
        return result;
    }
};
} // namespace

hoNDArray<vector_td<float, 2>>
Gadgetron::Registration::demons_step_ext(const hoNDArray<float>& fixed,
                                         const hoNDArray<float>& moving, float alpha, float beta,
                                         float noise_sigma) {

    auto stepper = DemonStep<float, 2>{moving, alpha, beta, noise_sigma};
    return stepper(fixed);
}

namespace {

template <class T, class R, int DIM>
hoNDArray<T> convolve2D(const hoNDArray<T>& input, const std::vector<R>& kernel) {
    using namespace std;

    hoNDArray<T> output(input.dimensions());
    const auto dims = from_std_vector<size_t, 2>(input.dimensions());
    const auto strides = vector_td<size_t, 2>(1, dims[0]);

    const int kernel_size = kernel.size();
    constexpr int dim = DIM;
    vector_td<int, 2> index;

    long long dim_len = dims[dim];

    for (int y = 0; y < dims[1]; y++) {
        index[1] = y;
        for (int x = 0; x < dims[0]; x++) {
            index[0] = x;
            T summation = T(0);
            for (int k = 0; k < kernel_size; k++) {
                int kl = k-kernel_size/2+int(index[dim]);
                if (kl >= 0 && kl < dims[dim]) {
//                    long long offset =
//                        (max(min(k - kernel_size / 2 + int(index[dim]), int(dims[dim] - 1)), 0) -
//                         index[dim]) *
//                        ((long long)(strides[dim]));

                    long long offset = (kl-index[dim])*((long long)(strides[dim]));
                    summation += kernel[k] * input[x + y * dims[0] + offset];
                }
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
    const auto dims = from_std_vector<size_t, 3>(input.dimensions());
    const auto strides = vector_td<size_t, 3>(1, dims[0], dims[0] * dims[2]);

    const int kernel_size = kernel.size();
    constexpr int dim = DIM;
    vector_td<size_t, 3> index;

    for (size_t z = 0; z < dims[2]; z++) {
        index[0] = z;
        for (size_t y = 0; y < dims[1]; y++) {
            index[1] = y;
            for (size_t x = 0; x < dims[0]; x++) {
                index[0] = x;
                T summation = T(0);
                for (int k = 0; k < kernel_size; k++) {
                    auto offset =
                        (max(min(k - kernel_size / 2 + int(index[dim]), int(dims[dim] - 1)), 0) -
                         index[dim]) *
                        ((long long)(strides[dim]));
                    summation += kernel[k] * input[x + y * dims[0] + offset];
                }
                output[x + y * dims[0]] = summation;
            }
        }
    }
    return output;
}

std::vector<float> calculate_gauss_kernel(float sigma) {
    int lw = int(sigma * 4 + 0.5f);
    auto kernel = std::vector<float>(std::max(lw * 2 + 1, 1));
    auto sigma2 = sigma * sigma;
    kernel[lw] = 1.0f;

    for (long long k = 1; k < lw + 1; k++) {
        float x = std::exp(-0.5 * float(k * k) / sigma2);
        kernel[lw + k] = x;
        kernel[lw - k] = x;
    }

    auto sum = std::accumulate(kernel.begin(), kernel.end(), 0.0f);
    std::transform(kernel.begin(), kernel.end(), kernel.begin(),
                   [sum](auto& val) { return val / sum; });
    return kernel;
}

} // namespace
template <class T>
hoNDArray<T> Gadgetron::Registration::gaussian_filter(const hoNDArray<T>& image, float sigma) {

    auto kernel = calculate_gauss_kernel(sigma);

    switch (image.get_number_of_dimensions()) {
    case 2: {
        auto result = convolve2D<T, float, 0>(image, kernel);
        result = convolve2D<T, float, 1>(result, kernel);
        return result;
    }
    case 3: {
        auto result = convolve3D<T, float, 0>(image, kernel);
        result = convolve3D<T, float, 1>(result, kernel);
        result = convolve3D<T, float, 2>(result, kernel);
        return result;
    }
    default:
        throw std::runtime_error("Gaussian filter only support 2 and 3 D images");
    }
}
template hoNDArray<float> Gadgetron::Registration::gaussian_filter(const hoNDArray<float>& image,
                                                                   float sigma);
template hoNDArray<double> Gadgetron::Registration::gaussian_filter(const hoNDArray<double>& image,
                                                                    float sigma);
template hoNDArray<vector_td<float, 2>>
Gadgetron::Registration::gaussian_filter(const hoNDArray<vector_td<float, 2>>& image, float sigma);
template hoNDArray<vector_td<float, 3>>
Gadgetron::Registration::gaussian_filter(const hoNDArray<vector_td<float, 3>>& image, float sigma);

template <class T, unsigned int D>
hoNDArray<vector_td<T, D>>
Gadgetron::Registration::compose_fields(const hoNDArray<vector_td<T, D>>& update_field,
                                        const hoNDArray<vector_td<T, D>>& vfield) {
    auto resulting_field = deform_image(vfield, update_field);
    resulting_field += update_field;
    return resulting_field;
}

template hoNDArray<vector_td<float, 2>>
Gadgetron::Registration::compose_fields(const hoNDArray<vector_td<float, 2>>&,
                                        const hoNDArray<vector_td<float, 2>>&);

template <class T, unsigned int D>
hoNDArray<vector_td<T, D>>
Gadgetron::Registration::vector_field_exponential(const hoNDArray<vector_td<T, D>>& vector_field) {

    auto maximum_vector_length = std::accumulate(
        vector_field.begin(), vector_field.end(), T(0),
        [](auto current, const auto& next) { return std::max(current, norm(next)); });
    auto fn_iteration =std::ceil(2.0 + 0.5f * std::log2(maximum_vector_length));
    int n_iteration = fn_iteration;
    auto field_exponent = hoNDArray<vector_td<T, D>>(vector_field.dimensions());
    std::transform(
        vector_field.begin(), vector_field.end(), field_exponent.begin(),
        [n_iteration](const auto& val) { return val * std::pow(T(2), T(-n_iteration)); });

    for (size_t i = 0; i < n_iteration; i++) {
        field_exponent = compose_fields(field_exponent, field_exponent);
    }
    return field_exponent;
}

template hoNDArray<vector_td<float, 2>> Gadgetron::Registration::vector_field_exponential(
    const hoNDArray<vector_td<float, 2>>& vector_field);
namespace {
using namespace Gadgetron::Registration;

template <class T, unsigned int D, class STEPPER>
hoNDArray<vector_td<T, D>>
base_diffeomorphic_demons_impl(const hoNDArray<T>& fixed, hoNDArray<vector_td<T, D>> vector_field,
                               hoNDArray<vector_td<T, D>> predictor_field, STEPPER stepper,
                               unsigned int iterations, float sigma) {

    predictor_field *= 0.5;
    for (size_t i = 0; i < iterations; i++) {
        GDEBUG("Starting iteration %d\n",i);
        auto current_fixed = deform_image(fixed, vector_field);
        auto update_field = stepper(current_fixed);
        update_field = compose_fields(predictor_field, update_field);
        update_field = vector_field_exponential(update_field);
        vector_field = compose_fields(update_field, vector_field);
        vector_field = gaussian_filter(vector_field, sigma);
        predictor_field = std::move(update_field);
        predictor_field *= 0.5;
    }
    return vector_field;
}

template <class T, class STEPPER>
auto base_diffeomorphic_demons(const hoNDArray<T>& fixed, STEPPER& stepper, unsigned int iterations,
                               float sigma) {
    auto vector_field = stepper(fixed);

    vector_field = vector_field_exponential(vector_field);
    vector_field = gaussian_filter(vector_field, sigma);
    return base_diffeomorphic_demons_impl(fixed, vector_field, vector_field, stepper,
                                          iterations - 1, sigma);
}

template <class T, unsigned int D, class STEPPER>
hoNDArray<vector_td<T, D>>
base_diffeomorphic_demons(const hoNDArray<T>& fixed, hoNDArray<vector_td<T, D>> vector_field,
                          STEPPER& stepper, unsigned int iterations, float sigma) {

    auto current_fixed = deform_image(fixed, vector_field);

    auto update_field = stepper(current_fixed);
    update_field = vector_field_exponential(update_field);
    vector_field = compose_fields(update_field, vector_field);
    vector_field = gaussian_filter(vector_field, sigma);
    return base_diffeomorphic_demons_impl(fixed, std::move(vector_field), std::move(update_field),
                                          stepper, iterations - 1, sigma);
}

} // namespace

template <class T, unsigned int D>
hoNDArray<vector_td<T, D>>
Gadgetron::Registration::diffeomorphic_demons(const hoNDArray<T>& fixed, const hoNDArray<T>& moving,
                                              unsigned int iterations, float sigma, float step_size,
                                              float noise_sigma) {
    if (fixed.dimensions() != moving.dimensions())
        throw std::runtime_error("Fixed and moving images have different sizes, "
                                 "which is not currently supported");
    auto stepper = DemonStep<T, D>{moving, 1 / step_size, 1e-6f, noise_sigma};
    return base_diffeomorphic_demons(fixed, stepper, iterations, sigma);
}

template <class T, unsigned int D>
hoNDArray<vector_td<T, D>>
Registration::diffeomorphic_demons(const hoNDArray<T>& fixed, const hoNDArray<T>& moving,
                                   hoNDArray<vector_td<T, D>> vector_field, unsigned int iterations,
                                   float sigma, float step_size, float noise_sigma) {

    if (fixed.dimensions() != moving.dimensions())
        throw std::runtime_error("Fixed and moving images have different sizes, "
                                 "which is not currently supported");
    if (fixed.dimensions() != vector_field.dimensions())
        throw std::runtime_error("Input vector field has mismatching dimensions for image size");

    auto stepper = DemonStep<T, D>{moving, 1 / step_size, 1e-6f, noise_sigma};
    return base_diffeomorphic_demons(fixed, std::move(vector_field), stepper, iterations, sigma);
}

template hoNDArray<vector_td<float, 2>>
Gadgetron::Registration::diffeomorphic_demons(const hoNDArray<float>&, const hoNDArray<float>&,
                                              unsigned int, float, float, float);
template hoNDArray<vector_td<float, 2>>
Gadgetron::Registration::diffeomorphic_demons(const hoNDArray<float>&, const hoNDArray<float>&,
                                              hoNDArray<vector_td<float, 2>>, unsigned int, float,
                                              float, float);
template <class T, unsigned int D>
hoNDArray<vector_td<T, D>> Gadgetron::Registration::ngf_diffeomorphic_demons(
    const hoNDArray<T>& fixed, const hoNDArray<T>& moving, unsigned int iterations, float sigma,
    float step_size, float gradient_eps) {
    if (fixed.dimensions() != moving.dimensions())
        throw std::runtime_error("Fixed and moving images have different sizes, "
                                 "which is not currently supported");
    auto stepper = NGFDemonsStep<T, D>(moving, 1 / step_size, 1e-6f, gradient_eps);
    return base_diffeomorphic_demons(fixed, stepper, iterations, sigma);
}

template <class T, unsigned int D>
hoNDArray<vector_td<T, D>> Gadgetron::Registration::ngf_diffeomorphic_demons(
    const hoNDArray<T>& fixed, const hoNDArray<T>& moving, hoNDArray<vector_td<T, D>> vector_field,
    unsigned int iterations, float sigma, float step_size, float gradient_eps) {
    if (fixed.dimensions() != moving.dimensions())
        throw std::runtime_error("Fixed and moving images have different sizes, "
                                 "which is not currently supported");
    auto stepper = NGFDemonsStep<T, D>(moving, 1 / step_size, 1e-6f, gradient_eps);
    return base_diffeomorphic_demons(fixed, std::move(vector_field),stepper, iterations, sigma);
}


template hoNDArray<vector_td<float, 2>>
Gadgetron::Registration::ngf_diffeomorphic_demons(const hoNDArray<float>&, const hoNDArray<float>&,
                                              unsigned int, float, float, float);
template hoNDArray<vector_td<float, 2>>
Gadgetron::Registration::ngf_diffeomorphic_demons(const hoNDArray<float>&, const hoNDArray<float>&,
                                              hoNDArray<vector_td<float, 2>>, unsigned int, float,
                                              float, float);

namespace {

template <class T, unsigned int D, class REGISTRATION>
hoNDArray<vector_td<T, D>> multi_scale_registration(const hoNDArray<T>& fixed,
                                                    const hoNDArray<T>& moving, unsigned int levels,
                                                    REGISTRATION reg, vector_td<T, D> tag) {
    if (levels <= 1)
        return reg(fixed, moving);
    auto fixed_pyramid = std::vector<hoNDArray<T>>{fixed};
    auto moving_pyramid = std::vector<hoNDArray<T>>{moving};

    for (int i = 0; i < int(levels) - 1; i++) {
        fixed_pyramid.push_back(downsample<T, D>(fixed_pyramid.back()));
        moving_pyramid.push_back(downsample<T, D>(moving_pyramid.back()));
    }

    auto current_vfield = reg(fixed_pyramid.back(), moving_pyramid.back());
    current_vfield = upsample<vector_td<T, D>, D>(current_vfield);
    current_vfield *= vector_td<T, D>(2);

    for (int i = fixed_pyramid.size() - 2; i > 0; i--) {
        auto update_field = reg(deform_image(fixed_pyramid[i], current_vfield), moving_pyramid[i]);
        current_vfield = compose_fields(update_field, current_vfield);
        current_vfield = upsample<vector_td<T, D>, D>(current_vfield);
        current_vfield *= vector_td<T, D>(2);
    }
    auto update_field = reg(deform_image(fixed, current_vfield), moving);
    return compose_fields(update_field, current_vfield);
}
template <class T, unsigned int D, class REGISTRATION>
hoNDArray<vector_td<T, D>> multi_scale_registration(const hoNDArray<T>& fixed,
                                                    const hoNDArray<T>& moving,
                                                    const hoNDArray<vector_td<T, D>>& vector_field,
                                                    unsigned int levels, REGISTRATION&& reg) {

    auto new_fixed = deform_image(fixed, vector_field);
    auto update_field = multi_scale_registration(
        new_fixed, moving, levels, std::forward<REGISTRATION>(reg), vector_td<T, D>{});

    return compose_fields(update_field, vector_field);
}

} // namespace

template <class T, unsigned int D>
hoNDArray<vector_td<T, D>> Gadgetron::Registration::multi_scale_diffeomorphic_demons(
    const hoNDArray<T>& fixed, const hoNDArray<T>& moving, unsigned int levels,
    unsigned int iterations, float sigma, float step_size, float noise_sigma) {
    return multi_scale_registration(
        fixed, moving, levels,
        [&](const auto& f, const auto& m) {
            return diffeomorphic_demons<T, D>(f, m, iterations, sigma, step_size, noise_sigma);
        },
        vector_td<T, D>{});
}

template hoNDArray<vector_td<float, 2>> Gadgetron::Registration::multi_scale_diffeomorphic_demons(
    const hoNDArray<float>& fixed, const hoNDArray<float>& moving, unsigned int levels,
    unsigned int iterations, float sigma, float step_size, float noise_sigma);


template <class T, unsigned int D>
hoNDArray<vector_td<T, D>> Gadgetron::Registration::multi_scale_ngf_diffeomorphic_demons(
    const hoNDArray<T>& fixed, const hoNDArray<T>& moving, unsigned int levels,
    unsigned int iterations, float sigma, float step_size, float gradient_eps) {
    return multi_scale_registration(
        fixed, moving, levels,
        [&](const auto& f, const auto& m) {
          return ngf_diffeomorphic_demons<T, D>(f, m, iterations, sigma, step_size, gradient_eps);
        },
        vector_td<T, D>{});
}

template hoNDArray<vector_td<float, 2>> Gadgetron::Registration::multi_scale_ngf_diffeomorphic_demons(
    const hoNDArray<float>& fixed, const hoNDArray<float>& moving, unsigned int levels,
    unsigned int iterations, float sigma, float step_size, float gradient_eps);
