//
// Created by dchansen on 5/24/18.
//

#include <ismrmrd/xml.h>
#include <unordered_map>
#include <numeric>
#include <boost/range/iterator_range.hpp>
#include <boost/iterator/transform_iterator.hpp>
#include "ImageAccumulatorGadget.h"


using namespace Gadgetron;


Gadgetron::ImageAccumulatorGadget::ImageAccumulatorGadget() {




}


static size_t image_dimension_from_string(std::string name){
    if (name == "X") return 0;
    if (name == "Y") return 1;
    if (name == "Z") return 2;
    if (name == "CHA") return 3;
    if (name == "N") return 4;
    if (name == "S") return 5;
    if (name == "LOC") return 6;

    throw std::runtime_error("Name " + name + " does not match an image dimension");

}

static size_t header_dimension_from_string(std::string name){
    if (name == "N") return 0;
    if (name == "S") return 1;
    if (name == "LOC") return 2;

    throw std::runtime_error("Name " + name + " does not match a header dimension");

}

template<class T> auto Gadgetron::ImageAccumulatorGadget::extract_value(T &val) {

    auto dimension = accumulate_dimension.value();
    if (dimension == "average") return val.average;
    if (dimension == "slice") return val.slice;
    if (dimension == "contrast") return val.contrast;
    if (dimension == "phase") return val.phase;
    if (dimension == "repetition") return val.repetition;
    if (dimension == "set") return val.set;

    throw std::runtime_error("Unknown dimension type " + dimension);
}





int Gadgetron::ImageAccumulatorGadget::process_config(ACE_Message_Block *mb) {

    ISMRMRD::IsmrmrdHeader h;
    ISMRMRD::deserialize(mb->rd_ptr(),h);

    auto limits_opt = extract_value(h.encoding[0].encodingLimits);
    if (!limits_opt) throw std::runtime_error("Encoding limits not set in data for dimension " + accumulate_dimension.value());
    ISMRMRD::Limit limits = limits_opt.get();



    required_values = std::vector<uint16_t>();
    for (auto val = limits.minimum; val <= limits.maximum; val++) required_values.push_back(val);




    return GADGET_OK;


}

//bool Gadgetron::ImageAccumulatorGadget::same_size(std::vector<Gadgetron::IsmrmrdImageArray>& values){
//
//    if (values.size() <= 1)
//        return true;
//
//    auto& im = values.front().data_;
//    return std::all_of(values.begin()+1,values.end(),[&](auto & im2 ){return im.dimensions_equal(im2.data_);});
//
//}

namespace {


    template<class T, int DIM> class combiner {

    public:
        static void
        combine_along(T *out, const T *in, const std::vector<size_t> dims, const std::vector<size_t> &out_stride,
                      const std::vector<size_t> in_stride, size_t combine_dim) {
            for (int i = 0; i < dims[DIM]; i++) {
                combiner<T,DIM-1>::combine_along(out + out_stride[DIM] * i, in + in_stride[DIM] * i, dims, out_stride,
                                          in_stride,
                                          combine_dim);
            }

        }

    };

    template<class T> class combiner<T,0> {
    public:
        static void
        combine_along(T *out, const T *in, const std::vector<size_t> dims, const std::vector<size_t> &out_stride,
                      const std::vector<size_t> in_stride, size_t combine_dim) {
            for (int i = 0; i < dims[0]; i++) {
                out[i * out_stride[0]] = in[i];
            }
        }
    };


            template<class T> std::vector<size_t> calculate_strides(const hoNDArray<T>& array) {
        auto dims = *array.get_dimensions();
        auto strides = std::vector<size_t>(dims.size(),1);

        std::partial_sum(dims.begin(),dims.end()-1,strides.begin()+1,std::multiplies<size_t>());

        return strides;
    }

    template<class RANGE, int DIMS> auto combine_arrays_along(RANGE &input_arrays,size_t combine_dim) {

        using T = typename decltype(*input_arrays.begin())::value_type;
        auto output_dims = *input_arrays.front().get_dimensions();
        output_dims[combine_dim] = 0;

        for (const auto & array : input_arrays) output_dims[combine_dim] += array.get_size(combine_dim);

        hoNDArray<T> out(output_dims);

        auto out_strides = calculate_strides(out);

        T* output_data = out.get_data_ptr();

        for (const auto & array : input_arrays) {
            auto dims = *array.get_dimensions();
            auto in_strides = calculate_strides(array);
            combiner<T,DIMS-1>::combine_along(output_data,array.get_data_ptr(),dims,out_strides,in_strides,combine_dim);
            output_data += out_strides[combine_dim-1]*dims[combine_dim];
        }

        return out;

    }

};






Gadgetron::IsmrmrdImageArray
Gadgetron::ImageAccumulatorGadget::combine_images(std::vector<Gadgetron::IsmrmrdImageArray>& images) {

//    if (!same_size(images)) throw std::runtime_error("Images do not have the same size");

    size_t combine_dimension_image = image_dimension_from_string(combine_along.value());
    size_t combine_dimension_header = header_dimension_from_string(combine_along.value());

    auto image_lambda = [](auto m) { return m.data_;};
    auto header_lambda = [](auto m) { return m.headers_;};

    //Using boost range and transform iterators to avoid copying from the vector

    auto image_range = boost::make_iterator_range(boost::make_transform_iterator(images.begin(),image_lambda),
            boost::make_transform_iterator(images.end(),image_lambda));

    auto header_range = boost::make_iterator_range(boost::make_transform_iterator(images.begin(),header_lambda),
                                             boost::make_transform_iterator(images.end(),header_lambda));

    auto ninjas = images[0].data_.get_number_of_dimensions();
    Gadgetron::IsmrmrdImageArray result =
            { combine_arrays_along<decltype(image_range), 7>(image_range,combine_dimension_image),
              combine_arrays_along<decltype(header_range),3>(header_range,combine_dimension_header),
               std::vector< ISMRMRD::MetaContainer>()};
    for (auto& im : images){
        result.meta_.insert(result.meta_.end(),im.meta_.begin(),im.meta_.end());
    }

    return result;
}

int Gadgetron::ImageAccumulatorGadget::process(Gadgetron::GadgetContainerMessage<Gadgetron::IsmrmrdImageArray> *m1) {


    auto& image_array = *m1->getObjectPtr();

    images.push_back(image_array);

    for (auto h : image_array.headers_)
        seen_values.emplace(extract_value(h));

    bool done = std::all_of(required_values.begin(),required_values.end(),
            [&](uint16_t val){return seen_values.count(val);});

    if (done){
        auto result_array = combine_images(images);
        images.clear();
        seen_values.clear();
        auto msg = new GadgetContainerMessage<IsmrmrdImageArray>(std::move(result_array));
        this->next()->putq(msg);

    }
    return GADGET_OK;
}


GADGET_FACTORY_DECLARE(ImageAccumulatorGadget)


