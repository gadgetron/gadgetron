#pragma once

#include <algorithm>
#include <range/v3/view.hpp>

namespace {

struct T1Data {
    hoNDArray<std::complex<float>> data;
    std::vector<float> TI_values;
};

std::vector<hoNDArray<vector_td<float, 2>>> group_t1_registration(const std::vector<T1Data>& t1data,
                                                                  int iterations) {
    using namespace ranges;
    auto vector_fields = t1data | views::transform([iterations](const T1Data& data) {
                             return T1::t1_moco_cmr(data.data, data.TI_values, iterations);
                         }) |
                         to<std::vector>();

    return vector_fields;
}

int find_max_TI_group(const std::vector<T1Data>& t1data) {
    using namespace ranges;
    auto max_TIs = t1data | view::transform([](const auto& data) {
                       return *std::max_element(data.TI_values.begin(), data.TI_values.end());
                   }) |
                   to<std::vector>();
    auto max_TI_group =
        std::max_element(ranges::begin(max_TIs), ranges::end(max_TIs)) - ranges::begin(max_TIs);

    return max_TI_group;
}

hoNDArray<vector_td<float, 2>>
register_groups(const std::vector<T1Data>& t1data,
                const std::vector<hoNDArray<vector_td<float, 2>>>& vector_fields) {
    using namespace ranges;
    using namespace Gadgetron::Indexing;

    const auto deformed_groups =
        view::zip_with([](const auto& data,
                          const auto& vfields) { return T1::deform_groups(data.data, vfields); },
                       t1data, vector_fields) |
        to<std::vector>();

    auto t1_values = t1data | view::transform([](const auto& t1_data) -> const std::vector<float>& {
                         return t1_data.TI_values;
                     });

    auto max_TI_frames =
        view::zip_with(
            [](const auto& deformed, const auto& TI) {
                auto max_TI_index = std::max_element(std::begin(TI), std::end(TI)) - std::begin(TI);
                return hoNDArray<std::complex<float>>(deformed(slice, slice, max_TI_index));
            },
            deformed_groups, t1_values) |
        to<std::vector>();

    int max_TI_group = find_max_TI_group(t1data);
    const auto reference_frame = abs(max_TI_frames[max_TI_group]);

    int total_number_of_frames = ranges::accumulate(
        t1data | view::transform([](const auto& data) { return data.data.get_size(2); }), 0,
        std::plus());

    hoNDArray<vector_td<float, 2>> result(t1data.front().data.get_size(0),
                                          t1data.front().data.get_size(1), total_number_of_frames);
    result.fill(vector_td<float, 2>{0, 0});
    size_t current_result_frame = 0;

    for (int i = 0; i < deformed_groups.size(); i++) {
                if (i == max_TI_group) {
//        if (true) {
            for (int inner = 0; inner < vector_fields[i].get_size(2); inner++) {
                result(slice, slice, current_result_frame) = vector_fields[i](slice, slice, inner);
                current_result_frame++;
            }
            continue;
        }

        auto vfield = T1::register_groups_CMR(abs(max_TI_frames[i]), reference_frame);
        for (int inner = 0; inner < vfield.get_size(2); inner++) {
            result(slice, slice, current_result_frame) = Registration::compose_fields<float, 2>(
                vfield(slice, slice, inner), vector_fields[i](slice, slice, inner));
            current_result_frame++;
        }
    }

    return result;
}

hoNDArray<std::complex<float>> t1_cmr_groupwise_registration(const std::vector<T1Data>& t1data,
                                                             int iterations) {

    auto vector_field_groups = group_t1_registration(t1data, iterations);

    auto composed_vfield = register_groups(t1data, vector_field_groups);

    auto combined_data = concat_along_dimension(
        t1data | ranges::view::transform([](const auto& t1) { return t1.data; }), 2);

    auto deformed = T1::deform_groups(combined_data, composed_vfield);
    return deformed;
}

} // namespace