//
// Created by david on 6/7/2018.
//

#include <random>
#include "ImageGraph.h"
#include <boost/graph/boykov_kolmogorov_max_flow.hpp>
#include "graph_cut.h"


namespace {
    using namespace Gadgetron;
    static std::mt19937 rng_state(4242);


    template<unsigned int D>
    void update_regularization_edge(ImageGraph<D> &graph, const hoNDArray<uint16_t> &field_map,
                                    const hoNDArray<uint16_t> &proposed_field_map,
                                    const hoNDArray<float> &second_deriv, const size_t idx, const size_t idx2,
                                    const size_t edge_idx, float scaling) {

        int f_value1 = field_map[idx];
        int pf_value1 = proposed_field_map[idx];
        int f_value2 = field_map[idx2];
        int pf_value2 = proposed_field_map[idx2];
        int a = std::norm(f_value1 - f_value2);
        int b = std::norm(f_value1 - pf_value2);
        int c = std::norm(pf_value1 - f_value2);
        int d = std::norm(pf_value1 - pf_value2);

        float weight = b + c - a - d;

        assert(weight >= 0);
        float lambda = std::max(std::min(second_deriv[idx], second_deriv[idx2]), 0.0f) * scaling;
        weight *= lambda;

        assert(lambda >= 0);

        auto &capacity_map = graph.edge_capacity_map;

        capacity_map[edge_idx] += weight;
        {
            float aq = lambda * (c - a);

            if (aq > 0) {
                capacity_map[graph.edge_from_source(idx)] += aq;

            } else {
                capacity_map[graph.edge_to_sink(idx)] -= aq;
            }
        }

        {
            float aj = lambda * (d - c);
            if (aj > 0) {
                capacity_map[graph.edge_from_source(idx2)] += aj;

            } else {
                capacity_map[graph.edge_to_sink(idx2)] -= aj;
            }
        }


    }

    template<unsigned int D>
    ImageGraph<D> make_graph(const hoNDArray<uint16_t> &field_map, const hoNDArray<uint16_t> &proposed_field_map,
                             const hoNDArray<float> &residual_diff_map, const hoNDArray<float> &second_deriv) {

        const auto dims = vector_td<int,3>(field_map.get_size(0),field_map.get_size(1),field_map.get_size(2));

        vector_td<int,D> graph_dims;
        for (int i = 0; i < D; i++) graph_dims[i] = dims[i];

        ImageGraph<D> graph = ImageGraph<D>(graph_dims);

        auto &capacity_map = graph.edge_capacity_map;
        //Add regularization edges

        for (size_t kz = 0; kz < dims[2]; kz++) {
            for (size_t ky = 0; ky < dims[1]; ky++) {
                for (size_t kx = 0; kx < dims[0]; kx++) {
                    size_t idx = kz*dims[1]*dims[0]+ky * dims[0] + kx;


                    if (kx < (dims[0] - 1)) {
                        size_t idx2 = idx + 1;

                        update_regularization_edge(graph, field_map, proposed_field_map, second_deriv, idx, idx2,
                                                   graph.edge(idx, idx2).first, 1);
                    }


                    if (ky < (dims[1] - 1)) {
                        size_t idx2 = idx + dims[0];
                        update_regularization_edge(graph, field_map, proposed_field_map, second_deriv, idx, idx2,
                                                   graph.edge(idx, idx2).first, 1);
                    }

                    if (kz < (dims[2] - 1)) {
                        size_t idx2 = idx + dims[0]*dims[1];
                        update_regularization_edge(graph, field_map, proposed_field_map, second_deriv, idx, idx2,
                                                   graph.edge(idx, idx2).first, 1);
                    }

                    float residual_diff = residual_diff_map[idx];

                    if (residual_diff > 0) {
                        capacity_map[graph.edge_to_sink(idx)] += int(residual_diff);

                    } else {
                        capacity_map[graph.edge_from_source(idx)] -= int(residual_diff);
                    }

                }
            }
        }

        return graph;
    }

    template<unsigned int DIMS>
    std::vector<boost::default_color_type>
    graph_cut(const hoNDArray<uint16_t> &field_map_index, const hoNDArray<uint16_t> &proposed_field_map_index,
              const hoNDArray<float> &lambda_map, const hoNDArray<float> &residual_diff_map) {

        ImageGraph<DIMS> graph = make_graph<DIMS>(field_map_index, proposed_field_map_index, residual_diff_map,
                                                  lambda_map);

        float flow = boost::boykov_kolmogorov_max_flow(graph, graph.source_vertex, graph.sink_vertex);

        return std::move(graph.color_map);
    }

}
namespace Gadgetron {


    hoNDArray<uint16_t>
    update_field_map(const hoNDArray<uint16_t> &field_map_index, const hoNDArray<uint16_t> &proposed_field_map_index,
                     const hoNDArray<float> &residuals_map, const hoNDArray<float> &lambda_map) {


        hoNDArray<float> residual_diff_map(field_map_index.dimensions());
        const auto X = field_map_index.get_size(0);
        const auto Y = field_map_index.get_size(1);
        const auto Z = field_map_index.get_size(2);

        for (size_t kz = 0; kz < Z; kz++) {
            for (size_t ky = 0; ky < Y; ky++) {
                for (size_t kx = 0; kx < X; kx++) {
                    residual_diff_map(kx, ky,kz) = residuals_map(field_map_index(kx, ky,kz), kx, ky,kz) -
                                                residuals_map(proposed_field_map_index(kx, ky,kz), kx, ky,kz);


                }
            }
        }


        std::vector<boost::default_color_type> color_map;
        if (Z == 1) {
            color_map = graph_cut<2>(field_map_index, proposed_field_map_index, lambda_map,
                                     residual_diff_map);
        } else {
            color_map = graph_cut<3>(field_map_index, proposed_field_map_index, lambda_map, residual_diff_map);
        }



        auto result = field_map_index;
        size_t updated_voxels = 0;
        for (size_t i = 0; i < field_map_index.get_number_of_elements(); i++) {
            if (color_map[i] != boost::default_color_type::black_color) {
                updated_voxels++;
                result[i] = proposed_field_map_index[i];
            }
        }

        return result;

    }

}