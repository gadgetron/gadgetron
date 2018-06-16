//
// Created by dch on 09/04/18.
//

#include "ImageGraph.h"


 template<unsigned int D>
    std::pair<typename Gadgetron::ImageGraph<D>::edge_descriptor, bool>
    edge(typename Gadgetron::ImageGraph<D>::vertex_descriptor v1, typename Gadgetron::ImageGraph<D>::vertex_descriptor v2,
         const Gadgetron::ImageGraph<D> &g) {
        return g.edge(v1, v2);

    };

namespace boost {
    template<unsigned int D>
    std::pair<typename Gadgetron::ImageGraph<D>::vertex_iterator, typename Gadgetron::ImageGraph<D>::vertex_iterator>
    vertices(const Gadgetron::ImageGraph<D> &g) {

        return std::make_pair(g.vertex_begin(), g.vertex_end());

    };




    template<unsigned int D>
    size_t num_vertices(const Gadgetron::ImageGraph<D> &g) {
        return g.num_vertices();
    }


    template<unsigned int D>
    size_t num_edges(const Gadgetron::ImageGraph<D> &g) {
        return g.num_edges();
    }


    template<unsigned int D>
    typename Gadgetron::ImageGraph<D>::vertex_descriptor
    source(typename Gadgetron::ImageGraph<D>::edge_descriptor e, const Gadgetron::ImageGraph<D> &g) {
        return g.source(e);
    }


    template<unsigned int D>
    typename Gadgetron::ImageGraph<D>::vertex_descriptor
    target(typename Gadgetron::ImageGraph<D>::edge_descriptor e, const Gadgetron::ImageGraph<D> &g) {
        return g.target(e);
    }


    template<unsigned int D>
    std::pair<typename Gadgetron::ImageGraph<D>::edge_iterator, typename Gadgetron::ImageGraph<D>::edge_iterator>
    out_edges(typename Gadgetron::ImageGraph<D>::vertex_descriptor v, const Gadgetron::ImageGraph<D> &g) {
        return g.out_edges(v);
    };


    template<unsigned int D>
    size_t out_degree(typename Gadgetron::ImageGraph<D>::vertex_descriptor v, const Gadgetron::ImageGraph<D> &g) {
        return g.out_degree(v);
    }


    template<unsigned int D>
    std::pair<typename Gadgetron::ImageGraph<D>::edge_iterator, typename Gadgetron::ImageGraph<D>::edge_iterator>
    edges(const Gadgetron::ImageGraph<D> &g) {
        return std::make_pair(g.edge_begin(), g.edge_end());
    };


    template<unsigned int D>
    float *get(edge_capacity_t, Gadgetron::ImageGraph<D> &g) {
        return g.edge_capacity_map.data();
    }

    template<unsigned int D>
    float *get(edge_residual_capacity_t, Gadgetron::ImageGraph<D> &g) {
        return g.edge_residual_capicty.data();
    }

    template<unsigned int D>
    default_color_type *get(vertex_color_t, Gadgetron::ImageGraph<D> &g) {
        return g.color_map.data();
    }

    template<unsigned int D>
    float *get(vertex_distance_t, Gadgetron::ImageGraph<D> &g) {
        return g.vertex_distance.data();
    }

    template<unsigned int D>
    const identity_property_map &get(vertex_index_t, const Gadgetron::ImageGraph<D> &g) {
        return g.vertex_index_map;
    }

    template<unsigned int D>
    size_t *get(edge_reverse_t, Gadgetron::ImageGraph<D> &g) {
        return g.reverse_edge_map.data();
    }

    template<unsigned int D>
    size_t *get(vertex_predecessor_t, Gadgetron::ImageGraph<D> &g) {
        return g.vertex_predecessor.data();

    }

    template<unsigned int D>
    const float *get(edge_capacity_t, const Gadgetron::ImageGraph<D> &g) {
        return g.edge_capacity_map.data();
    }

    template<unsigned int D>
    const float *get(edge_residual_capacity_t, const Gadgetron::ImageGraph<D> &g) {
        return g.edge_residual_capicty.data();
    }

    template<unsigned int D>
    const default_color_type *get(vertex_color_t, const Gadgetron::ImageGraph<D> &g) {
        return g.color_map.data();
    }

    template<unsigned int D>
    const float *get(vertex_distance_t, const Gadgetron::ImageGraph<D> &g) {
        return g.vertex_distance.data();
    }


    template<unsigned int D>
    const size_t *get(edge_reverse_t, const Gadgetron::ImageGraph<D> &g) {
        return g.reverse_edge_map.data();
    }

    template<unsigned int D>
    const size_t *get(vertex_predecessor_t, const Gadgetron::ImageGraph<D> &g) {
        return g.vertex_predecessor.data();
    }

    template<unsigned int D>
    float get(edge_capacity_t, Gadgetron::ImageGraph<D> &g, size_t i) {
        return g.edge_capacity_map[i];
    }
}

