//
// Created by dch on 09/04/18.
//

#ifndef IMAGEGRAPH_H
#define IMAGEGRAPH_H
#pragma once
#
#include <boost/graph/graph_traits.hpp>

#include <boost/iterator/counting_iterator.hpp>
#include <vector>

//template<class T> T& get(std::vector<T>& vec, size_t i){ return vec[i];}
//
//template<class T, class R>
//void put(std::vector<T> &vec, size_t i, R val) {
//    vec[i] = val;
//}

#include <boost/graph/properties.hpp>
#include <numeric>
#include "vector_td_utilities.h"

namespace Gadgetron {

    template<unsigned int D = 2>
    class ImageGraph {
    public:

        constexpr static char edges_per_vertex = 2 * D + 2; // 8 per vertex + source & sink
        constexpr static char source_offset = 2 * D;
        constexpr static char sink_offset = 2 * D + 1;
        static Gadgetron::vector_td<int, D> index_to_offset[2 * D];
        using vertex_descriptor = size_t;

        using vertices_size_type = size_t;

        using edges_size_type = size_t;
        using edge_descriptor = size_t;

        class traversal_category : public boost::vertex_list_graph_tag,
                                   public boost::edge_list_graph_tag,
                                   public boost::incidence_graph_tag {
        };


        using vertex_iterator = boost::counting_iterator<size_t>;
        using edge_iterator = boost::counting_iterator<size_t>;




        template<class... Args> ImageGraph(Args... args): ImageGraph(vector_td<int,D>(args...)) {

        }

        ImageGraph(const vector_td<int,D>& dims): dims_(dims){
            num_image_vertices_ = std::accumulate(std::begin(dims_),std::end(dims_),1,std::multiplies<int>());
            num_vertices_ = num_image_vertices_ + 2;
            num_edges_ = (edges_per_vertex + 2) * num_image_vertices_;
            source_vertex = num_image_vertices_;
            sink_vertex = source_vertex + 1;

            setup_reverse_edge_map();
            reset();
        }


        void reset() {
            edge_capacity_map = std::vector<float>(num_edges_, 0);
            edge_residual_capicty = std::vector<float>(num_edges_, 0);
            color_map = std::vector<boost::default_color_type>(num_vertices_, boost::default_color_type::gray_color);
            vertex_distance = std::vector<float>(num_vertices_, 0);
            vertex_predecessor = std::vector<vertex_descriptor>(num_vertices_, 0);

        }

        vertex_iterator vertex_begin() const {
            return vertex_iterator(0);
        }

        vertex_iterator vertex_end() const {
            return vertex_iterator(num_vertices_);
        }

        size_t num_vertices() const {
            return num_vertices_;
        }


        edge_iterator edge_begin() const {
            return edge_iterator(0);
        }

        edge_iterator edge_end() const {
            return edge_iterator(num_edges_);
        }


        vertex_descriptor source(edge_descriptor e) const {

            size_t normal_vertex_id = e / edges_per_vertex;
            if (normal_vertex_id < num_image_vertices_) {
                return normal_vertex_id;
            }

            size_t remainder = e - num_image_vertices_ * edges_per_vertex;
            if (remainder < num_image_vertices_) {
                return source_vertex;
            } else {
                return sink_vertex;
            }
        }


        const vertex_descriptor target(edge_descriptor e) const {


            size_t normal_vertex_id = e / edges_per_vertex;
            if (normal_vertex_id < num_image_vertices_) {
                size_t index_offset = e - normal_vertex_id * edges_per_vertex;
                if (index_offset < edges_per_vertex - 2) {
                    vector_td<int,D> offset = index_to_offset[index_offset];
                    vector_td<int,D> co = Gadgetron::idx_to_co<int>(normal_vertex_id, dims_);
                    vector_td<int,D> co2 = (co + offset + dims_) % dims_;
                    return Gadgetron::co_to_idx(co2, dims_);
                } else {
                    if (index_offset == edges_per_vertex - 2)
                        return source_vertex;
                    return sink_vertex;
                }
            } else {
                return e % num_image_vertices_;
            }
        }


        std::pair<edge_iterator, edge_iterator> out_edges(vertex_descriptor v) const {
            if (v < num_image_vertices_) {
                return std::make_pair(edge_iterator(v * edges_per_vertex), edge_iterator((v + 1) * edges_per_vertex));
            } else {
                if (v == source_vertex) {
                    edge_descriptor source_start = num_image_vertices_ * edges_per_vertex;
                    return std::make_pair(edge_iterator(source_start),
                                          edge_iterator(source_start + num_image_vertices_));
                } else {
                    edge_descriptor source_start = num_image_vertices_ * (edges_per_vertex + 1);
                    return std::make_pair(edge_iterator(source_start),
                                          edge_iterator(source_start + num_image_vertices_));
                }
            }

        }
//    size_t out_degree(vertex_descriptor v) const;

        size_t out_degree(vertex_descriptor v) const {
            if (v < num_image_vertices_) {
                return edges_per_vertex;
            } else {
                return num_image_vertices_;
            }
        }


        edge_descriptor reverse(edge_descriptor e) const {

            vertex_descriptor sv = source(e);
            vertex_descriptor tv = target(e);
            auto res = edge(tv,sv);
            assert(res.second);
            return res.first;

        }


        std::pair<edge_descriptor, bool> edge(vertex_descriptor u, vertex_descriptor v) const {



            if ((u == source_vertex || u == sink_vertex) && (v == source_vertex || v == sink_vertex)) {
                return std::make_pair(edge_descriptor(0), false);
            }
            edge_descriptor result;

            if (u == source_vertex) {
                result = num_image_vertices_ * edges_per_vertex + v;
            } else if (u == sink_vertex) {
                result = num_image_vertices_ * (edges_per_vertex + 1) + v;

            } else if (v == source_vertex) {
                result = u * edges_per_vertex + edges_per_vertex - 2;
            } else if (v == sink_vertex) {
                result = u * edges_per_vertex + edges_per_vertex - 1;
            } else {

                vector_td<int, D> co1 = Gadgetron::idx_to_co<int,D>(u, dims_);
                vector_td<int, D> co2 = Gadgetron::idx_to_co<int,D>(v, dims_);

                auto diff = co2 - co1;
                for (int i = 0; i < D; i++){ //Fix wrapping boundary conditions
                    if (diff[i] != 0) {
                        if (diff[i] == (dims_[i] - 1)) diff[i] = -1;
                        if (diff[i] == (-dims_[i] + 1)) diff[i] = 1;
                    }
                }

                auto res = std::accumulate(std::begin(diff),std::end(diff),0, [](auto v1, auto v2){ return v1+abs(v2);});
                if (res != 1) return std::make_pair(edge_descriptor(0),false);


                int offset;
                for (auto i = 0; i < D; i++){
                    if (diff[i] == -1) offset = i*2;
                    else if (diff[i] == 1 ) offset = i*2+1;
                }

                result = u * edges_per_vertex + offset;
            }

            return std::make_pair(result,true);

        };


        size_t num_edges() const {
            return num_edges_;
        }

        edge_descriptor edge_from_source(vertex_descriptor idx){
            return source_vertex*edges_per_vertex+idx;
        }

        edge_descriptor edge_to_sink(vertex_descriptor idx){
            return idx*edges_per_vertex+sink_offset;
        }




        std::vector<float> edge_capacity_map;
        std::vector<float> edge_residual_capicty;
        std::vector<boost::default_color_type> color_map;
        std::vector<float> vertex_distance;
        std::vector<vertex_descriptor> vertex_predecessor;
        std::vector<edge_descriptor> reverse_edge_map;

        boost::identity_property_map vertex_index_map;
//    ReverseEdgeMap reverse_edge_map;

        vertex_descriptor source_vertex;
        vertex_descriptor sink_vertex;
    private:

        void setup_reverse_edge_map() {
            reverse_edge_map = std::vector<edge_descriptor>(num_edges_);
            for (edge_descriptor edge = 0; edge < num_edges_; edge++) {
                reverse_edge_map[edge] = reverse(edge);
            }

        }

        size_t num_vertices_;
        size_t num_image_vertices_;
        size_t num_edges_;

        Gadgetron::vector_td<int, D> dims_;


    };
}















template<unsigned int D> std::pair<typename Gadgetron::ImageGraph<D>::edge_descriptor ,bool> edge(typename Gadgetron::ImageGraph<D>::vertex_descriptor v1, typename Gadgetron::ImageGraph<D>::vertex_descriptor v2, const Gadgetron::ImageGraph<D>& g);

namespace boost {
    template<unsigned int D> std::pair<typename Gadgetron::ImageGraph<D>::vertex_iterator, typename Gadgetron::ImageGraph<D>::vertex_iterator > vertices(const Gadgetron::ImageGraph<D>& g);

    template<unsigned int D> float* get(boost::edge_capacity_t, Gadgetron::ImageGraph<D> &g);

    template<unsigned int D>  float* get(boost::edge_residual_capacity_t, Gadgetron::ImageGraph<D> &g);

    template<unsigned int D>  boost::default_color_type* get(boost::vertex_color_t, Gadgetron::ImageGraph<D> &g);

    template<unsigned int D>  float* get(boost::vertex_distance_t, Gadgetron::ImageGraph<D> &g);

    template<unsigned int D>  const identity_property_map &get(boost::vertex_index_t, const Gadgetron::ImageGraph<D> &g);

    template<unsigned int D>  size_t* get(boost::edge_reverse_t, Gadgetron::ImageGraph<D> &g);

    template<unsigned int D>  size_t* get(boost::vertex_predecessor_t, Gadgetron::ImageGraph<D> &g);
    template<unsigned int D>  const float* get(boost::edge_capacity_t, const Gadgetron::ImageGraph<D> &g);

    template<unsigned int D> const float* get(boost::edge_residual_capacity_t, const Gadgetron::ImageGraph<D> &g);

    template<unsigned int D>  const boost::default_color_type* get(boost::vertex_color_t, const Gadgetron::ImageGraph<D> &g);

    template<unsigned int D> const float* get(boost::vertex_distance_t, const Gadgetron::ImageGraph<D> &g);



    template<unsigned int D> const size_t* get(boost::edge_reverse_t, const Gadgetron::ImageGraph<D> &g);
    template<unsigned int D> const size_t* get(boost::vertex_predecessor_t, const Gadgetron::ImageGraph<D> &g);

    template<unsigned int D>  float get(boost::edge_capacity_t, const Gadgetron::ImageGraph<D> &g, size_t i);

    template<unsigned int D> size_t num_vertices(const Gadgetron::ImageGraph<D>& g);

    template<unsigned int D> std::pair<typename Gadgetron::ImageGraph<D>::edge_iterator, typename Gadgetron::ImageGraph<D>::edge_iterator> edges(const Gadgetron::ImageGraph<D>& g);

    template<unsigned int D> std::pair<typename Gadgetron::ImageGraph<D>::edge_iterator,typename  Gadgetron::ImageGraph<D>::edge_iterator> out_edges(typename Gadgetron::ImageGraph<D>::vertex_descriptor v, const Gadgetron::ImageGraph<D>& g);
    template<unsigned int D> size_t out_degree(typename Gadgetron::ImageGraph<D>::vertex_descriptor v, const Gadgetron::ImageGraph<D>& g);
    template <unsigned int D> size_t num_edges(const Gadgetron::ImageGraph<D>& g);


    template<unsigned int D> typename Gadgetron::ImageGraph<D>::vertex_descriptor source(typename Gadgetron::ImageGraph<D>::edge_descriptor e, const Gadgetron::ImageGraph<D>& g);

    template<unsigned int D>  typename Gadgetron::ImageGraph<D>::vertex_descriptor target(typename Gadgetron::ImageGraph<D>::edge_descriptor e, const Gadgetron::ImageGraph<D>& g);

}


namespace boost {



//    typedef Gadgetron::ImageGraph* image_graph_ptr;
//    typedef const Gadgetron::ImageGraph* image_const_graph_ptr;



    template <unsigned  int D> struct graph_traits<Gadgetron::ImageGraph<D>> {
        typedef typename Gadgetron::ImageGraph<D>::vertex_descriptor vertex_descriptor;
        typedef typename Gadgetron::ImageGraph<D>::edge_descriptor edge_descriptor;
        typedef typename Gadgetron::ImageGraph<D>::edge_iterator out_edge_iterator;
        typedef void in_edge_iterator;

        typedef typename Gadgetron::ImageGraph<D>::vertex_iterator vertex_iterator;
        typedef typename Gadgetron::ImageGraph<D>::edge_iterator edge_iterator;
        typedef size_t vertices_size_type;
        typedef size_t edges_size_type;
        typedef size_t degree_size_type;
        typedef directed_tag directed_category;
        typedef typename Gadgetron::ImageGraph<D>::traversal_category traversal_category;

        typedef disallow_parallel_edge_tag edge_parallel_category;
        static edge_descriptor null_vertex(){ return std::numeric_limits<size_t>::max();}



    };
/*
    template <> struct graph_traits<image_const_graph_ptr> {
        typedef Gadgetron::ImageGraph::vertex_descriptor vertex_descriptor;
        typedef Gadgetron::ImageGraph::edge_descriptor edge_descriptor;
        typedef Gadgetron::ImageGraph::edge_iterator out_edge_iterator;
        typedef void in_edge_iterator;

        typedef Gadgetron::ImageGraph::vertex_iterator vertex_iterator;
        typedef Gadgetron::ImageGraph::edge_iterator edge_iterator;
        typedef size_t vertices_size_type;
        typedef size_t edge_size_type;
        typedef size_t degree_size_type;
        typedef directed_tag directed_category;
        typedef Gadgetron::ImageGraph::traversal_category traversal_category;


    };*/


    template<unsigned int D> struct property_map<Gadgetron::ImageGraph<D>,edge_capacity_t>{ typedef float* type; typedef const float* const_type; };
    template<unsigned int D> struct property_map<Gadgetron::ImageGraph<D>,edge_residual_capacity_t>{ typedef float* type; typedef const float* const_type;};
    template<unsigned int D> struct property_map<Gadgetron::ImageGraph<D>,edge_reverse_t >{ typedef size_t* type; typedef const size_t* const_type;};

    template<unsigned int D> struct property_map<Gadgetron::ImageGraph<D>,vertex_color_t >{ typedef default_color_type* type; typedef const default_color_type* const_type;};
    template<unsigned int D> struct property_map<Gadgetron::ImageGraph<D>,vertex_distance_t >{ typedef float* type; typedef const float* const_type;};
    template<unsigned int D> struct property_map<Gadgetron::ImageGraph<D>,vertex_index_t >{ typedef identity_property_map type; typedef identity_property_map const_type;};
    template<unsigned int D> struct property_map<Gadgetron::ImageGraph<D>,vertex_predecessor_t >{ typedef size_t* type; typedef const size_t* const_type;};

}

#include "ImageGraph.hxx"

#endif