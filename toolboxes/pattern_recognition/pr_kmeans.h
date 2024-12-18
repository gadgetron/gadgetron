/** \file   pr_kmeans.h
    \brief  Implement kmeans and related algorithm
    \author Hui Xue
*/

#pragma once

#include "GadgetronTimer.h"
#include "ImageIOAnalyze.h"
#include "hoNDArray.h"

namespace Gadgetron {

// ======================================================================================
// kmeans class
// input: X, [P N] data array, N samples with P dimensions
// K : number of clusters to group data into
// max_iter: number of maximal k-means iterations
// replicates : number of times to reinitialize kmenas and reperform the clustering
//
// initial_method : Initialization method,
// 'sample' : Select K samples from N data points as the initial cluster centroid
// 'uniform' : First detect the range of N data points and select K centroids from the range of X
// 'cluster' : First randomly selected 20% of all N samples and perform kmeans using 'sample' method
// then, the resulting centroids are used for whole data kmeans
// 'kmeans++': perform the kmeans++ method, http://ilpubs.stanford.edu:8090/778/1/2006-13.pdf
//
// online update: the kmeans can optionally use the so-called "online" update. In this process, every data point is reallocated to all clusters and the
// delta change of adding or removing this point is computed; those moves which will reduce the total sum cost will be performed.
//
// output
// IDX : [N 1] array, indicating to which clusters every data sample belongs (first cluster has index 0)
// C : [P K], K centroids
// sumD : [K 1], sum of intra-cluster distances
//
// The usage of this class is one of the following two options:
// 1) Call one of the get_initial_guess_XXX methods and compute initial centroid
// 2) Call the run function which performs kmeans iterations and generate three outputs

template <typename T>
class kmeans
{
public:

    typedef kmeans<T> Self;

    typedef hoNDArray<T> ArrayType;
    typedef std::vector<T> VectorType;
    typedef std::vector<size_t> ClusterType;

    kmeans();
    virtual ~kmeans();

    // ======================================================================================
    /// parameter for kmeans
    // ======================================================================================

    // maximal number of iterations
    size_t max_iter_;

    // number of repeated trials
    size_t replicates_;

    // whether to perform on-line update
    bool perform_online_update_;

    // ======================================================================================
    /// parameter for debugging
    // ======================================================================================

    // verbose mode for more output messages
    bool verbose_;
    // debug folder
    std::string debug_folder_;
    // whether to perform timing
    bool perform_timing_;

    // clock for timing
    Gadgetron::GadgetronTimer gt_timer_local_;
    Gadgetron::GadgetronTimer gt_timer_;

    // exporter
    Gadgetron::ImageIOAnalyze gt_exporter_;

    // ======================================================================================
    // perform every steps
    // ======================================================================================

    /// provide initial guess for centroid
    /// C_for_initial: initialization centroids[K P num]
    virtual void get_initial_guess_sample(const ArrayType& X, size_t K, ArrayType& C_for_initial);
    virtual void get_initial_guess_uniform(const ArrayType& X, size_t K, ArrayType& C_for_initial);
    virtual void get_initial_guess_cluster(const ArrayType& X, size_t K, ArrayType& C_for_initial);
    virtual void get_initial_guess_kmeansplusplus(const ArrayType& X, size_t K, ArrayType& C_for_initial);

    /// compute kmeans
    virtual void run_replicates(const ArrayType& X, size_t K, const ArrayType& C_for_initial, ClusterType& IDX, ArrayType& C, VectorType& sumD_rep, T& sumD);
    virtual void run(const ArrayType& X, size_t K, const ArrayType& C_for_initial, ClusterType& IDX, ArrayType& C, T& sumD);

    /// compute distance vector
    /// D: [P N] distance from a point to its closest centroid
    void compute_dist(const ArrayType& X, const ClusterType& IDX, const ArrayType& C, ArrayType& D);

    /// compute norm of distance vector
    void compute_norm_dist(const ArrayType& D, ArrayType& D_norm);

    /// given the current centroids, update the IDX
    /// norm_C is the norm of centroid, dot(C,C,1)
    void update_IDX(const ArrayType& X, const ArrayType& C, const VectorType& norm_C, ClusterType& IDX);

    /// update centroids, given the IDX
    void update_centroid(const ArrayType& X, const ClusterType& IDX, ArrayType& C, VectorType& norm_C);

    /// chech whether there are empty clusters
    void has_empty_cluster(const ArrayType& C, const ClusterType& IDX, std::vector<size_t>& cluster_size, std::vector<size_t>& empty_clusters);

    /// for all points, find the most furthest point from its OWN cluster
    void find_lonely_point(const ClusterType& IDX, const ArrayType& D_norm, const std::vector<size_t>& cluster_size, size_t& lonely);

    /// check whether clustering results changed
    /// return true, if clustering results are changed
    bool is_clustering_changed(const ClusterType&prev_IDX, const ClusterType& IDX);

    /// online update
    /// For input, IDX and C stores current clustering results and centroids
    /// On return, IDX and C may be updated
    /// max_iter_ is used for online update
    void perform_online_update(const ArrayType& X, ClusterType& IDX, ArrayType& C, T& sumD);
};

}
