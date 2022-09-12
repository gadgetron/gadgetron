/** \file   pr_kmeans.cpp
    \brief  Implement kmeans method with different initialization scheme
    \author Hui Xue
*/

#include "pr_kmeans.h"
#include "log.h"

#include "hoNDArray_reductions.h"
#include "hoNDArray_elemwise.h"
#include "hoNDArray_math.h"
#include "hoNDArray_linalg.h"

#include <boost/math/special_functions/sign.hpp>

#include <random>

namespace Gadgetron { 

template <typename T> 
kmeans<T>::kmeans()
{
    max_iter_ = 100;
    replicates_ = 10;
    perform_online_update_ = true;

    verbose_ = false;
    perform_timing_ = false;

    gt_timer_local_.set_timing_in_destruction(false);
    gt_timer_.set_timing_in_destruction(false);
}

template <typename T>
kmeans<T>::~kmeans()
{
}

template <typename T>
void kmeans<T>::get_initial_guess_sample(const ArrayType& X, size_t K, ArrayType& C_for_initial)
{
    try
    {
        if (this->perform_timing_) gt_timer_local_.start("get_initial_guess_sample");

        size_t P = X.get_size(0);
        size_t N = X.get_size(1);

        GADGET_CHECK_THROW(N>K);

        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dis(0, 1);

        C_for_initial.create(P, K, this->replicates_);
        Gadgetron::clear(C_for_initial);

        size_t n, k;
        for (n = 0; n < this->replicates_; n++)
        {
            for (k = 0; k < K; k++)
            {
                size_t ind = (size_t)(dis(gen)*N);
                if (ind >= N) ind = N - 1;
                memcpy(&C_for_initial(0, k, n), &X(0, ind), sizeof(T)*P);
            }
        }

        if (this->perform_timing_) gt_timer_local_.stop();
    }
    catch(...)
    {
        GERROR_STREAM("Exceptions happened in kmeans<T>::get_initial_guess_sample(...) ... ");
    }
}

template <typename T>
void kmeans<T>::get_initial_guess_uniform(const ArrayType& X, size_t K, ArrayType& C_for_initial)
{
    try
    {
        if (this->perform_timing_) gt_timer_local_.start("get_initial_guess_uniform");

        size_t P = X.get_size(0);
        size_t N = X.get_size(1);

        GADGET_CHECK_THROW(N>K);

        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dis(0, 1);

        C_for_initial.create(P, K, this->replicates_);
        Gadgetron::clear(C_for_initial);

        // find the range of X
        std::vector<T> xmin(P, 0);
        std::vector<T> xmax(P, 0);

        size_t p, n, k;

        for (p = 0; p<P; p++)
        {
            xmin[p] = X(p, 0);
            xmax[p] = xmin[p];
            for (n = 1; n < N; n++)
            {
                T v = X(p, n);
                if (v < xmin[p]) xmin[p] = v;
                if (v > xmax[p]) xmax[p] = v;
            }
        }

        for (n = 0; n < this->replicates_; n++)
        {
            for (k = 0; k < K; k++)
            {
                for (p = 0; p < P; p++)
                {
                    T v = xmin[p] + dis(gen) * (xmax[p] - xmin[p]);
                    C_for_initial(p, k, n) = v;
                }
            }
        }

        if (this->perform_timing_) gt_timer_local_.stop();
    }
    catch (...)
    {
        GERROR_STREAM("Exceptions happened in kmeans<T>::get_initial_guess_uniform(...) ... ");
    }
}

template <typename T>
void kmeans<T>::get_initial_guess_cluster(const ArrayType& X, size_t K, ArrayType& C_for_initial)
{
    try
    {
        if (this->perform_timing_) gt_timer_.start("get_initial_guess_cluster");

        size_t P = X.get_size(0);
        size_t N = X.get_size(1);

        GADGET_CHECK_THROW(N>K);

        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dis(0, 1);

        C_for_initial.create(P, K, this->replicates_);
        Gadgetron::clear(C_for_initial);

        // randomly pick 20% of data
        double ratio = 0.2;
        size_t M = (size_t)(ratio*N);
        while(M<K && ratio<=1.0)
        {
            ratio = ratio + 0.1;
            M = (size_t)(ratio*N);
        }

        if(M>=N)
        {
            GWARN_STREAM("Too few samples in the data array X ... ");
            this->get_initial_guess_sample(X, K, C_for_initial);
            return;
        }

        ArrayType X_subset;
        X_subset.create(P, M);

        ArrayType C_for_initial_subset;
        ClusterType IDX;
        ArrayType C;
        T sumD;

        size_t n, m;
        for (n = 0; n < this->replicates_; n++)
        {
            for (m = 0; m < M; m++)
            {
                size_t ind = (size_t)(dis(gen)*N);
                if (ind >= N) ind = N - 1;
                memcpy(&X_subset(0, m), &X(0, ind), sizeof(T)*P);
            }

            this->get_initial_guess_sample(X_subset, K, C_for_initial_subset);

            // call kmeans
            this->run(X_subset, K, C_for_initial_subset, IDX, C, sumD);

            memcpy(&C_for_initial(0, 0, n), C.begin(), sizeof(T)*K*P);
        }

        if (this->perform_timing_) gt_timer_.stop();
    }
    catch (...)
    {
        GERROR_STREAM("Exceptions happened in kmeans<T>::get_initial_guess_cluster(...) ... ");
    }
}

template <typename T>
void kmeans<T>::get_initial_guess_kmeansplusplus(const ArrayType& X, size_t K, ArrayType& C_for_initial)
{
    try
    {
        if (this->perform_timing_) gt_timer_.start("get_initial_guess_kmeansplusplus");

        size_t P = X.get_size(0);
        size_t N = X.get_size(1);

        GADGET_CHECK_THROW(N>K);

        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dis(0, 1);

        C_for_initial.create(P, K, this->replicates_);
        Gadgetron::clear(C_for_initial);

        ClusterType IDX;
        IDX.resize(N, 0);

        // find the first center
        ArrayType C;
        C.create(P, K);
        Gadgetron::clear(C);

        VectorType norm_C(K, 0);

        ArrayType D;
        D.create(P, N);

        ArrayType D_norm;
        D_norm.create(N);

        ArrayType cumsum_D_norm;
        cumsum_D_norm.create(N);

        ArrayType CX;

        size_t n, i, t, s;

        for (n = 0; n < this->replicates_; n++)
        {
            size_t ind = (size_t)(dis(gen)*N);
            if (ind >= N) ind = N - 1;
            memcpy(&C(0, 0), &X(0, ind), sizeof(T)*P);

            ArrayType aC;
            aC.create(P, C.begin());
            norm_C[0] = Gadgetron::dot(aC, aC,false);

            for (i = 1; i < K; i++)
            {
                // compute distance to the nearest centroid for all data points
                this->compute_dist(X, IDX, C, D);

                // compute norm of distance vector
                this->compute_norm_dist(D, D_norm);

                // compute accumulated distrance
                cumsum_D_norm(0) = D_norm(0);
                for (t = 1; t < N; t++)
                {
                    cumsum_D_norm(t) = cumsum_D_norm(t - 1) + D_norm(t);
                }

                if (std::abs(cumsum_D_norm(N - 1)) < FLT_EPSILON)
                {
                    GERROR_STREAM("std::abs(cumsum_D_norm(N-1))<FLT_EPSILON ... ");
                    // set centroid from i to K

                    for (s = i; s < K; s++)
                    {
                        size_t ind = (size_t)(dis(gen)*N);
                        if (ind >= N) ind = N - 1;
                        memcpy(&C(0, s), &X(0, ind), sizeof(T)*P);
                    }
                    break;
                }

                // convert to probability
                for (t = 0; t < N; t++)
                {
                    cumsum_D_norm(t) /= cumsum_D_norm(N - 1);
                }

                T v = dis(gen);
                for (t = 0; t < N; t++)
                {
                    if (cumsum_D_norm(t)>=v) break;
                }

                memcpy(&C(0, i), &X(0, t), sizeof(T)*P);

                aC.create(P, &C(0, i));
                norm_C[i] = Gadgetron::dot(aC, aC,false);

                // update the IDX
                ArrayType curr_C;
                curr_C.create(P, i+1, C.begin());

                this->update_IDX(X, curr_C, norm_C, IDX);
            }

            memcpy(&C_for_initial(0, 0, n), C.begin(), sizeof(T)*P*K);
        }

        if (this->perform_timing_) gt_timer_.stop();
    }
    catch (...)
    {
        GERROR_STREAM("Exceptions happened in kmeans<T>::get_initial_guess_kmeansplusplus(...) ... ");
    }
}

template <typename T>
void kmeans<T>::run_replicates(const ArrayType& X, size_t K, const ArrayType& C_for_initial, ClusterType& IDX, ArrayType& C, VectorType& sumD_rep, T& sumD)
{
    try
    {
        Gadgetron::GadgetronTimer timer;
        timer.set_timing_in_destruction(false);

        if (this->perform_timing_) gt_timer_.start("run_replicates");

        size_t P = X.get_size(0);
        size_t N = X.get_size(1);

        GADGET_CHECK_THROW(N>K);
        GADGET_CHECK_THROW(C_for_initial.get_size(0) == P);
        GADGET_CHECK_THROW(C_for_initial.get_size(1) == K);

        size_t R = C_for_initial.get_size(2);

        std::vector<ArrayType> C_rep(R);
        std::vector<ClusterType> IDX_rep(R);

        sumD_rep.resize(R, 0);

        size_t r;
        for (r=0; r<R; r++)
        {
            std::stringstream outs;
            outs << "-----> Kmeans, replicate " << r << " out of " << R;

            if (this->verbose_)
            {
                GDEBUG_STREAM(outs.str());
            }

            ArrayType curr_C_initial;
            curr_C_initial.create(P, K, const_cast<T*>(&C_for_initial(0, 0, r)) );

            if (this->perform_timing_) timer.start(outs.str().c_str());
            this->run(X, K, curr_C_initial, IDX_rep[r], C_rep[r], sumD_rep[r]);
            if (this->perform_timing_) timer.stop();

            if(this->verbose_)
            {
                GDEBUG_STREAM("Kmeans, replicate " << r << " out of " << R << " - " << sumD_rep[r]);
            }
        }

        size_t best_r = 0;
        sumD = sumD_rep[0];
        for (r = 1; r < R; r++)
        {
            if(sumD>sumD_rep[r])
            {
                sumD = sumD_rep[r];
                best_r = r;
            }
        }

        IDX = IDX_rep[best_r];
        C = C_rep[best_r];

        if (this->perform_timing_) gt_timer_.stop();
    }
    catch(...)
    {
        GERROR_STREAM("Exceptions happened in kmeans<T>::run_replicates(...) ... ");
    }
}

template <typename T>
void kmeans<T>::run(const ArrayType& X, size_t K, const ArrayType& C_for_initial, ClusterType& IDX, ArrayType& C, T& sumD)
{
    try
    {
        size_t P = X.get_size(0);
        size_t N = X.get_size(1);

        GADGET_CHECK_THROW(N>K);
        GADGET_CHECK_THROW(C_for_initial.get_size(0) == P);
        GADGET_CHECK_THROW(C_for_initial.get_size(1) == K);

        this->replicates_ = C_for_initial.get_size(2);

        IDX.resize(N, 0);
        C.create(P, K);
        Gadgetron::clear(C);
        sumD  = 0;

        C = C_for_initial;

        VectorType norm_C(K, 0);

        size_t k, p;
        for (k=0; k<K; k++)
        {
            T v = 0;
            for (p=0; p<P; p++)
            {
                v += C(p, k)*C(p, k);
            }

            norm_C[k] = v;
        }

        // first round of clustering
        this->update_IDX(X, C, norm_C, IDX);

        ClusterType prev_IDX;
        ArrayType D, D_norm;

        size_t num_iter = 0;
        T prev_sumD = std::numeric_limits<T>::max();

        std::vector<size_t> cluster_size;

        while (num_iter<=this->max_iter_ &&  this->is_clustering_changed(prev_IDX, IDX))
        {
            prev_IDX = IDX;

            // update the centroid
            this->update_centroid(X, IDX, C, norm_C);
            // update clustering
            this->update_IDX(X, C, norm_C, IDX);

            this->compute_dist(X, IDX, C, D);
            this->compute_norm_dist(D, D_norm);

            // if there are clusters having no member, find a point furthest away from its own cluster centroid
            // replace the empty cluster centroid with this point
            std::vector<size_t> empty_clusters;
            this->has_empty_cluster(C, IDX, cluster_size, empty_clusters);

            if (!empty_clusters.empty())
            {
                if (this->verbose_)
                {
                    GDEBUG_STREAM("Kmeas iteration, found empty cluster : iter - " << num_iter);
                }

                size_t e;
                for (e = 0; e < empty_clusters.size(); e++)
                {
                    // find the most "lonely" point
                    size_t lonely(0);
                    this->find_lonely_point(IDX, D_norm, cluster_size, lonely);

                    // replace this empty centroid with the lonely point, update IDX
                    IDX[lonely] = empty_clusters[e];

                    // updates centroid
                    this->update_centroid(X, IDX, C, norm_C);

                    // update distances
                    this->compute_dist(X, IDX, C, D);
                    this->compute_norm_dist(D, D_norm);
                }
            }

            sumD = 0;
            for (p = 0; p<N; p++)
            {
                sumD += D_norm(p)*D_norm(p);
            }

            if (sumD>prev_sumD)
            {
                if (this->verbose_)
                {
                    GDEBUG_STREAM("Kmeas iteration terminated due to increasted total distance : iter " << num_iter << " - " << sumD << " > " << prev_sumD);
                }

                IDX = prev_IDX;

                break;
            }
            else
            {
                prev_sumD = sumD;
            }

            num_iter++;
        }

        if (this->verbose_)
        {
            GDEBUG_STREAM("Kmeas iteration stopped : iter " << num_iter << " - " << sumD);
        }

        if(this->perform_online_update_)
        {
            this->perform_online_update(X, IDX, C, sumD);
            if (this->verbose_)
            {
                GDEBUG_STREAM("Kmeas online update : " << sumD);
            }
        }
    }
    catch (...)
    {
        GERROR_STREAM("Exceptions happened in kmeans<T>::run(...) ... ");
    }
}

template <typename T>
void kmeans<T>::compute_dist(const ArrayType& X, const ClusterType& IDX, const ArrayType& C, ArrayType& D)
{
    try
    {
        size_t P = X.get_size(0);
        size_t N = X.get_size(1);

        size_t K = C.get_size(1);

        GADGET_CHECK_THROW(N== IDX.size());

        D.create(P, N);
        Gadgetron::clear(D);

        const T* pX = X.begin();
        const T* pC = C.begin();
        T* pD = D.begin();

        long long p, n;

#pragma omp parallel for default(none) private(n, p) shared(N, IDX, K, P, pD, pX, pC)
        for (n=0; n<N; n++)
        {
            size_t nC = IDX[n];

            if (nC >= K)
            {
                // GWARN_STREAM("nC >= K :" << nC << " - " << K << " for " << n);
                continue;
            }

            for (p=0; p<P; p++)
            {
                pD[p + n*P] = pX[p + n*P] - pC[p + nC*P];
            }
        }
    }
    catch (...)
    {
        GERROR_STREAM("Exceptions happened in kmeans<T>::compute_dist(...) ... ");
    }
}

template <typename T>
void kmeans<T>::compute_norm_dist(const ArrayType& D, ArrayType& D_norm)
{
    try
    {
        size_t P = D.get_size(0);
        size_t N = D.get_size(1);

        D_norm.create(N);
        Gadgetron::clear(D_norm);

        const T* pD = D.begin();

        long long n, p;

#pragma omp parallel for default(none) private(n, p) shared(N, P, pD, D_norm)
        for(n=0; n<N; n++)
        {
            T v = 0;
            for (p = 0; p<P; p++)
            {
                v += pD[p + n*P] * pD[p + n*P];
            }

            D_norm(n) = std::sqrt(v);
        }
    }
    catch (...)
    {
        GERROR_STREAM("Exceptions happened in kmeans<T>::compute_norm_dist(...) ... ");
    }
}

template <typename T>
void kmeans<T>::update_IDX(const ArrayType& X, const ArrayType& C, const VectorType& norm_C, ClusterType& IDX)
{
    try
    {
        size_t P = X.get_size(0);
        size_t N = X.get_size(1);

        size_t K = C.get_size(1);

        IDX.resize(N);

        ArrayType CX;
        Gadgetron::gemm(CX, C, true, X, false);

        size_t t, s;

        for (t = 0; t < N; t++)
        {
            for (s = 0; s < K; s++)
            {
                CX(s, t) = 2* CX(s, t) - norm_C[s];
            }

            T maxCX = CX(0, t);
            IDX[t] = 0;
            for (s = 1; s < K; s++)
            {
                if (CX(s, t) > maxCX)
                {
                    maxCX = CX(s, t);
                    IDX[t] = s;
                }
            }
        }
    }
    catch (...)
    {
        GERROR_STREAM("Exceptions happened in kmeans<T>::update_IDX(...) ... ");
    }
}

template <typename T>
void kmeans<T>::update_centroid(const ArrayType& X, const ClusterType& IDX, ArrayType& C, VectorType& norm_C)
{
    try
    {
        size_t P = X.get_size(0);
        size_t N = X.get_size(1);
        const T* pX = X.begin();

        size_t K = C.get_size(1);

        norm_C.resize(K);

        std::vector<size_t> num_in_C(K, 0);

        Gadgetron::clear(C);
        T* pC = C.begin();

        size_t n, p;
        for (n=0; n<N; n++)
        {
            size_t currK = IDX[n];

            if(currK<K)
            {
                for (p=0; p<P; p++)
                {
                    pC[p+currK*P] += pX[p+n*P];
                }

                num_in_C[currK]++;
            }
            else
            {
                GERROR_STREAM("kmeans, currC>=K, in update_centroid : " << n);
            }
        }

        for (n = 0; n < K; n++)
        {
            T v = 0;
            for (p = 0; p<P; p++)
            {
                if (num_in_C[n] > 0)
                {
                    pC[p + n*P] /= num_in_C[n];
                    v += pC[p + n*P] * pC[p + n*P];
                }
            }

            norm_C[n] = v;
        }
    }
    catch (...)
    {
        GERROR_STREAM("Exceptions happened in kmeans<T>::update_centroid(...) ... ");
    }
}

template <typename T>
void kmeans<T>::has_empty_cluster(const ArrayType& C, const ClusterType& IDX, std::vector<size_t>& cluster_size, std::vector<size_t>& empty_clusters)
{
    try
    {
        size_t K = C.get_size(1);
        size_t N = IDX.size();

        cluster_size.resize(K, 0);

        size_t n;
        for (n=0; n<N; n++)
        {
            if (IDX[n] < K)
            {
                cluster_size[IDX[n]]++;
            }
        }

        empty_clusters.clear();
        for (n = 0; n < K; n++)
        {
            if(cluster_size[n]==0)
            {
                empty_clusters.push_back(n);
            }
        }
    }
    catch (...)
    {
        GERROR_STREAM("Exceptions happened in kmeans<T>::has_empty_cluster(...) ... ");
    }
}

template <typename T>
void kmeans<T>::find_lonely_point(const ClusterType& IDX, const ArrayType& D_norm, const std::vector<size_t>& cluster_size, size_t& lonely)
{
    try
    {
        size_t N = IDX.size();

        T maxD (0);
        Gadgetron::maxAbsolute(D_norm, maxD, lonely);

        size_t cluster = IDX[lonely];

        size_t n, c;

        if(cluster_size[cluster]<2)
        {
            // if the picked cluster has only one member ...
            std::vector<size_t> big_clusters;
            for (c=0; c<cluster_size.size(); c++)
            {
                if(cluster_size[c]>2)
                {
                    big_clusters.push_back(c);
                }
            }

            if(big_clusters.empty())
            {
                GADGET_THROW("All clusters have only one member ... ");
            }
            else
            {
                size_t cluster_picked = big_clusters[big_clusters.size() / 2]; // pick one cluster

                // find the furthest point in the picked cluster
                lonely = 0;
                maxD = 0;
                for (n=0; n<N; n++)
                {
                    if(IDX[n]==cluster_picked)
                    {
                        if(D_norm(n)>maxD)
                        {
                            maxD = D_norm(n);
                            lonely = n;
                        }
                    }
                }
            }
        }
    }
    catch (...)
    {
        GERROR_STREAM("Exceptions happened in kmeans<T>::find_lonely_point(...) ... ");
    }
}

template <typename T>
bool kmeans<T>::is_clustering_changed(const ClusterType&prev_IDX, const ClusterType& IDX)
{
    try
    {
        if (prev_IDX.size() != IDX.size()) return true;

        size_t N = IDX.size();

        size_t n;
        for (n=0; n<N; n++)
        {
            if (prev_IDX[n] != IDX[n])
            {
                return true;
            }
        }

        return false;
    }
    catch (...)
    {
        GERROR_STREAM("Exceptions happened in kmeans<T>::is_clustering_changed(...) ... ");
        return false;
    }
}

template <typename T>
void kmeans<T>::perform_online_update(const ArrayType& X, ClusterType& IDX, ArrayType& C, T& sumD)
{
    try
    {
        size_t P = X.get_size(0);
        size_t N = X.get_size(1);
        const T* pX = X.begin();

        size_t K = C.get_size(1);

        ArrayType del_cost; // store the delat change of sum cost if reassign a point
        del_cost.create(N, K);

        // count the number of points in each cluster
        std::vector<size_t> num_pt_clusters(K, 0);

        size_t n, p, k;

        num_pt_clusters.resize(K, 0);
        for (n = 0; n < N; n++)
        {
            num_pt_clusters[IDX[n]]++;
        }

        size_t iter(0);

        size_t lastmoved = 0;
        size_t nummoved = 0;
        ClusterType prevIDX, newIDX(IDX);

        while (iter < this->max_iter_)
        {
            // for every cluster K and every point N
            // compute change of delta sum cost
            for (k = 0; k < K; k++)
            {
                for (n = 0; n < N; n++)
                {
                    T v;
                    if (IDX[n] == k)
                    {
                        if (num_pt_clusters[k] > 1)
                            v = (T)num_pt_clusters[k] / (T)(num_pt_clusters[k] - 1);
                        else
                            v = 1;
                    }
                    else
                        v = (T)num_pt_clusters[k] / (T)(num_pt_clusters[k] + 1);

                    T t(0), d = 0;
                    for (p = 0; p < P; p++)
                    {
                        t = X(p, n) - C(p, k);
                        d += t*t;
                    }

                    del_cost(n, k) = v * d;
                }
            }

            prevIDX = IDX;

            // get the new IDX
            for (n = 0; n < N; n++)
            {
                newIDX[n] = 0;
                T min_del_cost = del_cost(n, 0);
                for (k = 1; k < K; k++)
                {
                    if(del_cost(n, k) < min_del_cost)
                    {
                        newIDX[n] = k;
                        min_del_cost = del_cost(n, k);
                    }
                }
            }

            // marked the moving points
            std::vector<size_t> moved;
            for (n = 0; n < N; n++)
            {
                if(prevIDX[n] != newIDX[n])
                {
                    moved.push_back(n);
                }
            }

            // if no candidates to move, stop
            if(moved.empty())
            {
                iter++;
                break;
            }

            // pick a point to move
            int moved_ind = N+1;
            int tt(0);
            for (size_t ii = 0; ii < moved.size(); ii++)
            {
                tt = (int)moved[ii] - (int)lastmoved - 1;
                if (tt < 0) tt += N;
                if (tt >= N) tt -= N;

                tt += lastmoved;

                if (tt < moved_ind) moved_ind = tt;
            }

            if (tt < 0) tt += N;
            if (tt >= N) tt -= N;
            moved_ind = (size_t)(tt);

            if(moved_ind<=lastmoved)
            {
                iter++;
                if (iter >= this->max_iter_) break;
                nummoved = 0;
            }

            nummoved++;
            lastmoved = moved_ind;

            size_t oidx = IDX[moved_ind];
            size_t nidx = newIDX[moved_ind];

            sumD = sumD + del_cost(moved_ind, nidx) - del_cost(moved_ind, oidx);

            IDX[moved_ind] = nidx;

            num_pt_clusters[oidx]--;
            num_pt_clusters[nidx]++;

            for (p=0; p<P; p++)
            {
                C(p, nidx) = C(p, nidx) + (X(p, moved_ind) - C(p, nidx)) / num_pt_clusters[nidx];
                C(p, oidx) = C(p, oidx) - (X(p, moved_ind) - C(p, oidx)) / num_pt_clusters[oidx];
            }
        }
    }
    catch (...)
    {
        GERROR_STREAM("Exceptions happened in kmeans<T>::perform_online_update(...) ... ");
    }
}

// ------------------------------------------------------------
// Instantiation
// ------------------------------------------------------------

template class EXPORTPR kmeans< float >;
template class EXPORTPR kmeans< double >;

}
