/** \file       pattern_recognition_test.cpp
    \brief      Test case for pattern recognition functionalities

    \author     Hui Xue
*/

#include "ImageIOAnalyze.h"
#include "GadgetronTimer.h"
#include "pr_kmeans.h"
#include <gtest/gtest.h>
#include <random>

using namespace Gadgetron;
using testing::Types;

template<typename T> class pattern_recognition_test : public ::testing::Test
{
protected:
    virtual void SetUp()
    {
    }
};

typedef Types<float> realImplementations;
TYPED_TEST_SUITE(pattern_recognition_test, realImplementations);

TYPED_TEST(pattern_recognition_test, kmeans_test)
{
    std::default_random_engine generator;
    std::normal_distribution<float> distribution(0.0f, 1.0f);

    Gadgetron::kmeans<float> km;

    km.max_iter_ = 100;
    km.replicates_ = 20;

    km.verbose_ = true;
    km.perform_timing_ = true;
    // km.debug_folder_ = "D:/gtuser/mrprogs/install/DebugOutputg/";

    // set up data array

    size_t P = 2;
    size_t N = 4000;

    hoNDArray<float> X;
    X.create(P, N);

    size_t n;

    for (n=0; n<N; n++)
    {
        if (n < N/2)
        {
            X(0, n) = distribution(generator) + 4;
            X(1, n) = distribution(generator) + 4;
        }
        else
        {
            X(0, n) = distribution(generator) - 4;
            X(1, n) = distribution(generator) - 4;
        }
    }

    //ImageIOAnalyze gt_io;
    //gt_io.export_array(X, km.debug_folder_ + "X");

    size_t K = 2;

    size_t k, p, r;

    hoNDArray<float> C_for_initial;
    km.get_initial_guess_sample(X, K, C_for_initial);
    std::cout << " ----------------------------------------" << std::endl;
    std::cout << "get_initial_guess_sample ... " << std::endl;

    for (r = 0; r < km.replicates_; r++)
    {
        std::cout << r << " - ";
        for (k = 0; k < K; k++)
        {
            std::cout << k << " - ";
            for (p = 0; p < P; p++)
            {
                std::cout << C_for_initial(p, k, r) << " , ";
            }
        }
        std::cout << std::endl;
    }

    std::cout << " ----------------------------------------" << std::endl;
    km.get_initial_guess_uniform(X, K, C_for_initial);
    std::cout << "get_initial_guess_uniform ... " << std::endl;
    for (r = 0; r < km.replicates_; r++)
    {
        std::cout << r << " - ";
        for (k = 0; k < K; k++)
        {
            std::cout << k << " - ";
            for (p = 0; p < P; p++)
            {
                std::cout << C_for_initial(p, k, r) << " , ";
            }
        }
        std::cout << std::endl;
    }
    std::cout << " ----------------------------------------" << std::endl;
    km.get_initial_guess_kmeansplusplus(X, K, C_for_initial);
    std::cout << "get_initial_guess_kmeansplusplus ... " << std::endl;
    for (r = 0; r < km.replicates_; r++)
    {
        std::cout << r << " - ";
        for (k = 0; k < K; k++)
        {
            std::cout << k << " - ";
            for (p = 0; p < P; p++)
            {
                std::cout << C_for_initial(p, k, r) << " , ";
            }
        }
        std::cout << std::endl;
    }
    std::cout << " ----------------------------------------" << std::endl;

    hoNDArray<float> C;
    C.create(2, K);
    C(0, 0) = 4;
    C(1, 0) = 4;

    C(0, 1) = -4;
    C(1, 1) = -4;

    std::vector<size_t> IDX;
    hoNDArray<float> C_res;
    float sumD;
    km.run(X, K, C, IDX, C_res, sumD);

    std::cout << " ----------------------------------------" << std::endl;
    std::cout << "kmeans results with known centroid " << std::endl;
    for (k = 0; k < K; k++)
    {
        std::cout << k << " - ";
        for (p = 0; p < P; p++)
        {
            std::cout << C_res(p, k) << " , ";
        }
    }
    std::cout << std::endl;
    std::cout << " ----------------------------------------" << std::endl;

    EXPECT_LE( std::sqrt(sumD) / N, 2.0);

    std::cout << " ----------------------------------------" << std::endl;
    std::cout << "kmeans with known centroid " << std::endl;
    for (k = 0; k < K; k++)
    {
        std::cout << k << " - ";
        for (p = 0; p < P; p++)
        {
            std::cout << C_res(p, k) << " , ";
        }
    }
    std::cout << std::endl;
    std::cout << " ----------------------------------------" << std::endl;

    std::vector<float> sumD_rep;
    km.run_replicates(X, K, C_for_initial, IDX, C_res, sumD_rep, sumD);

    std::cout << " ----------------------------------------" << std::endl;
    std::cout << "kmeans with replicated centroids " << std::endl;
    for (k = 0; k < K; k++)
    {
        std::cout << k << " - ";
        for (p = 0; p < P; p++)
        {
            std::cout << C_res(p, k) << " , ";
        }
    }
    std::cout << std::endl;
    std::cout << " ----------------------------------------" << std::endl;

    //hoNDArray<size_t> IDX_array;
    //IDX_array.create(IDX.size(), &IDX[0]);
    //gt_io.export_array(IDX_array, km.debug_folder_ + "IDX");

    EXPECT_LE( std::sqrt(sumD) / N, 2.0);
}
