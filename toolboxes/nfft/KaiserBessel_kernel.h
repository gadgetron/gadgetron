//
// Kaiser-Bessel convolution kernels
//
#pragma once

#include "core_defines.h"


namespace Gadgetron {

    __inline__ __device__ __host__

    double
    bessi0(double x) {
        double denominator;
        double numerator;
        double z;

        if (x == 0.0) {
            return 1.0;
        } else {
            z = x * x;
            numerator = (z * (z * (z * (z * (z * (z * (z * (z * (z * (z * (z * (z * (z *
                                                                                     (z * 0.210580722890567e-22 +
                                                                                      0.380715242345326e-19) +
                                                                                     0.479440257548300e-16) +
                                                                                0.435125971262668e-13) +
                                                                           0.300931127112960e-10) +
                                                                      0.160224679395361e-7) +
                                                                 0.654858370096785e-5) + 0.202591084143397e-2) +
                                                       0.463076284721000e0) + 0.754337328948189e2) +
                                             0.830792541809429e4) + 0.571661130563785e6) +
                                   0.216415572361227e8) + 0.356644482244025e9) +
                         0.144048298227235e10);

            denominator = (z * (z * (z - 0.307646912682801e4) +
                                0.347626332405882e7) - 0.144048298227235e10);
        }

        return -numerator / denominator;
    }

    __inline__ __device__ __host__

    float
    bessi0(float x) {
        float denominator;
        float numerator;
        float z;

        if (x == 0.0f) {
            return 1.0f;
        } else {
            z = x * x;
            numerator = (z * (z * (z * (z * (z * (z * (z * (z * (z * (z * (z * (z * (z *
                                                                                     (z * 0.210580722890567e-22f +
                                                                                      0.380715242345326e-19f) +
                                                                                     0.479440257548300e-16f) +
                                                                                0.435125971262668e-13f) +
                                                                           0.300931127112960e-10f) +
                                                                      0.160224679395361e-7f) +
                                                                 0.654858370096785e-5f) + 0.202591084143397e-2f) +
                                                       0.463076284721000e0f) + 0.754337328948189e2f) +
                                             0.830792541809429e4f) + 0.571661130563785e6f) +
                                   0.216415572361227e8f) + 0.356644482244025e9f) +
                         0.144048298227235e10f);

            denominator = (z * (z * (z - 0.307646912682801e4f) +
                                0.347626332405882e7f) - 0.144048298227235e10f);
        }

        return -numerator / denominator;
    }


// Kaiser Bessel according to Beatty et. al. IEEE TMI 2005;24(6):799-808.
// There is a slight difference wrt Jackson's formulation, IEEE TMI 1991;10(3):473-478.

    __inline__ __device__ __host__

    double
    KaiserBessel(double u, double matrix_size_os, double one_over_W, double beta) {
        double _tmp = 2.0 * u * one_over_W;
        double tmp = _tmp * _tmp;
        double arg = beta * std::sqrt(1.0 - tmp);
        double bessi = bessi0(arg);
        double ret = matrix_size_os * bessi * one_over_W;
        return ret;
    }

    __inline__ __device__ __host__

    float
    KaiserBessel(float u, float matrix_size_os, float one_over_W, float beta) {
        float _tmp = 2.0f * u * one_over_W;
        float tmp = _tmp * _tmp;
        float arg = beta * std::sqrt(1.0f - tmp);
        float bessi = bessi0(arg);
        float ret = matrix_size_os * bessi * one_over_W;
        return ret;
    }

//
// Below the intended interface
//

    template<class REAL> __inline__ __device__ __host__

    REAL
    KaiserBessel(const Gadgetron::vector_td<REAL, 1> &u, const Gadgetron::vector_td<REAL, 1> &matrix_size_os,
                 REAL one_over_W, const vector_td<REAL, 1> &beta) {
        REAL phi_x = KaiserBessel(u.vec[0], matrix_size_os.vec[0], one_over_W, beta[0]);
        return phi_x;
    }

    template<class REAL> __inline__ __device__ __host__

    REAL
    KaiserBessel(const Gadgetron::vector_td<REAL, 2> &u, const Gadgetron::vector_td<REAL, 2> &matrix_size_os,
                 REAL one_over_W, const vector_td<REAL, 2> &beta) {

        REAL phi_x = KaiserBessel(u.vec[0], matrix_size_os.vec[0], one_over_W, beta[0]);
        REAL phi_y = KaiserBessel(u.vec[1], matrix_size_os.vec[1], one_over_W, beta[1]);
        return phi_x * phi_y;
    }

    template<class REAL> __inline__ __device__ __host__

    REAL
    KaiserBessel(const Gadgetron::vector_td<REAL, 3> &u, const Gadgetron::vector_td<REAL, 3> &matrix_size_os,
                 REAL one_over_W, const vector_td<REAL, 3> &beta) {
        REAL phi_x = KaiserBessel(u.vec[0], matrix_size_os.vec[0], one_over_W, beta[0]);
        REAL phi_y = KaiserBessel(u.vec[1], matrix_size_os.vec[1], one_over_W, beta[1]);
        REAL phi_z = KaiserBessel(u.vec[2], matrix_size_os.vec[2], one_over_W, beta[2]);
        return phi_x * phi_y * phi_z;
    }

    template<class REAL> __inline__ __device__ __host__

    REAL
    KaiserBessel(const Gadgetron::vector_td<REAL, 4> &u, const Gadgetron::vector_td<REAL, 4> &matrix_size_os,
                 REAL one_over_W, const vector_td<REAL, 4> &beta) {
        REAL phi_x = KaiserBessel(u.vec[0], matrix_size_os.vec[0], one_over_W, beta[0]);
        REAL phi_y = KaiserBessel(u.vec[1], matrix_size_os.vec[1], one_over_W, beta[1]);
        REAL phi_z = KaiserBessel(u.vec[2], matrix_size_os.vec[2], one_over_W, beta[2]);
        REAL phi_w = KaiserBessel(u.vec[3], matrix_size_os.vec[3], one_over_W, beta[3]);
        return phi_x * phi_y * phi_z * phi_w;
    }

}