/** \file   cmr_stain_analysis.h
    \brief  Implement functionalities to handle cardiac time stamps
    \author Angela Gao
    \date   July 9, 2019
*/

#include "cmr_strain_analysis.h"
#include "hoNDArray_reductions.h"
#include "hoNDArray_elemwise.h"
#include "hoNDInterpolator.h"
#include "hoNDBoundaryHandler.h"
#define _USE_MATH_DEFINES
#include <math.h>

namespace Gadgetron {

    template <typename T>
    void compute_strain(const hoNDArray<double>& dx, const hoNDArray<double>& dy, const hoNDArray<T>& mask, const bool compare_mask, hoNDArray<T>& radial, hoNDArray<T>& circ, hoNDArray<T>& thetas)
    {
        try
        {
            size_t RO = dx.get_size(0);
            size_t E1 = dx.get_size(1);
            size_t N = dx.get_size(2);

            typedef hoNDArray<double> ArrayType;

            ArrayType& dx_used = const_cast<ArrayType&>(dx);

            hoNDBoundaryHandlerBorderValue<ArrayType> bhBorderValue(dx_used);
            hoNDInterpolatorLinear<ArrayType> interpLinear(dx_used, bhBorderValue);

            radial.create(RO, E1, N);
            Gadgetron::clear(radial);

            circ.create(RO, E1, N);
            Gadgetron::clear(circ);

            thetas.create(RO, E1, N);
            Gadgetron::clear(thetas);

            size_t centroidR = 0;
            size_t centroidE = 0;
            size_t counter = 0;
            // find centroid
            for (size_t e1 = 0; e1 < E1; e1++)
            {

                for (size_t ro = 0; ro < RO; ro++)
                {
                    if (mask(ro, e1) > 0)
                    {
                        centroidR += ro;
                        centroidE += e1;
                        counter += 1;
                    }
                }
            }

            double Cr = (double)centroidR / counter;
            double Ce = (double)centroidE / counter;

            int phs;

#pragma omp parallel for private(phs) shared(N, RO, E1, dx, dy, radial, circ, mask, Cr, Ce, thetas)
            for (phs = 0; phs < N; phs++)
            {
                ArrayType dx_2D(RO, E1, const_cast<double*>(&dx(0, 0, phs)));
                ArrayType dy_2D(RO, E1, const_cast<double*>(&dy(0, 0, phs)));

                hoNDBoundaryHandlerBorderValue<ArrayType> bv_dx(dx_2D);
                hoNDInterpolatorLinear<ArrayType> interp_dx(dx_2D, bv_dx);

                hoNDBoundaryHandlerBorderValue<ArrayType> bv_dy(dy_2D);
                hoNDInterpolatorLinear<ArrayType> interp_dy(dy_2D, bv_dy);

                // compute strain
                for (size_t e1 = 0; e1 < E1; e1++)
                {
                    for (size_t ro = 0; ro < RO; ro++)
                    {
                        // compute strain for this point
                        double x = ro - (double)Cr;
                        double y = (double)Ce - e1;

                        double theta = atan(y / (x + FLT_EPSILON)) + M_PI * (x < 0) + M_PI * 2 * (x >= 0) * (y < 0);

                        double x_in = ro + 0.5 * cos(theta);
                        double x_out = ro - 0.5 * cos(theta);
                        double y_in = e1 - 0.5 * sin(theta);
                        double y_out = e1 + 0.5 * sin(theta);

                        double distances = 1;

                        double dx_in = interp_dx(x_in, y_in);
                        double dx_out = interp_dx(x_out, y_out);
                        double dy_in = interp_dy(x_in, y_in);
                        double dy_out = interp_dy(x_out, y_out);

                        thetas(ro, e1) = theta;

                        double x_prime_in = dx_in + x_in;
                        double x_prime_out = dx_out + x_out;
                        double y_prime_in = dy_in + y_in;
                        double y_prime_out = dy_out + y_out;

                        double a_x = x_prime_out - x_prime_in;
                        double a_y = y_prime_out - y_prime_in;
                        double b_x = x_out - x_in;
                        double b_y = y_out - y_in;
                        double comp_ab_rad = a_x * b_x + a_y * b_y;

                        double theta_rot = theta + M_PI / 2;

                        double x_in_rot = ro + 0.5 * cos(theta_rot);
                        double x_out_rot = ro - 0.5 * cos(theta_rot);
                        double y_in_rot = e1 - 0.5 * sin(theta_rot);
                        double y_out_rot = e1 + 0.5 * sin(theta_rot);

                        double dx_in_rot = interp_dx(x_in_rot, y_in_rot);
                        double dx_out_rot = interp_dx(x_out_rot, y_out_rot);
                        double dy_in_rot = interp_dy(x_in_rot, y_in_rot);
                        double dy_out_rot = interp_dy(x_out_rot, y_out_rot);

                        double x_prime_in_rot = dx_in_rot + x_in_rot;
                        double x_prime_out_rot = dx_out_rot + x_out_rot;
                        double y_prime_in_rot = dy_in_rot + y_in_rot;
                        double y_prime_out_rot = dy_out_rot + y_out_rot;

                        double a_x_rot = x_prime_out_rot - x_prime_in_rot;
                        double a_y_rot = y_prime_out_rot - y_prime_in_rot;
                        double b_x_rot = x_out_rot - x_in_rot;
                        double b_y_rot = y_out_rot - y_in_rot;
                        double comp_ab_circ = a_x_rot * b_x_rot + a_y_rot * b_y_rot;

                        if (compare_mask == true)
                        {
                            radial(ro, e1, phs) = (comp_ab_rad - distances) / distances * mask(ro, e1);
                            circ(ro, e1, phs) = (comp_ab_circ - distances) / distances * mask(ro, e1);
                        }
                        else
                        {
                            radial(ro, e1, phs) = (comp_ab_rad - distances) / distances;
                            circ(ro, e1, phs) = (comp_ab_circ - distances) / distances;
                        }
                    }
                }
            }
        }
        catch (...)
        {
            GADGET_THROW("Exceptions happened in compute_strain(...) ... ");
        }
    }

    template EXPORTCMR void compute_strain(const hoNDArray<double>& dx, const hoNDArray<double>& dy, const hoNDArray<float>& mask, const bool compare_mask, hoNDArray<float>& radial, hoNDArray<float>& circ, hoNDArray<float>& thetas);
    template EXPORTCMR void compute_strain(const hoNDArray<double>& dx, const hoNDArray<double>& dy, const hoNDArray<double>& mask, const bool compare_mask, hoNDArray<double>& radial, hoNDArray<double>& circ, hoNDArray<double>& thetas);
}
