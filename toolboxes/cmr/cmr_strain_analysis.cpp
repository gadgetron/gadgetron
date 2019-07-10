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
    void compute_strain(const hoNDArray<T>& dx, const hoNDArray<T>& dy, const hoNDArray<T>& mask, hoNDArray<T>& radial, hoNDArray<T>& circ)
    {
        try
        {
            size_t RO = dx.get_size(0);
            size_t E1 = dx.get_size(1);
            size_t N = dx.get_size(2);

            typedef hoNDArray<T> ArrayType;

            radial.create(RO, E1, N);
            Gadgetron::clear(radial);

            circ.create(RO, E1, N);
            Gadgetron::clear(circ);

            size_t centroidX = 0;
            size_t centroidY = 0;
            size_t counter = 0;
            // find centroid
            for (size_t e1 = 0; e1 < E1; e1++)
            {

                for (size_t ro = 0; ro < RO; ro++)
                {
                    if (mask(ro, e1) > 0)
                    {
                        centroidX += ro;
                        centroidY += e1;
                        counter += 1;
                    }
                }
            }

            double Cx = (double)centroidX / counter;
            double Cy = (double)centroidY / counter;

            int phs;

            for (phs = 0; phs < N; phs++)
            {
                ArrayType dx_2D(RO, E1, const_cast<T*>(&dx(0, 0, phs)));
                ArrayType dy_2D(RO, E1, const_cast<T*>(&dy(0, 0, phs)));

                hoNDBoundaryHandlerBorderValue<ArrayType> bv_dx(dx_2D);
                hoNDInterpolatorLinear<ArrayType> interp_dx(dx_2D, bv_dx);

                hoNDBoundaryHandlerBorderValue<ArrayType> bv_dy(dy_2D);
                hoNDInterpolatorLinear<ArrayType> interp_dy(dy_2D, bv_dy);

                // compute strain
                for(size_t e1=0; e1<E1; e1++)
                {
                    for (size_t ro = 0; ro < RO; ro++)
                    {
                        // compute strain for this point
                        double x = e1 - (double)Cx;
                        double y = (double)Cy - ro;

                        double theta = atan(y / (x+FLT_EPSILON) ) + M_PI * (x < 0) + M_PI * 2 * (x >= 0)*(y < 0);

                        double x_in = e1 + 0.5 * cos(theta);
                        double x_out = e1 - 0.5 * cos(theta);
                        double y_in = ro - 0.5 * sin(theta);
                        double y_out = ro + 0.5 * sin(theta);

                        double distances = 1;

                        double dx_in = interp_dy(x_in, y_in);
                        double dx_out = interp_dy(x_out, y_out);
                        double dy_in = interp_dx(x_in, y_in);
                        double dy_out = interp_dx(x_out, y_out);

                        double x_prime_in = dx_in + x_in;
                        double x_prime_out = dx_out + x_out;
                        double y_prime_in = dy_in + y_in;
                        double y_prime_out = dy_out + y_out;

                        double a_x = x_prime_out - x_prime_in;
                        double a_y = y_prime_out - y_prime_in;
                        double b_x = x_out - x_in;
                        double b_y = y_out - y_in;
                        double comp_ab_rad = a_x*b_x + a_y*b_y;
                        radial(ro, e1, phs) = (comp_ab_rad - distances) / distances;

                        // double b_x_rot = b_x*cos(M_PI / 2) + b_y*sin(M_PI / 2);
                        double b_x_rot = b_y;
                        // double b_y_rot = -b_x*sin(M_PI / 2) + b_y*cos(M_PI / 2);
                        double b_y_rot = -b_x;
                        double comp_ab_circ = a_x*b_x_rot + a_y*b_y_rot;
                        circ(ro, e1, phs) = (comp_ab_circ - distances) / distances;
                    }
                }
            }
        }
        catch(...)
        {
            GADGET_THROW("Exceptions happened in compute_strain(...) ... ");
        }
    }

    template EXPORTCMR void compute_strain(const hoNDArray<float>& dx, const hoNDArray<float>& dy, const hoNDArray<float>& mask, hoNDArray<float>& radial, hoNDArray<float>& circ);
    template EXPORTCMR void compute_strain(const hoNDArray<double>& dx, const hoNDArray<double>& dy, const hoNDArray<double>& mask, hoNDArray<double>& radial, hoNDArray<double>& circ);
}
