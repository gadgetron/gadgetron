/** \file   cmr_radial_thickening.h
    \brief  Implement functionalities to handle cardiac time stamps
    \author Angela Gao
    \date   July 9, 2019
*/

#include "cmr_analytical_strain.h"
#include "hoNDArray_reductions.h"
#include "hoNDArray_elemwise.h"
#include "hoNDInterpolator.h"
#include "hoNDBoundaryHandler.h"
#define _USE_MATH_DEFINES
#include <math.h>

namespace Gadgetron {

    template <typename T>
    void compute_analytical_strain(const hoNDArray<double>& dx, const hoNDArray<double>& dy, const hoNDArray<T>& mask, hoNDArray<T>& rad_strain, hoNDArray<T>& circ_strain)
    {
        try
        {
            size_t RO = mask.get_size(0);
            size_t E1 = mask.get_size(1);
            size_t PHS = mask.get_size(2);

            typedef hoNDArray<T> ArrayType;

            hoNDArray<T> centroidE(PHS), centroidR(PHS);
            Gadgetron::clear(centroidE);
            Gadgetron::clear(centroidR);

            hoNDArray<T> counter(PHS);
            Gadgetron::clear(counter);

            rad_strain.create(RO, E1, PHS);
            Gadgetron::clear(rad_strain);

            circ_strain.create(RO, E1, PHS);
            Gadgetron::clear(circ_strain);

            // find centroid
            for (size_t e1 = 0; e1 < E1; e1++)
            {

                for (size_t ro = 0; ro < RO; ro++)
                {
                    for (size_t p = 0; p < PHS; p++)
                    {
                        if (mask(ro, e1, p) > 0)
                        {
                            centroidR(p) += ro;
                            centroidE(p) += e1;
                            counter(p) += 1;
                        }
                    }
                }
            }

            Gadgetron::divide(centroidR, counter, centroidR);
            Gadgetron::divide(centroidE, counter, centroidE);

            size_t p, ro, e1;
// #pragma omp parallel for default(none) private(p, e1, ro) shared(PHS, RO, E1, rad_strain, circ_strain, centroidE, centroidR, mask, dx, dy)
            for (p = 0; p < PHS; p++)
            {
                for (e1 = 1; e1 < (E1-1); e1++)
                {
                    for (ro = 1; ro < (RO-1); ro++)
                    {
                        if (mask(ro, e1, p) == 1)
                        {
                            // compute strain for this point
                            double x = ro - (double)centroidR(p);
                            double y = (double)centroidE(p) - e1;

                            double theta = atan(y / (x + FLT_EPSILON)) + M_PI * (x < 0) + M_PI * 2 * (x >= 0) * (y < 0);

                            double ddrdr = (dx(ro + 1, e1, p) - dx(ro - 1, e1, p)) / 2;
                            double ddrde = -(dx(ro, e1 + 1, p) - dx(ro, e1 - 1, p)) / 2;
                            double ddedr = -(dy(ro + 1, e1, p) - dy(ro - 1, e1, p)) / 2;
                            double ddede = (dy(ro, e1 + 1, p) - dy(ro, e1 - 1, p)) / 2;

                            double f00 = (1 + ddrdr);
                            double f01 = ddrde;
                            double f10 = ddedr;
                            double f11 = (1 + ddede);

                            double e00 = 1.0 / 2.0 * (f00 * f00 + f10 * f10 - 1);
                            double e01 = 1.0 / 2.0 * (f00 * f01 + f10 * f11);
                            double e10 = 1.0 / 2.0 * (f00 * f01 + f10 * f11);
                            double e11 = 1.0 / 2.0 * (f01 * f01 + f11 * f11 - 1);

                            rad_strain(ro, e1, p) = e00 * (std::pow(cos(theta), 2)) + e11 * (std::pow(sin(theta), 2)) + 2 * e01 * sin(theta) * cos(theta);
                            circ_strain(ro, e1, p) = e00 * (std::pow(sin(theta), 2)) + e11 * (std::pow(cos(theta), 2)) - 2 * e01 * sin(theta) * cos(theta);
                        }
                    }
                }
            }
        }
        catch (...)
        {
            GADGET_THROW("Exceptions happened in compute_analytical_strain(...) ... ");
        }
    }

    template void compute_analytical_strain(const hoNDArray<double>& dx, const hoNDArray<double>& dy, const hoNDArray<float>& mask, hoNDArray<float>& rad_strain, hoNDArray<float>& circ_strain);
    template void compute_analytical_strain(const hoNDArray<double>& dx, const hoNDArray<double>& dy, const hoNDArray<double>& mask, hoNDArray<double>& rad_strain, hoNDArray<double>& circ_strain);
}