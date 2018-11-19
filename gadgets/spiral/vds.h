#pragma once

#include "gadgetron_spiral_export.h"
#include "hoNDArray.h"
#include "vector_td.h"

namespace Gadgetron {

    void EXPORTGADGETS_SPIRAL
    calc_vds(double slewmax, double gradmax, double Tgsample, double Tdsample, int Ninterleaves,
             double *fov, int numfov, double krmax,
             int ngmax, double **xgrad, double **ygrad, int *numgrad);

    void EXPORTGADGETS_SPIRAL
    calc_traj(double *xgrad, double *ygrad, int ngrad, int Nints, double Tgsamp, double krmax,
              double **x_trajectory, double **y_trajectory, double **weights);


    EXPORTGADGETS_SPIRAL
    hoNDArray<floatd2> calculate_trajectories(const hoNDArray<floatd2> &gradients, float sample_time, float krmax);

    EXPORTGADGETS_SPIRAL
    hoNDArray<float> calculate_weights(const hoNDArray<floatd2> &gradients, const hoNDArray<floatd2> &trajectories);


    /**
     * Calculates the density compensation weights according to Hoge, R. D., Kwan, R. K. and Bruce Pike, G. (1997),
     * Density compensation functions for spiral MRI. Magn. Reson. Med., 38: 117-128. doi:10.1002/mrm.1910380117
     * @param gradients
     * @param trajectories
     * @return
     */
    EXPORTGADGETS_SPIRAL
    hoNDArray<float>
    calculate_weights_Hoge(const hoNDArray<floatd2> &gradients, const hoNDArray<floatd2> &trajectories);

    EXPORTGADGETS_SPIRAL
    hoNDArray<floatd2> calculate_vds(double slewmax, double gradmax, double Tgsample, double Tdsample, int Ninterleaves,
                                     double *fov, int numfov, double krmax, int ngmax, int max_nsamples);

    EXPORTGADGETS_SPIRAL
    hoNDArray<floatd2> create_rotations(const hoNDArray<floatd2> &trajectories, int Nints);

}
