/** \file   cmr_spirit_recon.h
    \brief  Implement some functionalities commonly used in cmr applications for spirit recon
    \author Hui Xue
*/

#pragma once

#include "hoNDArray.h"

namespace Gadgetron {

    /// perform spirit linear 2DT recon
    /// kspace: input kspace [RO E1 CHA N S]
    /// startE1, endE1 : sampling range along E1
    /// kerIm: kspace kernel image domain [RO E1 CHA CHA 1orN 1orS]
    /// kspaceInitial: initial guess for the solution; if empty, zero initialization will be used
    /// res: recon results
    /// iter_max: maximal number of iterations
    /// iter_thres: iteration stop threshold
    /// print_iter: whether to print iteration information
    template <typename T> void perform_spirit_recon_linear_2DT(const Gadgetron::hoNDArray<T>& kspace, size_t startE1, size_t endE1, const Gadgetron::hoNDArray<T>& kerIm, const Gadgetron::hoNDArray<T>& kspaceInitial, Gadgetron::hoNDArray<T>& res, size_t iter_max=90, double iter_thres=0.0015, bool print_iter=false);

    /// perform spirit nonlinear 2DT recon
    /// kspace: [RO E1 CHA N S]
    /// kerIm: kspace kernel image domain [RO E1 CHA CHA 1orN 1orS]
    /// coilMap: [RO E1 CHA 1orN 1orS]
    template <typename T> void perform_spirit_recon_non_linear_2DT(const Gadgetron::hoNDArray<T>& kspace, const Gadgetron::hoNDArray<T>& kerIm, const Gadgetron::hoNDArray<T>& coilMap, const Gadgetron::hoNDArray<T>& kspaceInitial, Gadgetron::hoNDArray<T>& res,
        size_t iter_max=15, double iter_thres=0.004, double data_fidelity_lamda=1.0, double image_reg_lamda=0.0005, double reg_N_weighting_ratio=10.0, bool reg_use_coil_sen_map=false, bool reg_with_approx_coeff=false, const std::string& wav_name="db1", bool print_iter=false);
}
