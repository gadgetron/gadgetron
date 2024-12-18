
/** \file   mri_core_grappa_python.h
    \brief  Python binding for mri_core grappa functinalities
    \author Hui Xue
*/

#include <complex>
#include <vector>
#include "hoNDArray.h"

namespace Gadgetron
{
    class grappa2D
    {
    public:

        typedef std::complex<float> T;

        grappa2D();
        ~grappa2D();

        int initialize(size_t accelFactor, size_t kRO, size_t kNE1, bool fitItself, double thres);
        int calib(hoNDArray<T> acsSrc, hoNDArray<T> acsDst);

        // kspace [RO, E1, srcCHA], perform recon using ker_
        hoNDArray<T> recon(hoNDArray<T> kspace, bool periodic_boundary_condition);

        // if kspace is already converted to data matrix (dataA and dataAInd), perform recon using ker_
        hoNDArray<T> recon_data_matrix(hoNDArray<T> dataA, hoNDArray<unsigned short> dataAInd, size_t RO, size_t E1, bool periodic_boundary_condition);
        // if recon = dataA*ker_, fill back recon points to result kspace
        hoNDArray<T> recon_fill_back_kspace(hoNDArray<T> recon, hoNDArray<unsigned short> dataAInd, size_t RO, size_t E1);

        hoNDArray<T> get_A() 
        {
            return A_;
        }

        hoNDArray<T> get_B()
        {
            return B_;
        }

        hoNDArray<T> get_ker()
        {
            return ker_;
        }

        hoNDArray<T> get_data_A()
        {
            return data_A_;
        }

        hoNDArray<unsigned short> get_data_A_index()
        {
            return data_AInd_;
        }

        void help();
        void status();

    protected:

        hoNDArray<T> A_;
        hoNDArray<T> B_;
        hoNDArray<T> ker_;

        hoNDArray<T> data_A_;
        hoNDArray<unsigned short> data_AInd_;

        size_t kRO_;
        size_t accelFactor_;
        std::vector<int> kE1_;
        std::vector<int> oE1_;
        size_t convKRO_;
        size_t convKE1_;
        double thres_;
    };
}
