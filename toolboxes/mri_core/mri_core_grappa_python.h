
/** \file   mri_core_grappa_python.h
    \brief  Python binding for mri_core grappa functinalities
    \author Hui Xue
*/

#include "mri_core_export.h"
#include <complex>
#include <vector>
#include "hoNDArray.h"

namespace Gadgetron
{
    class EXPORTMRICORE grappa2D
    {
    public:

        typedef std::complex<float> T;

        grappa2D();
        ~grappa2D();

        int initialize(size_t accelFactor, size_t kRO, size_t kNE1, bool fitItself, double thres);
        int calib(hoNDArray<T> acsSrc, hoNDArray<T> acsDst);
        hoNDArray<T> recon(hoNDArray<T> kspace, bool periodic_boundary_condition);

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
