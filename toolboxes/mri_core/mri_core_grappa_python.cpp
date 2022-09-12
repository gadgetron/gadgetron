
/** \file   mri_core_grappa_python.cpp
    \brief  Python binding implementation for mri_core grappa
    \author Hui Xue
*/

#include "mri_core_grappa_python.h"
#include "mri_core_grappa.h"
#include "hoNDArray_utils.h"
#include "hoNDArray_elemwise.h"
#include "hoNDArray_reductions.h"
#include "python_toolbox.h"

#ifdef USE_OMP
    #include "omp.h"
#endif // USE_OMP

namespace Gadgetron
{
    grappa2D::grappa2D()
    {
        // ensure boost can convert between hoNDArrays and NumPy arrays automatically
        register_converter<hoNDArray<std::complex<float> > >();
        register_converter<hoNDArray< float > >();
        register_converter<hoNDArray< uint32_t > >();
        register_converter<hoNDArray< unsigned short > >();

        this->thres_ = 1e-4;
        this->kRO_ = 5;
        this->accelFactor_ = 0;
        this->convKRO_ = 0;
        this->convKE1_ = 0;
    }

    grappa2D::~grappa2D()
    {
    }

    int grappa2D::initialize(size_t accelFactor, size_t kRO, size_t kNE1, bool fitItself, double thres)
    {
        try
        {
            Gadgetron::grappa2d_kerPattern(this->kE1_, this->oE1_, this->convKRO_, this->convKE1_, accelFactor, kRO, kNE1, fitItself);
            this->kRO_ = kRO;
            this->accelFactor_ = accelFactor;
            this->thres_ = thres;
        }
        catch(...)
        {
            GERROR_STREAM("Exceptions happened in grappa2D::initialize(...) ... ");
            return 1;
        }

        return 0;
    }

    int grappa2D::calib(hoNDArray< std::complex<float> > acsSrc, hoNDArray< std::complex<float> > acsDst)
    {
        try
        {
            size_t startRO = 0;
            size_t endRO = acsSrc.get_size(0)-1;
            size_t startE1 = 0;
            size_t endE1 = acsSrc.get_size(1)-1;

            Gadgetron::grappa2d_prepare_calib(acsSrc, acsDst, this->kRO_, this->kE1_, this->oE1_, startRO, endRO, startE1, endE1, this->A_, this->B_);
            Gadgetron::grappa2d_perform_calib(this->A_, this->B_, this->kRO_, this->kE1_, this->oE1_, this->thres_, this->ker_);
        }
        catch(...)
        {
            GERROR_STREAM("Exceptions happened in grappa2D::calib(...) ... ");
            return 1;
        }

        return 0;
    }

    hoNDArray< std::complex<float> > grappa2D::recon(hoNDArray< std::complex<float> > kspace, bool periodic_boundary_condition)
    {
        hoNDArray<T> res;
        try
        {
            Gadgetron::grappa2d_prepare_recon(kspace, this->kRO_, this->kE1_, this->oE1_, periodic_boundary_condition, this->data_A_, this->data_AInd_);
            Gadgetron::grappa2d_perform_recon(this->data_A_, this->ker_, this->data_AInd_, this->oE1_, kspace.get_size(0), kspace.get_size(1), res);
        }
        catch (...)
        {
            GERROR_STREAM("Exceptions happened in grappa2D::recon(...) ... ");
            throw;
        }

        return res;
    }

    hoNDArray< std::complex<float> > grappa2D::recon_data_matrix(hoNDArray<T> dataA, hoNDArray<unsigned short> dataAInd, size_t RO, size_t E1, bool periodic_boundary_condition)
    {
        hoNDArray<T> res;
        try
        {
            Gadgetron::grappa2d_perform_recon(dataA, this->ker_, dataAInd, this->oE1_, RO, E1, res);
        }
        catch (...)
        {
            GERROR_STREAM("Exceptions happened in grappa2D::recon_data_matrix(...) ... ");
            throw;
        }

        return res;
    }

    hoNDArray< std::complex<float> > grappa2D::recon_fill_back_kspace(hoNDArray<T> recon, hoNDArray<unsigned short> dataAInd, size_t RO, size_t E1)
    {
        hoNDArray<T> res;
        try
        {
            Gadgetron::grappa2d_fill_reconed_kspace(dataAInd, recon, this->oE1_, RO, E1, res);
        }
        catch (...)
        {
            GERROR_STREAM("Exceptions happened in grappa2D::recon_fill_back_kspace(...) ... ");
            throw;
        }

        return res;
    }

    void grappa2D::help()
    {
        GDEBUG_STREAM("This class is to provide python binding for 2D grappa functionalities in Gadgetron");
        GDEBUG_STREAM("Please call initialize(...) first to set up grappa recon");
        GDEBUG_STREAM("Then call calib(...) and recon(...) for kspace grappa 2D reconstruction");
    }

    void grappa2D::status()
    {
        std::cout << "grappa2D of Gadgetron python binding" << std::endl;
        std::cout << "accelFactor is " << accelFactor_ << std::endl;
        std::cout << "kRO is " << kRO_ << std::endl;

        std::cout << "kE1 is [";
        size_t n;
        for (n=0; n<kE1_.size(); n++)
        {
            std::cout << kE1_[n];
        }
        std::cout << "]" << std::endl;

        std::cout << "oE1 is [";
        for (n = 0; n<oE1_.size(); n++)
        {
            std::cout << oE1_[n];
        }
        std::cout << "]" << std::endl;

        std::cout << "thres is " << thres_ << std::endl;
    }
}
