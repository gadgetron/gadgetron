
/** \file   mri_core_parallel_imaging.cpp
    \brief  Implement some functionalities for 2D and 3D MRI parallel imaging
    \author Hui Xue
*/

#include "mri_core_parallel_imaging.h"
#include "GadgetronTimer.h"
#include "mri_core_utility.h"

#ifdef USE_OMP
    #include "omp.h"
#endif // USE_OMP

namespace Gadgetron
{

template <typename T> 
void parallelImaging<T>::printInfo(std::ostream& os)
{
    using namespace std;

    os << "---------------------------- MRI core parallel imaging ---------------" << endl;
    os << "Implement some functionalities for 2D and 3D MRI parallel imaging" << endl;
    os << "----------------------------------------------------------------------" << endl;
}

template <typename T> 
void parallelImaging<T>::unmixCoeff(const hoNDArray<T>& kerIm, const hoNDArray<T>& coilMap, hoNDArray<T>& unmixCoeff, hoNDArray<T>& gFactor)
{
    try
    {
        size_t RO = kerIm.get_size(0);
        size_t E1 = kerIm.get_size(1);
        size_t srcCHA = kerIm.get_size(2);
        size_t dstCHA = kerIm.get_size(3);

        GADGET_CHECK_THROW(coilMap.get_size(0)==RO);
        GADGET_CHECK_THROW(coilMap.get_size(1)==E1);
        GADGET_CHECK_THROW(coilMap.get_size(2)==dstCHA);

        unmixCoeff.create(RO, E1, srcCHA);
        Gadgetron::clear(&unmixCoeff);

        gFactor.create(RO, E1);
        Gadgetron::clear(&gFactor);

        int src;

        T* pKerIm = const_cast<T*>(kerIm.begin());
        T* pCoilMap = const_cast<T*>(coilMap.begin());
        T* pCoeff = unmixCoeff.begin();

        std::vector<size_t> dim(2);
        dim[0] = RO;
        dim[1] = E1;

        #pragma omp parallel default(none) private(src) shared(RO, E1, srcCHA, dstCHA, pKerIm, pCoilMap, pCoeff, dim)
        {
            hoNDArray<T> coeff2D, coeffTmp(&dim);
            hoNDArray<T> coilMap2D;
            hoNDArray<T> kerIm2D;

            #pragma omp for
            for ( src=0; src<(int)srcCHA; src++ )
            {
                coeff2D.create(&dim, pCoeff+src*RO*E1);

                for ( size_t dst=0; dst<dstCHA; dst++ )
                {
                    kerIm2D.create(&dim, pKerIm+src*RO*E1+dst*RO*E1*srcCHA);
                    coilMap2D.create(&dim, pCoilMap+dst*RO*E1);
                    Gadgetron::multiplyConj(kerIm2D, coilMap2D, coeffTmp);
                    Gadgetron::add(coeff2D, coeffTmp, coeff2D);
                }
            }
        }

        hoNDArray<T> conjUnmixCoeff(unmixCoeff);
        Gadgetron::multiplyConj(unmixCoeff, conjUnmixCoeff, conjUnmixCoeff);
        Gadgetron::sumOverLastDimension(conjUnmixCoeff, gFactor);
        Gadgetron::sqrt(gFactor, gFactor);
    }
    catch(...)
    {
        GADGET_THROW("Errors in gtPlusReconWorker2DT<T>::unmixCoeff(const hoNDArray<T>& kerIm, const hoNDArray<T>& coilMap, hoNDArray<T>& unmixCoeff, hoNDArray<T>& gFactor) ... ");
    }
}

template <typename T> 
void parallelImaging<T>::applyUnmixCoeff(const hoNDArray<T>& kspace, const hoNDArray<T>& unmixCoeff, hoNDArray<T>& complexIm)
{
    try
    {
        GADGET_CHECK_THROW(kspace.get_size(0)==unmixCoeff.get_size(0));
        GADGET_CHECK_THROW(kspace.get_size(1)==unmixCoeff.get_size(1));
        GADGET_CHECK_THROW(kspace.get_size(2)==unmixCoeff.get_size(2));

        hoNDArray<T> buffer2DT(kspace);

        GADGET_CHECK_THROW(Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft2c(kspace, buffer2DT));
        this->applyUnmixCoeffImage(buffer2DT, unmixCoeff, complexIm);
    }
    catch(...)
    {
        GADGET_THROW("Errors in gtPlusReconWorker2DT<T>::applyUnmixCoeff(const hoNDArray<T>& kspace, const hoNDArray<T>& unmixCoeff, hoNDArray<T>& complexIm) ... ");
    }
}

template <typename T> 
void parallelImaging<T>::applyUnmixCoeffImage(const hoNDArray<T>& aliasedIm, const hoNDArray<T>& unmixCoeff, hoNDArray<T>& complexIm)
{
    try
    {
        GADGET_CHECK_THROW(aliasedIm.get_size(0)==unmixCoeff.get_size(0));
        GADGET_CHECK_THROW(aliasedIm.get_size(1)==unmixCoeff.get_size(1));
        GADGET_CHECK_THROW(aliasedIm.get_size(2)==unmixCoeff.get_size(2));

        boost::shared_ptr< std::vector<size_t> > dim = aliasedIm.get_dimensions();

        std::vector<size_t> dimIm(*dim);
        dimIm[2] = 1;

        if ( !complexIm.dimensions_equal(&dimIm) )
        {
            complexIm.create(&dimIm);
        }
        Gadgetron::clear(&complexIm);

        // hoNDArray<T> tmp(aliasedIm);
        hoNDArray<T> buffer2DT(aliasedIm);

        Gadgetron::multipleMultiply(unmixCoeff, aliasedIm, buffer2DT);
        Gadgetron::sumOver3rdDimension(buffer2DT, complexIm);
    }
    catch(...)
    {
        GADGET_THROW("Errors in gtPlusReconWorker2DT<T>::applyUnmixCoeffImage(const hoNDArray<T>& aliasedIm, const hoNDArray<T>& unmixCoeff, hoNDArray<T>& complexIm) ... ");
    }
}

template <typename T> 
void parallelImaging<T>::unmixCoeff3D(const hoNDArray<T>& kerIm, const hoNDArray<T>& coilMap, hoNDArray<T>& unmixCoeff, hoNDArray<T>& gFactor)
{
    try
    {
        size_t RO = kerIm.get_size(0);
        size_t E1 = kerIm.get_size(1);
        size_t E2 = kerIm.get_size(2);
        size_t srcCHA = kerIm.get_size(3);
        size_t dstCHA = kerIm.get_size(4);

        GADGET_CHECK_THROW(coilMap.get_size(0)==RO);
        GADGET_CHECK_THROW(coilMap.get_size(1)==E1);
        GADGET_CHECK_THROW(coilMap.get_size(2)==E2);
        GADGET_CHECK_THROW(coilMap.get_size(3)==dstCHA);

        unmixCoeff.create(RO, E1, E2, srcCHA);
        Gadgetron::clear(&unmixCoeff);
        gFactor.create(RO, E1, E2);
        Gadgetron::clear(&gFactor);

        int src;

        T* pKerIm = const_cast<T*>(kerIm.begin());
        T* pCoilMap = const_cast<T*>(coilMap.begin());
        T* pCoeff = unmixCoeff.begin();

        std::vector<size_t> dim(3);
        dim[0] = RO;
        dim[1] = E1;
        dim[2] = E2;

        #pragma omp parallel default(none) private(src) shared(RO, E1, E2, srcCHA, dstCHA, pKerIm, pCoilMap, pCoeff, dim)
        {
            hoNDArray<T> coeff2D, coeffTmp(&dim);
            hoNDArray<T> coilMap2D;
            hoNDArray<T> kerIm2D;

            #pragma omp for
            for ( src=0; src<(int)srcCHA; src++ )
            {
                coeff2D.create(&dim, pCoeff+src*RO*E1*E2);

                for ( size_t dst=0; dst<dstCHA; dst++ )
                {
                    kerIm2D.create(&dim, pKerIm+src*RO*E1*E2+dst*RO*E1*E2*srcCHA);
                    coilMap2D.create(&dim, pCoilMap+dst*RO*E1*E2);
                    Gadgetron::multiplyConj(kerIm2D, coilMap2D, coeffTmp);
                    Gadgetron::add(coeff2D, coeffTmp, coeff2D);
                }
            }
        }

        hoNDArray<T> conjUnmixCoeff(unmixCoeff);
        Gadgetron::multiplyConj(unmixCoeff, conjUnmixCoeff, conjUnmixCoeff);
        Gadgetron::sumOverLastDimension(conjUnmixCoeff, gFactor);
        Gadgetron::sqrt(gFactor, gFactor);
    }
    catch(...)
    {
        GADGET_THROW("Errors in gtPlusReconWorker3DT<T>::unmixCoeff(const hoNDArray<T>& kerIm, const hoNDArray<T>& coilMap, hoNDArray<T>& unmixCoeff, hoNDArray<T>& gFactor) ... ");
    }
}

template <typename T> 
void parallelImaging<T>::applyUnmixCoeff3D(const hoNDArray<T>& kspace, const hoNDArray<T>& unmixCoeff, hoNDArray<T>& complexIm)
{
    try
    {
        GADGET_CHECK_THROW(kspace.get_size(0)==unmixCoeff.get_size(0));
        GADGET_CHECK_THROW(kspace.get_size(1)==unmixCoeff.get_size(1));
        GADGET_CHECK_THROW(kspace.get_size(2)==unmixCoeff.get_size(2));
        GADGET_CHECK_THROW(kspace.get_size(3)==unmixCoeff.get_size(3));

        // buffer3DT_unwrapping_ = kspace;
        hoNDArray<T> buffer3DT(kspace.get_dimensions());

        GADGET_CHECK_THROW(Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft3c(kspace, buffer3DT));
        this->applyUnmixCoeffImage(buffer3DT, unmixCoeff, complexIm);
    }
    catch(...)
    {
        GADGET_THROW("Errors in gtPlusReconWorker3DT<T>::applyUnmixCoeff(const hoNDArray<T>& kspace, const hoNDArray<T>& unmixCoeff, hoNDArray<T>& complexIm) ... ");
    }
}

template <typename T> 
void parallelImaging<T>::applyUnmixCoeffImage3D(const hoNDArray<T>& aliasedIm, const hoNDArray<T>& unmixCoeff, hoNDArray<T>& complexIm)
{
    try
    {
        GADGET_CHECK_THROW(aliasedIm.get_size(0)==unmixCoeff.get_size(0));
        GADGET_CHECK_THROW(aliasedIm.get_size(1)==unmixCoeff.get_size(1));
        GADGET_CHECK_THROW(aliasedIm.get_size(2)==unmixCoeff.get_size(2));
        GADGET_CHECK_THROW(aliasedIm.get_size(3)==unmixCoeff.get_size(3));

        boost::shared_ptr< std::vector<size_t> > dim = aliasedIm.get_dimensions();
        std::vector<size_t> dimIm(*dim);
        dimIm[3] = 1;

        if ( !complexIm.dimensions_equal(&dimIm) )
        {
            complexIm.create(&dimIm);
        }
        Gadgetron::clear(&complexIm);

        // hoNDArray<T> tmp(aliasedIm);
        // buffer3DT_unwrapping_ = aliasedIm;

        hoNDArray<T> buffer3DT(aliasedIm.get_dimensions());

        Gadgetron::multipleMultiply(unmixCoeff, aliasedIm, buffer3DT);
        Gadgetron::sumOver4thDimension(buffer3DT, complexIm);
    }
    catch(...)
    {
        GADGET_THROW("Errors in gtPlusReconWorker3DT<T>::applyUnmixCoeffImage(const hoNDArray<T>& aliasedIm, const hoNDArray<T>& unmixCoeff, hoNDArray<T>& complexIm) ... ");
    }
}

// 
// Template instantiation
//

template class EXPORTMRICORE parallelImaging< std::complex<float> >;
template class EXPORTMRICORE parallelImaging< std::complex<double> >;

}
