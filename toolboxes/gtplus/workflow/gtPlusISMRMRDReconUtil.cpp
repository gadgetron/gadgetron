
#include "gtPlusISMRMRDReconUtil.h"

namespace Gadgetron { namespace gtPlus {

//
// Instantiation
//

template EXPORTGTPLUS class gtPlusISMRMRDReconUtil<float>;
template EXPORTGTPLUS class gtPlusISMRMRDReconUtil<double>;
template EXPORTGTPLUS class gtPlusISMRMRDReconUtil<GT_Complex8>;
template EXPORTGTPLUS class gtPlusISMRMRDReconUtil<GT_Complex16>;

template EXPORTGTPLUS class gtPlusISMRMRDReconUtilComplex<GT_Complex8>;
template EXPORTGTPLUS class gtPlusISMRMRDReconUtilComplex<GT_Complex16>;

}}

namespace Gadgetron {

    // ----------------------------------------------------------------------------------------
    // templated functions
    // ----------------------------------------------------------------------------------------

    template<typename T> 
    bool sumOverLastDimension(const hoNDArray<T>& x, hoNDArray<T>& r)
    {
        try
        {
            boost::shared_ptr< std::vector<size_t> > dim = x.get_dimensions();
            size_t NDim = dim->size();

            std::vector<size_t> dimR(NDim-1);

            size_t d;
            for ( d=0; d<NDim-1; d++ )
            {
                dimR[d] = (*dim)[d];
            }

            if ( !r.dimensions_equal(&dimR) )
            {
                r.create(&dimR);
            }

            // Gadgetron::clear(&r);

            if ( x.get_size(NDim-1) <= 1 )
            {
                memcpy(r.begin(), x.begin(), x.get_number_of_bytes());
                return true;
            }

            size_t lastDim = x.get_size(NDim-1);
            size_t NR = r.get_number_of_elements();
            T* pA = const_cast<T*>(x.begin());
            T* pR = r.begin();

            memcpy(pR, pA, sizeof(T)*NR);

            // sum over the last dim
            hoNDArray<T> tmp;
            for ( d=1; d<lastDim; d++ )
            {
                tmp.create(&dimR, pA+d*NR);
                add(tmp, r, r);
            }
        }
        catch (...)
        {
            GADGET_ERROR_MSG("Errors in sumOverLastDimension(const hoNDArray<T>& x, hoNDArray<T>& r) ... ");
            return false;
        }
        return true;
    }

    template<typename T> 
    bool sumOverSecondLastDimension(const hoNDArray<T>& x, hoNDArray<T>& r)
    {
        try
        {
            boost::shared_ptr< std::vector<size_t> > dim = x.get_dimensions();
            size_t NDim = dim->size();

            if ( NDim < 2 ) return true;

            std::vector<size_t> dimR(NDim-1);
            std::vector<size_t> dimRInternal(NDim-2);

            size_t d;
            for ( d=0; d<NDim-2; d++ )
            {
                dimR[d] = (*dim)[d];
                dimRInternal[d] = (*dim)[d];
            }
            dimR[NDim-2] = (*dim)[NDim-1];

            if ( !r.dimensions_equal(&dimR) )
            {
                r.create(&dimR);
            }

            if ( x.get_size(NDim-2) <= 1 )
            {
                memcpy(r.begin(), x.begin(), x.get_number_of_bytes());
                return true;
            }

            size_t lastDim = x.get_size(NDim-1);
            size_t secondLastDim = x.get_size(NDim-2);
            size_t NS = x.get_number_of_elements()/lastDim;
            size_t NR = r.get_number_of_elements()/lastDim;
            T* pA = const_cast<T*>(x.begin());
            T* pR = r.begin();

            int l;
            #pragma omp parallel default(none) private(l) shared(lastDim, secondLastDim, NS, NR, pA, pR, dimRInternal)
            {
                hoNDArray<T> tmp, tmp2;

                #pragma omp for
                for ( l=0; l<(int)lastDim; l++ )
                {
                    memcpy(pR+l*NR, pA+l*NS, sizeof(T)*NR);
                    tmp.create(&dimRInternal, pR+l*NR);
                    for ( size_t s=1; s<secondLastDim; s++ )
                    {
                        tmp2.create(&dimRInternal, pA+l*NS+s*NR);
                        add(tmp, tmp2, tmp);
                    }
                }
            }
        }
        catch (...)
        {
            GADGET_ERROR_MSG("Errors in sumOverSecondLastDimension(const hoNDArray<T>& x, hoNDArray<T>& r) ... ");
            return false;
        }
        return true;
    }

    // e.g. x is 3D and y is 4D array, r(:,:,:,n) = y(:,:,:,n) .* x
    template<typename T> 
    bool multiplyOverLastDimension(const hoNDArray<T>& x, const hoNDArray<T>& y, hoNDArray<T>& r)
    {
        try
        {
            boost::shared_ptr< std::vector<size_t> > dimX = x.get_dimensions();
            boost::shared_ptr< std::vector<size_t> > dimY = y.get_dimensions();

            size_t NDim = dimY->size();

            GADGET_CHECK_RETURN_FALSE(dimX->size()==NDim-1);

            if ( !r.dimensions_equal(dimY.get()) )
            {
                r.create(dimY);
            }

            if ( y.get_size(NDim-1) <= 1 )
            {
                GADGET_CHECK_RETURN_FALSE(multiply(x, y, r));
                return true;
            }

            size_t lastDim = y.get_size(NDim-1);
            size_t N = x.get_number_of_elements();
            const T* pX = x.begin();
            const T* pY = y.begin();
            T* pR = r.begin();

            int d;

            #ifdef GCC_OLD_FLAG
                #pragma omp parallel default(none) private(d) shared(dimX, lastDim, N, pY, pR)
            #else
                #pragma omp parallel default(none) private(d) shared(x, dimX, lastDim, N, pY, pR)
            #endif
            {
                hoNDArray<T> tmpY, tmpR;

                #pragma omp for
                for ( d=0; d<(int)lastDim; d++ )
                {
                    tmpY.create(dimX.get(), const_cast<T*>(pY+d*N));
                    tmpR.create(dimX.get(), pR+d*N);
                    multiply(x, tmpY, tmpR);
                }
            }
        }
        catch (...)
        {
            GADGET_ERROR_MSG("Errors in multiplyOverLastDimension(const hoNDArray<float>& x, const hoNDArray<float>& y, hoNDArray<float>& r) ... ");
            return false;
        }
        return true;
    }

    // e.g. x is 3D and y is 4D array, r(:,:,:,n) = y(:,:,:,n) ./ x
    template<typename T> 
    bool divideOverLastDimension(const hoNDArray<T>& x, const hoNDArray<T>& y, hoNDArray<T>& r)
    {
        try
        {
            boost::shared_ptr< std::vector<size_t> > dimX = x.get_dimensions();
            boost::shared_ptr< std::vector<size_t> > dimY = y.get_dimensions();

            size_t NDim = dimY->size();

            GADGET_CHECK_RETURN_FALSE(dimX->size()==NDim-1);

            if ( !r.dimensions_equal(dimY.get()) )
            {
                r.create(dimY);
            }

            if ( y.get_size(NDim-1) <= 1 )
            {
                GADGET_CHECK_RETURN_FALSE(divide(y, x, r));
                return true;
            }

            size_t lastDim = y.get_size(NDim-1);
            size_t N = x.get_number_of_elements();
            T* pY = const_cast<T*>(y.begin());
            T* pR = r.begin();

            int d;

            #ifdef GCC_OLD_FLAG
                #pragma omp parallel default(none) private(d) shared(dimX, lastDim, N, pY, pR)
            #else
                #pragma omp parallel default(none) private(d) shared(x, dimX, lastDim, N, pY, pR)
            #endif
            {
                hoNDArray<T> tmpY, tmpR;

                #pragma omp for
                for ( d=0; d<(int)lastDim; d++ )
                {
                    tmpY.create(dimX, pY+d*N);
                    tmpR.create(dimX, pR+d*N);
                    divide(tmpY, x, tmpR);
                }
            }
        }
        catch (...)
        {
            GADGET_ERROR_MSG("Errors in divideOverLastDimension(const hoNDArray<T>& x, const hoNDArray<T>& y, hoNDArray<T>& r) ... ");
            return false;
        }
        return true;
    }

    template<typename T> 
    bool sumOver1stDimension(const hoNDArray<T>& x, hoNDArray<T>& r)
    {
        try
        {
            size_t RO = x.get_size(0);
            size_t num = x.get_number_of_elements()/(RO);

            boost::shared_ptr< std::vector<size_t> > dim = x.get_dimensions();

            std::vector<size_t> dimAve(*dim);
            dimAve[0] = 1;
            r.create(&dimAve);

            const T* pX = x.begin();
            T* pR = r.begin();

            int n;
            #pragma omp parallel for default(none) private(n) shared(RO, num, pX, pR)
            for ( n=0; n<(int)num; n++ )
            {
                T xsum = pX[n*RO];
                for (size_t ro=1; ro<RO; ro++ )
                {
                    xsum += pX[n*RO+ro];
                }

                pR[n] = xsum;
            }
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Errors in sumOver1stDimension(...) ... ");
            return false;
        }

        return true;
    }

    template<typename T> 
    bool sumOver2ndDimension(const hoNDArray<T>& x, hoNDArray<T>& r)
    {
        try
        {
            size_t NDim = x.get_number_of_dimensions();

            if ( NDim < 2 )
            {
                r = x;
                return true;
            }

            size_t RO = x.get_size(0);
            size_t E1 = x.get_size(1);

            size_t num = x.get_number_of_elements()/(RO*E1);

            boost::shared_ptr< std::vector<size_t> > dim = x.get_dimensions();

            std::vector<size_t> dimAve(*dim);
            dimAve[1] = 1;
            r.create(&dimAve);

            int n;
            #ifdef GCC_OLD_FLAG
                #pragma omp parallel for default(none) private(n) shared(RO, E1, num)
            #else
                #pragma omp parallel for default(none) private(n) shared(RO, E1, num, x, r)
            #endif
            for ( n=0; n<(int)num; n++ )
            {
                hoNDArray<T> xsum(RO, const_cast<T*>(r.begin()+n*RO));
                memcpy(xsum.begin(), x.begin()+n*RO*E1, xsum.get_number_of_bytes());

                for (size_t e1=1; e1<E1; e1++ )
                {
                    hoNDArray<T> x1D(RO, const_cast<T*>(x.begin()+n*RO*E1+e1*RO));
                    Gadgetron::add(x1D, xsum, xsum);
                }
            }
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Errors in sumOver2ndDimension(...) ... ");
            return false;
        }

        return true;
    }

    template<typename T> 
    bool sumOver3rdDimension(const hoNDArray<T>& x, hoNDArray<T>& r)
    {
        try
        {
            size_t NDim = x.get_number_of_dimensions();

            if ( NDim < 3 )
            {
                r = x;
                return true;
            }

            size_t RO = x.get_size(0);
            size_t E1 = x.get_size(1);
            size_t CHA = x.get_size(2);

            size_t num = x.get_number_of_elements()/(RO*E1*CHA);

            boost::shared_ptr< std::vector<size_t> > dim = x.get_dimensions();

            std::vector<size_t> dimAve(*dim);
            dimAve[2] = 1;
            r.create(&dimAve);

            int n;
            #ifdef GCC_OLD_FLAG
                #pragma omp parallel default(none) private(n) shared(RO, E1, CHA, num) if (num>1)
            #else
                #pragma omp parallel default(none) private(n) shared(RO, E1, CHA, num, x, r) if (num>1)
            #endif
            {
                hoNDArray<T> xsum;
                hoNDArray<T> x2D;

                #pragma omp for
                for ( n=0; n<(int)num; n++ )
                {
                    xsum.create(RO, E1, const_cast<T*>(r.begin()+n*RO*E1));
                    memcpy(xsum.begin(), x.begin()+n*RO*E1*CHA, xsum.get_number_of_bytes());

                    for (size_t cha=1; cha<CHA; cha++ )
                    {
                        x2D.create(RO, E1, const_cast<T*>(x.begin()+n*RO*E1*CHA+cha*RO*E1));
                        Gadgetron::add(x2D, xsum, xsum);
                    }
                }
            }
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Errors in sumOver3rdDimension(...) ... ");
            return false;
        }

        return true;
    }

    template<typename T> bool sumOver4thDimension(const hoNDArray<T>& x, hoNDArray<T>& r)
    {
        try
        {
            size_t NDim = x.get_number_of_dimensions();

            if ( NDim < 4 )
            {
                r = x;
                return true;
            }

            size_t RO = x.get_size(0);
            size_t E1 = x.get_size(1);
            size_t CHA = x.get_size(2);
            size_t N = x.get_size(3);

            size_t num = x.get_number_of_elements()/(RO*E1*CHA*N);

            boost::shared_ptr< std::vector<size_t> > dim = x.get_dimensions();

            std::vector<size_t> dimAve(*dim);
            dimAve[3] = 1;
            r.create(&dimAve);

            int n;
            #ifdef GCC_OLD_FLAG
                #pragma omp parallel for default(none) private(n) shared(RO, E1, CHA, N, num)
            #else
                #pragma omp parallel for default(none) private(n) shared(RO, E1, CHA, N, num, x, r)
            #endif
            for ( n=0; n<(int)num; n++ )
            {
                hoNDArray<T> xsum(RO, E1, CHA, const_cast<T*>(r.begin()+n*RO*E1*CHA));
                memcpy(xsum.begin(), x.begin()+n*RO*E1*CHA*N, xsum.get_number_of_bytes());

                for (size_t nn=1; nn<N; nn++ )
                {
                    hoNDArray<T> x3D(RO, E1, CHA, const_cast<T*>(x.begin()+n*RO*E1*CHA*N+nn*RO*E1*CHA));
                    Gadgetron::add(x3D, xsum, xsum);
                }
            }
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Errors in sumOver4thDimension(const hoNDArray<T>& x, hoNDArray<T>& r) ... ");
            return false;
        }

        return true;
    }

    template<typename T> bool sumOver5thDimension(const hoNDArray<T>& x, hoNDArray<T>& r)
    {
        try
        {
            size_t NDim = x.get_number_of_dimensions();

            if ( NDim < 5 )
            {
                r = x;
                return true;
            }

            size_t RO = x.get_size(0);
            size_t E1 = x.get_size(1);
            size_t CHA = x.get_size(2);
            size_t N = x.get_size(3);
            size_t S = x.get_size(4);

            size_t num = x.get_number_of_elements()/(RO*E1*CHA*N*S);

            boost::shared_ptr< std::vector<size_t> > dim = x.get_dimensions();

            std::vector<size_t> dimAve(*dim);
            dimAve[4] = 1;
            r.create(&dimAve);

            int n;
            #ifdef GCC_OLD_FLAG
                #pragma omp parallel for default(none) private(n) shared(RO, E1, CHA, N, S, num) if (num > 4)
            #else
                #pragma omp parallel for default(none) private(n) shared(RO, E1, CHA, N, S, num, x, r) if (num > 4)
            #endif
            for ( n=0; n<(int)num; n++ )
            {
                hoNDArray<T> xsum(RO, E1, CHA, N, const_cast<T*>(r.begin()+n*RO*E1*CHA*N));
                memcpy(xsum.begin(), x.begin()+n*RO*E1*CHA*N*S, xsum.get_number_of_bytes());

                for (size_t s=1; s<S; s++ )
                {
                    hoNDArray<T> x4D(RO, E1, CHA, N, const_cast<T*>(x.begin()+n*RO*E1*CHA*N*S+s*RO*E1*CHA*N));
                    Gadgetron::add(x4D, xsum, xsum);
                }
            }
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Errors in sumOver5thDimension(const hoNDArray<T>& x, hoNDArray<T>& r) ... ");
            return false;
        }

        return true;
    }

    // e.g. x is 3D and y is 4D array, r(:,:,n,:) = y(:,:,n,:) .* x3D
    template<typename T> 
    bool multiplyOver3rdDimension(const hoNDArray<T>& x3D, const hoNDArray<T>& y4D, hoNDArray<T>& r)
    {
        try
        {
            boost::shared_ptr< std::vector<size_t> > dimX = x3D.get_dimensions();
            boost::shared_ptr< std::vector<size_t> > dimY = y4D.get_dimensions();

            size_t NDim = dimY->size();

            GADGET_CHECK_RETURN_FALSE(dimX->size()>=3);
            GADGET_CHECK_RETURN_FALSE(NDim>=4);
            GADGET_CHECK_RETURN_FALSE((*dimX)[0]==(*dimY)[0]);
            GADGET_CHECK_RETURN_FALSE((*dimX)[1]==(*dimY)[1]);

            if ( !r.dimensions_equal(dimY.get()) )
            {
                r.create(dimY);
            }

            long long t, N2D = (long long)x3D.get_size(0)*x3D.get_size(1);
            long long sz = (long long)y4D.get_size(2);
            long long st = (long long)y4D.get_number_of_elements()/(N2D*sz);

            if ( sz == 1 )
            {
                GADGET_CHECK_RETURN_FALSE(multiply(x3D, y4D, r));
                return true;
            }

            const T* pX = x3D.begin();
            const T* pY = y4D.begin();
            T* pR = r.begin();

            std::vector<size_t> dim2D(2);
            dim2D[0] = (*dimY)[0];
            dim2D[1] = (*dimY)[1];

            #pragma omp parallel for default(none) private(t) shared(N2D, sz, st, dim2D, pX, pY, pR)
            for ( t=0; t<st; t++ )
            {
                hoNDArray<T> tmpX, tmpY, tmpR;
                tmpX.create(&dim2D, const_cast<T*>(pX+t*N2D));

                for ( long long z=0; z<sz; z++ )
                {
                    tmpY.create(&dim2D, const_cast<T*>(pY+t*N2D*sz+z*N2D));
                    tmpR.create(&dim2D, pR+t*N2D*sz+z*N2D);
                    multiply(tmpX, tmpY, tmpR);
                }
            }
        }
        catch (...)
        {
            GADGET_ERROR_MSG("Errors in multiplyOver3rdDimension(const hoNDArray<float>& x3D, const hoNDArray<float>& y4D, hoNDArray<float>& r) ... ");
            return false;
        }
        return true;
    }

    template<typename T> 
    bool multiplyOver4thDimension(const hoNDArray<T>& x4D, const hoNDArray<T>& y5D, hoNDArray<T>& r)
    {
        try
        {
            boost::shared_ptr< std::vector<size_t> > dimX = x4D.get_dimensions();
            boost::shared_ptr< std::vector<size_t> > dimY = y5D.get_dimensions();

            size_t NDim = dimY->size();

            GADGET_CHECK_RETURN_FALSE(dimX->size()>=4);
            GADGET_CHECK_RETURN_FALSE((*dimX)[0]==(*dimY)[0]);
            GADGET_CHECK_RETURN_FALSE((*dimX)[1]==(*dimY)[1]);
            GADGET_CHECK_RETURN_FALSE((*dimX)[2]==(*dimY)[2]);

            if ( !r.dimensions_equal(dimY.get()) )
            {
                r.create(dimY);
            }

            size_t RO = (*dimX)[0];
            size_t E1 = (*dimX)[1];
            size_t CHA = (*dimX)[2];

            long long t, N3D = (long long)RO*E1*CHA;

            size_t N = (*dimY)[3];
            size_t num = x4D.get_number_of_elements()/(RO*E1*CHA);

            const T* pX = x4D.begin();
            const T* pY = y5D.begin();
            T* pR = r.begin();

            std::vector<size_t> dim3D(3);
            dim3D[0] = RO;
            dim3D[1] = E1;
            dim3D[2] = CHA;

            #pragma omp parallel for default(none) private(t) shared(N3D, N, dim3D, pX, pY, pR, num)
            for ( t=0; t<(int)num; t++ )
            {
                hoNDArray<T> tmpX, tmpY, tmpR;
                tmpX.create(&dim3D, const_cast<T*>(pX+t*N3D));

                for ( int n=0; n<N; n++ )
                {
                    tmpY.create(&dim3D, const_cast<T*>(pY+t*N3D*N+n*N3D));
                    tmpR.create(&dim3D, pR+t*N3D*N+n*N3D);
                    multiply(tmpX, tmpY, tmpR);
                }
            }
        }
        catch (...)
        {
            GADGET_ERROR_MSG("Errors in multiplyOver4thDimension(const hoNDArray<float>& x4D, const hoNDArray<float>& y5D, hoNDArray<float>& r) ... ");
            return false;
        }
        return true;
    }

    template<typename T> 
    bool multiplyOver4thDimensionExcept(const hoNDArray<T>& x4D, const hoNDArray<T>& y5D, size_t n, hoNDArray<T>& r, bool copyY2R)
    {
        try
        {
            boost::shared_ptr< std::vector<size_t> > dimX = x4D.get_dimensions();
            boost::shared_ptr< std::vector<size_t> > dimY = y5D.get_dimensions();

            size_t NDim = dimY->size();

            GADGET_CHECK_RETURN_FALSE(dimX->size()>=4);
            GADGET_CHECK_RETURN_FALSE((*dimX)[0]==(*dimY)[0]);
            GADGET_CHECK_RETURN_FALSE((*dimX)[1]==(*dimY)[1]);
            GADGET_CHECK_RETURN_FALSE((*dimX)[2]==(*dimY)[2]);

            const T* pX = x4D.begin();
            const T* pY = y5D.begin();
            T* pR = r.begin();

            if ( (pR!=pY) && (!r.dimensions_equal(dimY.get())) )
            {
                r.create(dimY);
                pR = r.begin();
            }

            size_t RO = (*dimX)[0];
            size_t E1 = (*dimX)[1];
            size_t CHA = (*dimX)[2];

            long long t, N3D = (long long)RO*E1*CHA;

            long long N = (long long)(*dimY)[3];
            long long num = (long long)x4D.get_number_of_elements()/(RO*E1*CHA);

            std::vector<size_t> dim3D(3);
            dim3D[0] = RO;
            dim3D[1] = E1;
            dim3D[2] = CHA;

            #pragma omp parallel for default(none) private(t) shared(N3D, N, dim3D, pX, pY, pR, num, n, copyY2R)
            for ( t=0; t<num; t++ )
            {
                hoNDArray<T> tmpX, tmpY, tmpR;
                tmpX.create(&dim3D, const_cast<T*>(pX+t*N3D));

                for ( long long z=0; z<N; z++ )
                {
                    if ( z != n )
                    {
                        tmpY.create(&dim3D, const_cast<T*>(pY+t*N3D*N+z*N3D));
                        tmpR.create(&dim3D, pR+t*N3D*N+z*N3D);
                        multiply(tmpX, tmpY, tmpR);
                    }
                    else
                    {
                        if ( pR != pY )
                        {
                            if ( copyY2R )
                            {
                                memcpy(pR+t*N3D*N+z*N3D, const_cast<T*>(pY+t*N3D*N+z*N3D), sizeof(T)*N3D);
                            }
                        }
                    }
                }
            }
        }
        catch (...)
        {
            GADGET_ERROR_MSG("Errors in multiplyOver4thDimensionExcept(const hoNDArray<float>& x4D, const hoNDArray<float>& y5D, size_t n, hoNDArray<float>& r, bool copyY2R) ... ");
            return false;
        }
        return true;
    }

    template<typename T>
    bool multiplyOver5thDimension(const hoNDArray<T>& x, const hoNDArray<T>& y, hoNDArray<T>& r)
    {
        try
        {
            boost::shared_ptr< std::vector<size_t> > dimX = x.get_dimensions();
            boost::shared_ptr< std::vector<size_t> > dimY = y.get_dimensions();

            size_t NDim = dimY->size();

            GADGET_CHECK_RETURN_FALSE(dimX->size()>=5);
            GADGET_CHECK_RETURN_FALSE((*dimX)[0]==(*dimY)[0]);
            GADGET_CHECK_RETURN_FALSE((*dimX)[1]==(*dimY)[1]);
            GADGET_CHECK_RETURN_FALSE((*dimX)[2]==(*dimY)[2]);
            GADGET_CHECK_RETURN_FALSE((*dimX)[3]==(*dimY)[3]);

            if ( !r.dimensions_equal(dimY.get()) )
            {
                r.create(dimY);
            }

            size_t RO = (*dimX)[0];
            size_t E1 = (*dimX)[1];
            size_t E2 = (*dimX)[2];
            size_t CHA = (*dimX)[3];

            int t;
            size_t N4D = RO*E1*E2*CHA;

            size_t N = (*dimY)[4];
            size_t num = x.get_number_of_elements()/N4D;

            const T* pX = x.begin();
            const T* pY = y.begin();
            T* pR = r.begin();

            std::vector<size_t> dim4D(4);
            dim4D[0] = RO;
            dim4D[1] = E1;
            dim4D[2] = E2;
            dim4D[3] = CHA;

            #pragma omp parallel for default(none) private(t) shared(N4D, N, dim4D, pX, pY, pR, num)
            for ( t=0; t<(int)num; t++ )
            {
                hoNDArray<T> tmpX, tmpY, tmpR;
                tmpX.create(&dim4D, const_cast<T*>(pX+t*N4D));

                for ( int n=0; n<N; n++ )
                {
                    tmpY.create(&dim4D, const_cast<T*>(pY+t*N4D*N+n*N4D));
                    tmpR.create(&dim4D, pR+t*N4D*N+n*N4D);
                    multiply(tmpX, tmpY, tmpR);
                }
            }
        }
        catch (...)
        {
            GADGET_ERROR_MSG("Errors in multiplyOver5thDimension(const hoNDArray<float>& x, const hoNDArray<float>& y, hoNDArray<float>& r) ... ");
            return false;
        }
        return true;
    }

    template<typename T> 
    bool multiplyOver5thDimensionExcept(const hoNDArray<T>& x, const hoNDArray<T>& y, size_t n, hoNDArray<T>& r, bool copyY2R)
    {
        try
        {
            boost::shared_ptr< std::vector<size_t> > dimX = x.get_dimensions();
            boost::shared_ptr< std::vector<size_t> > dimY = y.get_dimensions();

            size_t NDim = dimY->size();

            GADGET_CHECK_RETURN_FALSE(dimX->size()>=5);
            GADGET_CHECK_RETURN_FALSE((*dimX)[0]==(*dimY)[0]);
            GADGET_CHECK_RETURN_FALSE((*dimX)[1]==(*dimY)[1]);
            GADGET_CHECK_RETURN_FALSE((*dimX)[2]==(*dimY)[2]);
            GADGET_CHECK_RETURN_FALSE((*dimX)[3]==(*dimY)[3]);

            const T* pX = x.begin();
            const T* pY = y.begin();
            T* pR = r.begin();

            if ( (pR!=pY) && (!r.dimensions_equal(dimY.get())) )
            {
                r.create(dimY);
                pR = r.begin();
            }

            size_t RO = (*dimX)[0];
            size_t E1 = (*dimX)[1];
            size_t E2 = (*dimX)[2];
            size_t CHA = (*dimX)[3];

            int t;
            size_t N4D = RO*E1*E2*CHA;

            size_t N = (*dimY)[4];
            size_t num = x.get_number_of_elements()/N4D;

            std::vector<size_t> dim4D(4);
            dim4D[0] = RO;
            dim4D[1] = E1;
            dim4D[2] = E2;
            dim4D[3] = CHA;

            #pragma omp parallel for default(none) private(t) shared(N4D, dim4D, pX, pY, pR, num, n, N, copyY2R)
            for ( t=0; t<(int)num; t++ )
            {
                hoNDArray<T> tmpX, tmpY, tmpR;
                tmpX.create(&dim4D, const_cast<T*>(pX+t*N4D));

                for ( int z=0; z<N; z++ )
                {
                    if ( z != n )
                    {
                        tmpY.create(&dim4D, const_cast<T*>(pY+t*N4D*N+z*N4D));
                        tmpR.create(&dim4D, pR+t*N4D*N+z*N4D);
                        multiply(tmpX, tmpY, tmpR);
                    }
                    else
                    {
                        if ( pR != pY )
                        {
                            if ( copyY2R )
                            {
                                memcpy(pR+t*N4D*N+z*N4D, const_cast<T*>(pY+t*N4D*N+z*N4D), sizeof(T)*N4D);
                            }
                        }
                    }
                }
            }
        }
        catch (...)
        {
            GADGET_ERROR_MSG("Errors in multiplyOver5thDimensionExcept(const hoNDArray<T>& x, const hoNDArray<T>& y, size_t n, hoNDArray<T>& r, bool copyY2R) ... ");
            return false;
        }
        return true;
    }

    template <typename T> 
    bool multipleAdd(const hoNDArray<T>& x, const hoNDArray<T>& y, hoNDArray<T>& r)
    {
        GADGET_DEBUG_CHECK_RETURN_FALSE(x.get_number_of_elements()<=y.get_number_of_elements());
        if ( r.get_number_of_elements()!=y.get_number_of_elements())
        {
            r = y;
        }

        long long Nx = x.get_number_of_elements();
        long long N = y.get_number_of_elements() / Nx;

        long long n;

        if ( typeid(T)==typeid(float) )
        {
            #ifdef GCC_OLD_FLAG
                #pragma omp parallel for default(none) private(n) shared(Nx, N)
            #else
                #pragma omp parallel for default(none) private(n) shared(x, y, r, Nx, N)
            #endif
            for ( n=0; n<N; n++ )
            {
                //#ifdef USE_MKL
                //    vsAdd(x.get_number_of_elements(), reinterpret_cast<const float*>(x.begin()), reinterpret_cast<const float*>(y.begin()+n*Nx), reinterpret_cast<float*>(r.begin()+n*Nx));
                //#else
                    Gadgetron::math::add(x.get_number_of_elements(), x.begin(), y.begin()+n*Nx, r.begin()+n*Nx);
                //#endif // USE_MKL
            }
        }
        else if ( typeid(T)==typeid(double) )
        {
            #ifdef GCC_OLD_FLAG
                #pragma omp parallel for default(none) private(n) shared(Nx, N)
            #else
                #pragma omp parallel for default(none) private(n) shared(x, y, r, Nx, N)
            #endif
            for ( n=0; n<N; n++ )
            {
                //#ifdef USE_MKL
                //    vdAdd(x.get_number_of_elements(), reinterpret_cast<const double*>(x.begin()), reinterpret_cast<const double*>(y.begin()+n*Nx), reinterpret_cast<double*>(r.begin()+n*Nx));
                //#else
                    Gadgetron::math::add(x.get_number_of_elements(), x.begin(), y.begin()+n*Nx, r.begin()+n*Nx);
                //#endif // USE_MKL
            }
        }
        else if ( typeid(T)==typeid(GT_Complex8) )
        {
            #ifdef GCC_OLD_FLAG
                #pragma omp parallel for default(none) private(n) shared(Nx, N)
            #else
                #pragma omp parallel for default(none) private(n) shared(x, y, r, Nx, N)
            #endif
            for ( n=0; n<N; n++ )
            {
                //#ifdef USE_MKL
                //    vcAdd(x.get_number_of_elements(), reinterpret_cast<const MKL_Complex8*>(x.begin()), reinterpret_cast<const MKL_Complex8*>(y.begin()+n*Nx), reinterpret_cast<MKL_Complex8*>(r.begin()+n*Nx));
                //#else
                    Gadgetron::math::add(x.get_number_of_elements(), x.begin(), y.begin()+n*Nx, r.begin()+n*Nx);
                //#endif // USE_MKL
            }
        }
        else if ( typeid(T)==typeid(GT_Complex16) )
        {
            #ifdef GCC_OLD_FLAG
                #pragma omp parallel for default(none) private(n) shared(Nx, N)
            #else
                #pragma omp parallel for default(none) private(n) shared(x, y, r, Nx, N)
            #endif
            for ( n=0; n<N; n++ )
            {
                //#ifdef USE_MKL
                //    vzAdd(x.get_number_of_elements(), reinterpret_cast<const MKL_Complex16*>(x.begin()), reinterpret_cast<const MKL_Complex16*>(y.begin()+n*Nx), reinterpret_cast<MKL_Complex16*>(r.begin()+n*Nx));
                //#else
                    Gadgetron::math::add(x.get_number_of_elements(), x.begin(), y.begin()+n*Nx, r.begin()+n*Nx);
                //#endif // USE_MKL
            }
        }
        else
        {
            GADGET_ERROR_MSG("multipleAdd : unsupported type " << typeid(T).name());
            return false;
        }

        //#ifdef USE_MKL
        //    GADGET_CHECK_RETURN_FALSE(vmlGetErrStatus()==0);
        //#endif // USE_MKL

        return true;
    }

    template <typename T> 
    bool multipleMultiply(const hoNDArray<T>& x, const hoNDArray<T>& y, hoNDArray<T>& r)
    {
        GADGET_DEBUG_CHECK_RETURN_FALSE(x.get_number_of_elements()<=y.get_number_of_elements());
        if ( r.get_number_of_elements()!=y.get_number_of_elements())
        {
            r = y;
        }

        long long Nx = x.get_number_of_elements();
        long long N = y.get_number_of_elements() / Nx;

        long long n;

        if ( typeid(T)==typeid(float) )
        {
            #ifdef GCC_OLD_FLAG
                #pragma omp parallel for default(none) private(n) shared(Nx, N)
            #else
                #pragma omp parallel for default(none) private(n) shared(x, y, r, Nx, N)
            #endif
            for ( n=0; n<N; n++ )
            {
                //#ifdef USE_MKL
                //    vsMul(x.get_number_of_elements(), reinterpret_cast<const float*>(x.begin()), reinterpret_cast<const float*>(y.begin()+n*Nx), reinterpret_cast<float*>(r.begin()+n*Nx));
                //#else
                    Gadgetron::math::multiply(x.get_number_of_elements(), x.begin(), y.begin()+n*Nx, r.begin()+n*Nx);
                //#endif // USE_MKL
            }
        }
        else if ( typeid(T)==typeid(double) )
        {
            #ifdef GCC_OLD_FLAG
                #pragma omp parallel for default(none) private(n) shared(Nx, N)
            #else
                #pragma omp parallel for default(none) private(n) shared(x, y, r, Nx, N)
            #endif
            for ( n=0; n<N; n++ )
            {
                //#ifdef USE_MKL
                //    vdMul(x.get_number_of_elements(), reinterpret_cast<const double*>(x.begin()), reinterpret_cast<const double*>(y.begin()+n*Nx), reinterpret_cast<double*>(r.begin()+n*Nx));
                //#else
                    Gadgetron::math::multiply(x.get_number_of_elements(), x.begin(), y.begin()+n*Nx, r.begin()+n*Nx);
                //#endif // USE_MKL
            }
        }
        else if ( typeid(T)==typeid(GT_Complex8) )
        {
            #ifdef GCC_OLD_FLAG
                #pragma omp parallel for default(none) private(n) shared(Nx, N)
            #else
                #pragma omp parallel for default(none) private(n) shared(x, y, r, Nx, N)
            #endif
            for ( n=0; n<N; n++ )
            {
                //#ifdef USE_MKL
                //    vcMul(x.get_number_of_elements(), reinterpret_cast<const MKL_Complex8*>(x.begin()), reinterpret_cast<const MKL_Complex8*>(y.begin()+n*Nx), reinterpret_cast<MKL_Complex8*>(r.begin()+n*Nx));
                //#else
                    Gadgetron::math::multiply(x.get_number_of_elements(), x.begin(), y.begin()+n*Nx, r.begin()+n*Nx);
                //#endif // USE_MKL
            }
        }
        else if ( typeid(T)==typeid(GT_Complex16) )
        {
            #ifdef GCC_OLD_FLAG
                #pragma omp parallel for default(none) private(n) shared(Nx, N)
            #else
                #pragma omp parallel for default(none) private(n) shared(x, y, r, Nx, N)
            #endif
            for ( n=0; n<N; n++ )
            {
                //#ifdef USE_MKL
                //    vzMul(x.get_number_of_elements(), reinterpret_cast<const MKL_Complex16*>(x.begin()), reinterpret_cast<const MKL_Complex16*>(y.begin()+n*Nx), reinterpret_cast<MKL_Complex16*>(r.begin()+n*Nx));
                //#else
                    Gadgetron::math::multiply(x.get_number_of_elements(), x.begin(), y.begin()+n*Nx, r.begin()+n*Nx);
                //#endif // USE_MKL
            }
        }
        else
        {
            GADGET_ERROR_MSG("multipleMultiply : unsupported type " << typeid(T).name());
            return false;
        }

        //#ifdef USE_MKL
        //    GADGET_CHECK_RETURN_FALSE(vmlGetErrStatus()==0);
        //#endif // USE_MKL

        return true;
    }

    template <typename T> 
    bool cropUpTo11DArray(const hoNDArray<T>& x, hoNDArray<T>& r, const std::vector<size_t>& startND, std::vector<size_t>& size)
    {
        GADGET_CHECK_RETURN_FALSE( startND.size() == size.size() );
        GADGET_CHECK_RETURN_FALSE( startND.size() <= 11 );

        r.create(&size);
        if ( r.get_number_of_elements() == x.get_number_of_elements() )
        {
            r = x;
            return true;
        }

        std::vector<size_t> start(11, 0);
        std::vector<size_t> end(11, 0);

        size_t ii;
        for ( ii=0; ii<startND.size(); ii++ )
        {
            start[ii] = startND[ii];
            end[ii] = start[ii] + size[ii] - 1;
            GADGET_CHECK_RETURN_FALSE(end[ii] < x.get_size(ii));
        }

        // [Ro E1 Cha Slice E2 Con Phase Rep Set Seg Ave]
        size_t e1, cha, n, s, con, phs, rep, set, seg, ave;

        std::vector<size_t> srcInd(11), dstInd(11);

        for ( ave=start[10]; ave<=end[10]; ave++ )
        {
            srcInd[10] = ave; dstInd[10] = ave-start[10];

            for ( seg=start[9]; seg<=end[9]; seg++ )
            {
                srcInd[9] = seg; dstInd[9] = seg-start[9];

                for ( set=start[8]; set<=end[8]; set++ )
                {
                    srcInd[8] = set; dstInd[8] = set-start[8];

                    for ( rep=start[7]; rep<=end[7]; rep++ )
                    {
                        srcInd[7] = rep; dstInd[7] = rep-start[7];

                        for ( phs=start[6]; phs<=end[6]; phs++ )
                        {
                            srcInd[6] = phs; dstInd[6] = phs-start[6];

                            for ( con=start[5]; con<=end[5]; con++ )
                            {
                                srcInd[5] = con; dstInd[5] = con-start[5];

                                for ( s=start[4]; s<=end[4]; s++ )
                                {
                                    srcInd[4] = s; dstInd[4] = s-start[4];

                                    for ( n=start[3]; n<=end[3]; n++ )
                                    {
                                        srcInd[3] = n; dstInd[3] = n-start[3];

                                        for ( cha=start[2]; cha<=end[2]; cha++ )
                                        {
                                            srcInd[2] = cha; dstInd[2] = cha-start[2];

                                            for ( e1=start[1]; e1<=end[1]; e1++ )
                                            {
                                                srcInd[1] = e1; dstInd[1] = e1-start[1];

                                                srcInd[0] = start[0];
                                                dstInd[0] = 0;

                                                size_t offsetSrc = x.calculate_offset(srcInd);
                                                size_t offsetDst = r.calculate_offset(dstInd);

                                                memcpy(r.begin()+offsetDst, x.begin()+offsetSrc, sizeof(T)*(end[0]-start[0]+1));

                                                /*for ( ro=start[0]; ro<=end[0]; ro++ )
                                                {
                                                    srcInd[0] = ro;
                                                    dstInd[0] = ro-start[0];

                                                    int offsetSrc = x.calculate_offset(srcInd);
                                                    int offsetDst = r.calculate_offset(dstInd);

                                                    r(offsetDst) = x(offsetSrc);
                                                }*/
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        return true;
    }

    template <typename T> 
    bool setSubArrayUpTo11DArray(const hoNDArray<T>& x, hoNDArray<T>& r, const std::vector<size_t>& startND, std::vector<size_t>& size)
    {
        GADGET_CHECK_RETURN_FALSE( startND.size() == size.size() );
        GADGET_CHECK_RETURN_FALSE( startND.size() <= 11 );

        if ( r.get_number_of_elements() == x.get_number_of_elements() )
        {
            r = x;
            return true;
        }

        std::vector<size_t> start(11, 0);
        std::vector<size_t> end(11, 0);

        size_t ii;
        for ( ii=0; ii<startND.size(); ii++ )
        {
            start[ii] = startND[ii];
            end[ii] = start[ii] + size[ii] - 1;
            GADGET_CHECK_RETURN_FALSE(end[ii] < r.get_size(ii));
        }

        // [Ro E1 Cha Slice E2 Con Phase Rep Set Seg Ave]
        size_t e1, cha, n, s, con, phs, rep, set, seg, ave;

        std::vector<size_t> srcInd(11), dstInd(11);

        for ( ave=start[10]; ave<=end[10]; ave++ )
        {
            dstInd[10] = ave; srcInd[10] = ave-start[10];

            for ( seg=start[9]; seg<=end[9]; seg++ )
            {
                dstInd[9] = seg; srcInd[9] = seg-start[9];

                for ( set=start[8]; set<=end[8]; set++ )
                {
                    dstInd[8] = set; srcInd[8] = set-start[8];

                    for ( rep=start[7]; rep<=end[7]; rep++ )
                    {
                        dstInd[7] = rep; srcInd[7] = rep-start[7];

                        for ( phs=start[6]; phs<=end[6]; phs++ )
                        {
                            dstInd[6] = phs; srcInd[6] = phs-start[6];

                            for ( con=start[5]; con<=end[5]; con++ )
                            {
                                dstInd[5] = con; srcInd[5] = con-start[5];

                                for ( s=start[4]; s<=end[4]; s++ )
                                {
                                    dstInd[4] = s; srcInd[4] = s-start[4];

                                    for ( n=start[3]; n<=end[3]; n++ )
                                    {
                                        dstInd[3] = n; srcInd[3] = n-start[3];

                                        for ( cha=start[2]; cha<=end[2]; cha++ )
                                        {
                                            dstInd[2] = cha; srcInd[2] = cha-start[2];

                                            for ( e1=start[1]; e1<=end[1]; e1++ )
                                            {
                                                dstInd[1] = e1; srcInd[1] = e1-start[1];

                                                dstInd[0] = start[0];
                                                srcInd[0] = 0;

                                                size_t offsetSrc = x.calculate_offset(srcInd);
                                                size_t offsetDst = r.calculate_offset(dstInd);

                                                memcpy(r.begin()+offsetDst, x.begin()+offsetSrc, sizeof(T)*(end[0]-start[0]+1));
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        return true;
    }

    template<typename T> 
    bool extractSampledLinesUpTo11DArray(const hoNDArray<T>& x, hoNDArray<T>& r, const hoNDArray<float>& timeStamp, double acceFactorE1, double acceFactorE2)
    {
        try
        {
            std::vector<size_t> dim;
            x.get_dimensions(dim);

            size_t RO = x.get_size(0);
            size_t E1 = x.get_size(1);
            size_t CHA = x.get_size(2);
            size_t SLC = x.get_size(3);
            size_t E2 = x.get_size(4);
            size_t CON = x.get_size(5);
            size_t PHS = x.get_size(6);
            size_t REP = x.get_size(7);
            size_t SET = x.get_size(8);
            size_t SEG = x.get_size(9);
            size_t AVE = x.get_size(10);

            size_t Num = AVE*SEG*SET*REP*PHS*SLC;

            std::vector<size_t> dimRes(dim);

            if ( acceFactorE1>1 && E1>1 )
            {
                dimRes[1] = (size_t)(E1/acceFactorE1) + 1;
            }

            size_t dstE1 = dimRes[1];

            if ( acceFactorE2>1 && E2>1 )
            {
                dimRes[4] = (size_t)(E2/acceFactorE2) + 1;
            }

            r.create(&dimRes);
            Gadgetron::clear(r);

            // [Ro E1 Cha Slice E2 Con Phase Rep Set Seg Ave]

            size_t ROLen = sizeof(T)*RO;
            hoNDArray<T> dummyArray(SLC, CON, PHS, REP, SET, SEG, AVE);

            long long n;
            #pragma omp parallel default(none) private(n) shared(Num, dummyArray, RO, E1, CHA, SLC, E2, CON, PHS, REP, SET, SEG, AVE, timeStamp, x, r, ROLen, dstE1)
            {

                std::vector<size_t> indN;
                std::vector<size_t> srcInd(11, 0), dstInd(11, 0);
                size_t e1, cha, slc, e2, con, rep, phs, set, seg, ave;

                #pragma omp for
                for ( n=0; n<(long long)Num; n++ )
                {
                    indN = dummyArray.calculate_index(n);

                    ave = indN[6];
                    seg = indN[5];
                    set = indN[4];
                    rep = indN[3];
                    phs = indN[2];
                    con = indN[1];
                    slc = indN[0];

                    srcInd[10] = ave; dstInd[10] = ave;
                    srcInd[9] = seg; dstInd[9] = seg;
                    srcInd[8] = set; dstInd[8] = set;
                    srcInd[7] = rep; dstInd[7] = rep;
                    srcInd[6] = phs; dstInd[6] = phs;
                    srcInd[5] = con; dstInd[5] = con;
                    srcInd[4] = slc; dstInd[4] = slc;

                    // ------------------
                    size_t indE2(0);
                    size_t prevE2(0);
                    for ( e2=0; e2<E2; e2++ )
                    {
                        srcInd[3] = e2; dstInd[3] = indE2;

                        size_t indE1(0);
                        for ( e1=0; e1<E1; e1++ )
                        {
                            srcInd[1] = e1; 

                            srcInd[2] = 0;
                            if ( timeStamp(srcInd) > 0 )
                            {
                                dstInd[1] = indE1;
                                indE1++;

                                if ( e2 != prevE2 )
                                {
                                    prevE2 = e2;
                                    indE2++;
                                }

                                // ------------------
                                srcInd[2] = 0; dstInd[2] = 0;
                                size_t offsetSrc = x.calculate_offset(srcInd);
                                size_t offsetDst = r.calculate_offset(dstInd);

                                for ( cha=0; cha<CHA; cha++ )
                                {
                                    memcpy(r.begin()+offsetDst, x.begin()+offsetSrc, ROLen);

                                    offsetSrc += RO*E1;
                                    offsetDst += RO*dstE1;
                                }
                                // ------------------
                            }
                        }
                        // ------------------
                    }
                }
            }
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Errors in extractSampledLinesUpTo11DArray(...) ... ");
            return false;
        }

        return true;
    }

    template<typename T> 
    bool fillSampledLinesUpTo11DArray(const hoNDArray<T>& x, hoNDArray<T>& r, const hoNDArray<float>& timeStamp)
    {
        try
        {
            size_t RO = x.get_size(0);
            size_t E1 = timeStamp.get_size(1);
            size_t CHA = x.get_size(2);
            size_t SLC = timeStamp.get_size(3);
            size_t E2 = timeStamp.get_size(4);
            size_t CON = timeStamp.get_size(5);
            size_t PHS = timeStamp.get_size(6);
            size_t REP = timeStamp.get_size(7);
            size_t SET = timeStamp.get_size(8);
            size_t SEG = timeStamp.get_size(9);
            size_t AVE = timeStamp.get_size(10);

            size_t srcE1 = x.get_size(1);

            size_t Num = AVE*SEG*SET*REP*PHS*SLC;

            std::vector<size_t> dimRes;
            timeStamp.get_dimensions(dimRes);

            dimRes[0] = RO;
            dimRes[2] = CHA;
            r.create(&dimRes);
            Gadgetron::clear(r);

            size_t ROLen = sizeof(T)*RO;
            hoNDArray<T> dummyArray(SLC, CON, PHS, REP, SET, SEG, AVE);

            long long n;
            #pragma omp parallel default(none) private(n) shared(Num, dummyArray, RO, E1, CHA, SLC, E2, CON, PHS, REP, SET, SEG, AVE, timeStamp, x, r, ROLen, srcE1)
            {

                std::vector<size_t> indN;
                std::vector<size_t> srcInd(11, 0), dstInd(11, 0);
                size_t e1, cha, slc, e2, con, rep, phs, set, seg, ave;

                #pragma omp for
                for ( n=0; n<(long long)Num; n++ )
                {
                    indN = dummyArray.calculate_index(n);

                    ave = indN[6];
                    seg = indN[5];
                    set = indN[4];
                    rep = indN[3];
                    phs = indN[2];
                    con = indN[1];
                    slc = indN[0];

                    srcInd[10] = ave; dstInd[10] = ave;
                    srcInd[9] = seg; dstInd[9] = seg;
                    srcInd[8] = set; dstInd[8] = set;
                    srcInd[7] = rep; dstInd[7] = rep;
                    srcInd[6] = phs; dstInd[6] = phs;
                    srcInd[5] = con; dstInd[5] = con;
                    srcInd[4] = slc; dstInd[4] = slc;

                    // ------------------
                    size_t indE2(0);
                    size_t prevE2(0);
                    for ( e2=0; e2<E2; e2++ )
                    {
                        srcInd[3] = indE2; dstInd[3] = e2;

                        size_t indE1(0);
                        for ( e1=0; e1<E1; e1++ )
                        {
                            dstInd[1] = e1; 

                            dstInd[2] = 0;
                            if ( timeStamp(dstInd) > 0 )
                            {
                                srcInd[1] = indE1;
                                indE1++;

                                if ( e2 != prevE2 )
                                {
                                    prevE2 = e2;
                                    indE2++;
                                }

                                // ------------------

                                srcInd[2] = 0; dstInd[2] = 0;
                                size_t offsetSrc = x.calculate_offset(srcInd);
                                size_t offsetDst = r.calculate_offset(dstInd);

                                for ( cha=0; cha<CHA; cha++ )
                                {
                                    memcpy(r.begin()+offsetDst, x.begin()+offsetSrc, ROLen);

                                    offsetSrc += RO*srcE1;
                                    offsetDst += RO*E1;
                                }
                                // ------------------
                            }
                        }
                        // ------------------
                    }
                }
            }
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Errors in fillSampledLinesUpTo11DArray(...) ... ");
            return false;
        }

        return true;
    }

    template<typename T> 
    bool stdOver3rdDimension(const hoNDArray<T>& x, hoNDArray<T>& std, bool NMinusOne)
    {
        try
        {
            typedef typename realType<T>::Type value_type;

            GADGET_DEBUG_CHECK_RETURN_FALSE(x.get_number_of_dimensions() >= 3);

            size_t RO = x.get_size(0);
            size_t E1 = x.get_size(1);
            size_t CHA = x.get_size(2);

            long long num = (long long)x.get_number_of_elements() / (RO*E1*CHA);

            boost::shared_ptr< std::vector<size_t> > dim = x.get_dimensions();

            std::vector<size_t> dimStd(*dim);
            dimStd.erase(dimStd.begin()+2);
            std.create(&dimStd);

            std::vector<size_t> dim3D(3);
            dim3D[0] = RO;
            dim3D[1] = E1;
            dim3D[2] = CHA;

            T S( (value_type)CHA );
            if ( NMinusOne )
            {
                S = T( (value_type)CHA-1);
            }

            T v(0), v1(0);
            T S2 = T( (value_type)1.0 )/S;
            T S3 = T( (value_type)1.0 )/T( (value_type)CHA );

            long long n;

            #ifdef GCC_OLD_FLAG
                #pragma omp parallel for default(none) private(n) shared(num, RO, E1, CHA, S, S2, S3, v, v1)
            #else
                #pragma omp parallel for default(none) private(n) shared(num, RO, E1, CHA, x, std, S, S2, S3, v, v1)
            #endif
            for ( n=0; n<num; n++ )
            {
                hoNDArray<T> xTmp(RO, E1, CHA, const_cast<T*>(x.begin()+n*RO*E1*CHA));
                hoNDArray<T> mean(RO, E1);
                Gadgetron::clear(mean);

                size_t ro, e1, cha;
                for ( cha=0; cha<CHA; cha++ )
                {
                    for ( e1=0; e1<E1; e1++ )
                    {
                        for ( ro=0; ro<RO; ro++ )
                        {
                            mean(ro+e1*RO) += xTmp(cha*RO*E1+e1*RO+ro)*S3;
                        }
                    }
                }

                for ( e1=0; e1<E1; e1++ )
                {
                    for ( ro=0; ro<RO; ro++ )
                    {
                        size_t ind = e1*RO+ro;

                        v = 0; v1 = 0;
                        for ( cha=0; cha<CHA; cha++ )
                        {
                            v1 = std::abs(xTmp(cha*RO*E1+ind)-mean(ind));
                            v += v1*v1;
                        }

                        v /= S;
                        std(ind+n*RO*E1) = std::sqrt(v);
                    }
                }
            }
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Errors in stdOver3rdDimension(const hoNDArray<T>& x, hoNDArray<T>& std, bool NMinusOne) ... ");
            return false;
        }

        return true;
    }

    template<typename T> 
    bool cropOver3rdDimension(const hoNDArray<T>& x, hoNDArray<T>& r, size_t start, size_t end)
    {
        try
        {
            boost::shared_ptr< std::vector<size_t> > dimX = x.get_dimensions();

            size_t NDim = dimX->size();

            if ( NDim <= 2 )
            {
                r = x;
                return true;
            }

            size_t RO = x.get_size(0);
            size_t E1 = x.get_size(1);
            size_t E2 = x.get_size(2);

            size_t E2_R = end-start+1;

            if ( E2 <= E2_R )
            {
                r = x;
                return true;
            }

            std::vector<size_t> dimR(*dimX);
            dimR[2] = E2_R;

            r.create(&dimR);

            size_t N2D = RO*E1;
            size_t N3D = RO*E1*E2;
            size_t N3D_R = RO*E1*E2_R;

            size_t N = x.get_number_of_elements()/N3D;

            const T* pX = x.begin();
            T* pR = r.begin();

            size_t n;
            for ( n=0; n<N; n++ )
            {
                long long e2;
                #pragma omp parallel for default(none) private(e2) shared(N2D, N3D, N3D_R, pX, pR, RO, E1, E2, n, start, end)
                for ( e2=(long long)start; e2<=(long long)end; e2++ )
                {
                    memcpy(pR+n*N3D_R+(e2-start)*N2D, pX+n*N3D+e2*N2D, sizeof(T)*N2D);
                }
            }
        }
        catch (...)
        {
            GADGET_ERROR_MSG("Errors in cropOver3rdDimension(const hoNDArray<T>& x, hoNDArray<T>& r, size_t start, size_t end) ... ");
            return false;
        }
        return true;
    }

    template<typename T> bool setSubArrayOver3rdDimension(const hoNDArray<T>& x, hoNDArray<T>& r, size_t start, size_t end)
    {
        try
        {
            boost::shared_ptr< std::vector<size_t> > dimR = r.get_dimensions();

            size_t NDim = dimR->size();

            if ( NDim <= 2 )
            {
                r = x;
                return true;
            }

            size_t RO = r.get_size(0);
            size_t E1 = r.get_size(1);
            size_t E2 = r.get_size(2);

            size_t E2_X = end-start+1;
            GADGET_CHECK_RETURN_FALSE( E2_X == x.get_size(2) );

            if ( E2_X >= E2 )
            {
                r = x;
                return true;
            }

            size_t N2D = RO*E1;
            size_t N3D = RO*E1*E2;
            size_t N3D_X = RO*E1*E2_X;

            size_t N = r.get_number_of_elements()/N3D;

            const T* pX = x.begin();
            T* pR = r.begin();

            size_t n;
            for ( n=0; n<N; n++ )
            {
                long long e2;
                #pragma omp parallel for default(none) private(e2) shared(N2D, N3D, N3D_X, pX, pR, RO, E1, E2, n, start, end)
                for ( e2=(long long)start; e2<=(long long)end; e2++ )
                {
                    memcpy(pR+n*N3D+e2*N2D, pX+n*N3D_X+(e2-start)*N2D, sizeof(T)*N2D);
                }
            }
        }
        catch (...)
        {
            GADGET_ERROR_MSG("Errors in setSubArrayOver3rdDimension(const hoNDArray<T>& x, hoNDArray<T>& r, size_t start, size_t end) ... ");
            return false;
        }
        return true;
    }

    template<typename T> 
    bool permuteE2To3rdDimension(const hoNDArray<T>& x, hoNDArray<T>& r)
    {
        try
        {
            boost::shared_ptr< std::vector<size_t> > dimX = x.get_dimensions();

            size_t NDim = dimX->size();

            if ( NDim <= 5 )
            {
                r = x;
                return true;
            }

            size_t RO = x.get_size(0);
            size_t E1 = x.get_size(1);
            size_t CHA = x.get_size(2);
            size_t SLC = x.get_size(3);
            size_t E2 = x.get_size(4);

            std::vector<size_t> dimR(*dimX);
            dimR[2] = E2;
            dimR[3] = CHA;
            dimR[4] = SLC;

            r.create(&dimR);

            size_t N2D = RO*E1;
            size_t N5D = RO*E1*CHA*E2*SLC;

            size_t N = x.get_number_of_elements()/N5D;

            const T* pX = x.begin();
            T* pR = r.begin();

            size_t n;
            for ( n=0; n<N; n++ )
            {
                int e2;
                #pragma omp parallel for default(none) private(e2) shared(N5D, N2D, pX, pR, CHA, SLC, E2, n)
                for ( e2=0; e2<E2; e2++ )
                {
                    for ( size_t slc=0; slc<SLC; slc++ )
                    {
                        for ( size_t cha=0; cha<CHA; cha++ )
                        {
                            memcpy(pR+n*N5D+slc*CHA*E2*N2D+cha*E2*N2D+e2*N2D, pX+n*N5D+e2*SLC*CHA*N2D+slc*CHA*N2D+cha*N2D, sizeof(T)*N2D);
                        }
                    }
                }
            }
        }
        catch (...)
        {
            GADGET_ERROR_MSG("Errors in permuteE2To3rdDimension(const hoNDArray<T>& x, hoNDArray<T>& r) ... ");
            return false;
        }
        return true;
    }

    template<typename T> 
    bool permuteE2To5thDimension(const hoNDArray<T>& x, hoNDArray<T>& r)
    {
        try
        {
            boost::shared_ptr< std::vector<size_t> > dimX = x.get_dimensions();

            size_t NDim = dimX->size();

            if ( NDim < 5 )
            {
                r = x;
                return true;
            }

            size_t RO = x.get_size(0);
            size_t E1 = x.get_size(1);
            size_t E2 = x.get_size(2);
            size_t CHA = x.get_size(3);
            size_t SLC = x.get_size(4);

            std::vector<size_t> dimR(*dimX);
            dimR[2] = CHA;
            dimR[3] = SLC;
            dimR[4] = E2;

            r.create(&dimR);

            size_t N2D = RO*E1;
            size_t N5D = RO*E1*CHA*E2*SLC;

            size_t N = x.get_number_of_elements()/N5D;

            const T* pX = x.begin();
            T* pR = r.begin();

            size_t n;
            for ( n=0; n<N; n++ )
            {
                int e2;
                #pragma omp parallel for default(none) private(e2) shared(N5D, N2D, pX, pR, CHA, SLC, E2, n)
                for ( e2=0; e2<E2; e2++ )
                {
                    for ( size_t slc=0; slc<SLC; slc++ )
                    {
                        for ( size_t cha=0; cha<CHA; cha++ )
                        {
                            memcpy(pR+n*N5D+e2*SLC*CHA*N2D+slc*CHA*N2D+cha*N2D, pX+n*N5D+slc*CHA*E2*N2D+cha*E2*N2D+e2*N2D, sizeof(T)*N2D);
                        }
                    }
                }
            }
        }
        catch (...)
        {
            GADGET_ERROR_MSG("Errors in permuteE2To5thDimension(const hoNDArray<T>& x, hoNDArray<T>& r) ... ");
            return false;
        }
        return true;
    }

    template<typename T> 
    bool permuteROTo3rdDimensionFor3DRecon(const hoNDArray<T>& x, hoNDArray<T>& r)
    {
        try
        {
            boost::shared_ptr< std::vector<size_t> > dimX = x.get_dimensions();

            size_t NDim = dimX->size();

            if ( NDim < 3 )
            {
                r = x;
                return true;
            }

            size_t RO = x.get_size(0);
            size_t E1 = x.get_size(1);
            size_t E2 = x.get_size(2);

            std::vector<size_t> dimR(*dimX);
            dimR[0] = E1;
            dimR[1] = E2;
            dimR[2] = RO;

            r.create(&dimR);

            size_t N3D = RO*E1*E2;

            size_t N = x.get_number_of_elements()/N3D;

            const T* pX = x.begin();
            T* pR = r.begin();

            long long n;

            #pragma omp parallel for default(none) private(n) shared(RO, E1, E2, N, pR, N3D, pX)
            for ( n=0; n<(long long)N; n++ )
            {
                T* pRn = pR + n*N3D;
                T* pXn = const_cast<T*>(pX) + n*N3D;

                for ( size_t e2=0; e2<E2; e2++ )
                {
                    for ( size_t e1=0; e1<E1; e1++ )
                    {
                        for ( size_t ro=0; ro<RO; ro++ )
                        {
                            pRn[e1+e2*E1+ro*E1*E2] = pXn[ro+e1*RO+e2*RO*E1];
                        }
                    }
                }
            }
        }
        catch (...)
        {
            GADGET_ERROR_MSG("Errors in permuteROTo3rdDimensionFor3DRecon(const hoNDArray<T>& x, hoNDArray<T>& r) ... ");
            return false;
        }
        return true;
    }

    template<typename T> 
    bool permuteROTo4thDimensionFor3DRecon(const hoNDArray<T>& x, hoNDArray<T>& r)
    {
        try
        {
            boost::shared_ptr< std::vector<size_t> > dimX = x.get_dimensions();

            size_t NDim = dimX->size();

            if ( NDim < 4 )
            {
                r = x;
                return true;
            }

            size_t RO = x.get_size(0);
            size_t E1 = x.get_size(1);
            size_t E2 = x.get_size(2);
            size_t CHA = x.get_size(3);

            std::vector<size_t> dimR(*dimX);
            dimR[0] = E1;
            dimR[1] = E2;
            dimR[2] = CHA;
            dimR[3] = RO;

            r.create(&dimR);

            size_t N4D = RO*E1*E2*CHA;

            size_t N = x.get_number_of_elements()/N4D;

            const T* pX = x.begin();
            T* pR = r.begin();

            long long n;
            for ( n=0; n<(long long)N; n++ )
            {
                T* pRn = pR + n*N4D;
                T* pXn = const_cast<T*>(pX) + n*N4D;

                long long cha;

                #pragma omp parallel for default(none) private(cha) shared(RO, E1, E2, CHA, pXn, pRn)
                for ( cha=0; cha<(long long)CHA; cha++ )
                {
                    for ( size_t e2=0; e2<E2; e2++ )
                    {
                        for ( size_t e1=0; e1<E1; e1++ )
                        {
                            for ( size_t ro=0; ro<RO; ro++ )
                            {
                                pRn[e1+e2*E1+cha*E1*E2+ro*E1*E2*CHA] = pXn[ro+e1*RO+e2*RO*E1+cha*RO*E1*E2];
                            }
                        }
                    }
                }
            }
        }
        catch (...)
        {
            GADGET_ERROR_MSG("Errors in permuteROTo4thDimensionFor3DRecon(const hoNDArray<T>& x, hoNDArray<T>& r) ... ");
            return false;
        }
        return true;
    }

    template<typename T> 
    bool permuteROTo1stDimensionFor3DRecon(const hoNDArray<T>& x, hoNDArray<T>& r)
    {
        try
        {
            boost::shared_ptr< std::vector<size_t> > dimX = x.get_dimensions();

            size_t NDim = dimX->size();

            if ( NDim < 4 )
            {
                r = x;
                return true;
            }

            size_t E1 = x.get_size(0);
            size_t E2 = x.get_size(1);
            size_t CHA = x.get_size(2);
            size_t RO = x.get_size(3);

            std::vector<size_t> dimR(*dimX);
            dimR[0] = RO;
            dimR[1] = E1;
            dimR[2] = E2;
            dimR[3] = CHA;

            r.create(&dimR);

            size_t N4D = RO*E1*E2*CHA;

            size_t N = x.get_number_of_elements()/N4D;

            const T* pX = x.begin();
            T* pR = r.begin();

            long long n;
            for ( n=0; n<(long long)N; n++ )
            {
                T* pRn = pR + n*N4D;
                T* pXn = const_cast<T*>(pX) + n*N4D;

                long long cha;

                #pragma omp parallel for default(none) private(cha) shared(RO, E1, E2, CHA, pXn, pRn)
                for ( cha=0; cha<(long long)CHA; cha++ )
                {
                    for ( size_t e2=0; e2<E2; e2++ )
                    {
                        for ( size_t e1=0; e1<E1; e1++ )
                        {
                            size_t indRn = e1*RO+e2*RO*E1+cha*RO*E1*E2;
                            size_t indXn = e1+e2*E1+cha*E1*E2;
                            for ( size_t ro=0; ro<RO; ro++ )
                            {
                                pRn[ro+indRn] = pXn[indXn+ro*E1*E2*CHA];
                            }
                        }
                    }
                }
            }
        }
        catch (...)
        {
            GADGET_ERROR_MSG("Errors in permuteROTo1stDimensionFor3DRecon(const hoNDArray<T>& x, hoNDArray<T>& r) ... ");
            return false;
        }
        return true;
    }

    template<typename T> 
    bool permute3rdDimensionTo1stDimension(const hoNDArray<T>& x, hoNDArray<T>& r)
    {
        try
        {
            boost::shared_ptr< std::vector<size_t> > dimX = x.get_dimensions();

            size_t NDim = dimX->size();

            if ( NDim < 3 )
            {
                r = x;
                return true;
            }

            size_t RO = x.get_size(0);
            size_t E1 = x.get_size(1);
            size_t E2 = x.get_size(2);

            std::vector<size_t> dimR(*dimX);
            dimR[0] = E2;
            dimR[1] = RO;
            dimR[2] = E1;

            r.create(&dimR);

            size_t N3D = RO*E1*E2;

            size_t N = x.get_number_of_elements()/N3D;

            const T* pX = x.begin();
            T* pR = r.begin();

            long long n, e2;
            for ( n=0; n<(long long)N; n++ )
            {
                T* pRn = pR + n*N3D;
                T* pXn = const_cast<T*>(pX) + n*N3D;

                #pragma omp parallel for default(none) private(e2) shared(RO, E1, E2, pXn, pRn)
                for ( e2=0; e2<(long long)E2; e2++ )
                {
                    for ( size_t e1=0; e1<E1; e1++ )
                    {
                        size_t indRn = e2+e1*E2*RO;
                        size_t indXn = e1*RO+e2*RO*E1;
                        for ( size_t ro=0; ro<RO; ro++ )
                        {
                            pRn[ro*E2+indRn] = pXn[ro+indXn];
                        }
                    }
                }
            }
        }
        catch (...)
        {
            GADGET_ERROR_MSG("Errors in permute3rdDimensionTo1stDimension(const hoNDArray<T>& x, hoNDArray<T>& r) ... ");
            return false;
        }
        return true;
    }

    template<typename T> 
    bool permuteROTo5thDimensionFor3DRecon(const hoNDArray<T>& x, hoNDArray<T>& r)
    {
        try
        {
            boost::shared_ptr< std::vector<size_t> > dimX = x.get_dimensions();

            size_t NDim = dimX->size();

            if ( NDim < 5 )
            {
                r = x;
                return true;
            }

            size_t RO = x.get_size(0);
            size_t E1 = x.get_size(1);
            size_t E2 = x.get_size(2);
            size_t srcCHA = x.get_size(3);
            size_t dstCHA = x.get_size(4);

            std::vector<size_t> dimR(*dimX);
            dimR[0] = E1;
            dimR[1] = E2;
            dimR[2] = srcCHA;
            dimR[3] = dstCHA;
            dimR[4] = RO;

            r.create(&dimR);

            size_t N5D = RO*E1*E2*srcCHA*dstCHA;

            size_t N = x.get_number_of_elements()/N5D;

            const T* pX = x.begin();
            T* pR = r.begin();

            long long n;
            for ( n=0; n<(long long)N; n++ )
            {
                T* pRn = pR + n*N5D;
                T* pXn = const_cast<T*>(pX) + n*N5D;

                long long dcha;

                #pragma omp parallel for default(none) private(dcha) shared(RO, E1, E2, srcCHA, dstCHA, pXn, pRn)
                for ( dcha=0; dcha<(long long)dstCHA; dcha++ )
                {
                    for ( size_t scha=0; scha<(int)srcCHA; scha++ )
                    {
                        for ( size_t e2=0; e2<E2; e2++ )
                        {
                            for ( size_t e1=0; e1<E1; e1++ )
                            {
                                size_t indRn = e1+e2*E1+scha*E1*E2+dcha*E1*E2*srcCHA;
                                size_t indXn = e1*RO+e2*RO*E1+scha*RO*E1*E2+dcha*RO*E1*E2*srcCHA;
                                for ( size_t ro=0; ro<RO; ro++ )
                                {
                                    pRn[indRn+ro*E1*E2*srcCHA*dstCHA] = pXn[ro+indXn];
                                }
                            }
                        }
                    }
                }
            }
        }
        catch (...)
        {
            GADGET_ERROR_MSG("Errors in permuteROTo5thDimensionFor3DRecon(const hoNDArray<T>& x, hoNDArray<T>& r) ... ");
            return false;
        }
        return true;
    }

    template<typename T> 
    bool imageDomainUnwrapping2D(const hoNDArray<T>& x, const hoNDArray<T>& kernel, hoNDArray<T>& buf, hoNDArray<T>& y)
    {
        try
        {
            T* pX = const_cast<T*>(x.begin());
            T* ker = const_cast<T*>(kernel.begin());
            T* pY = y.begin();

            size_t ro = x.get_size(0);
            size_t e1 = x.get_size(1);
            size_t srcCHA = x.get_size(2);
            size_t dstCHA = kernel.get_size(3);

            if ( buf.get_number_of_elements() < ro*e1*srcCHA )
            {
                buf.create(ro, e1, srcCHA);
            }
            T* pBuf = buf.begin();

            size_t dCha;

            //#pragma omp parallel default(shared)
            {
                //#ifdef WIN32
                //    int tid = omp_get_thread_num();
                //    DWORD_PTR mask = (1 << tid);
                //    // GADGET_MSG("thread id : " << tid << " - mask : " << mask);
                //    SetThreadAffinityMask( GetCurrentThread(), mask );
                //#endif // WIN32

                //#pragma omp for

                if ( typeid(T)==typeid(GT_Complex8) )
                {
                    for ( dCha=0; dCha<dstCHA; dCha++ )
                    {
                        Gadgetron::math::multiply(ro*e1*srcCHA, pX, ker+dCha*ro*e1*srcCHA, pBuf);

                        memcpy(pY+dCha*ro*e1, pBuf, sizeof(T)*ro*e1);
                        for ( size_t sCha=1; sCha<srcCHA; sCha++ )
                        {
                            Gadgetron::math::add(ro*e1, pY+dCha*ro*e1, pBuf+sCha*ro*e1, pY+dCha*ro*e1);
                        }
                    }
                }
                else if ( typeid(T)==typeid(GT_Complex16) )
                {
                    for ( dCha=0; dCha<dstCHA; dCha++ )
                    {
                        Gadgetron::math::multiply(ro*e1*srcCHA, pX, ker+dCha*ro*e1*srcCHA, pBuf);

                        memcpy(pY+dCha*ro*e1, pBuf, sizeof(T)*ro*e1);
                        for ( size_t sCha=1; sCha<srcCHA; sCha++ )
                        {
                            Gadgetron::math::add(ro*e1, pY+dCha*ro*e1, pBuf+sCha*ro*e1, pY+dCha*ro*e1);
                        }
                    }
                }
            }
        }
        catch (...)
        {
            GADGET_ERROR_MSG("Errors in imageDomainUnwrapping2D(const hoNDArray<T>& x, const hoNDArray<T>& ker, hoNDArray<T>& buf, hoNDArray<T>& y) ... ");
            return false;
        }
        return true;
    }

    template<typename T> 
    bool imageDomainUnwrapping2DT(const hoNDArray<T>& x, const hoNDArray<T>& kernel, hoNDArray<T>& buf, hoNDArray<T>& y)
    {
        try
        {
            long long ro = (long long)x.get_size(0);
            long long e1 = (long long)x.get_size(1);
            long long srcCHA = (long long)x.get_size(2);
            long long N = (long long)x.get_size(3);

            long long dstCHA = (long long)kernel.get_size(3);
            long long kerN = (long long)kernel.get_size(4);

            if ( (long long)buf.get_number_of_elements() < ro*e1*srcCHA )
            {
                buf.create(ro, e1, srcCHA);
            }
            T* pBuf = buf.begin();

            long long n, dCha;

            //#pragma omp parallel default(shared)
            {
                //#ifdef WIN32
                //    int tid = omp_get_thread_num();
                //    DWORD_PTR mask = (1 << tid);
                //    // GADGET_MSG("thread id : " << tid << " - mask : " << mask);
                //    SetThreadAffinityMask( GetCurrentThread(), mask );
                //#endif // WIN32

                //#pragma omp for

                if ( typeid(T)==typeid(GT_Complex8) )
                {
                    const T* pXN = x.begin();
                    T* pYN = y.begin();
                    T* pBufN = buf.begin();
                    const T* pKerN = kernel.begin();

                    // #pragma omp parallel for default(none) private(dCha, n) shared(N, ro, e1, srcCHA, dstCHA, kerN, pXN, pYN, pBufN, pKerN)
                    for ( dCha=0; dCha<(long long)dstCHA; dCha++ )
                    {
                        for ( n=0; n<N; n++  )
                        {
                            const T* ker = pKerN + n*ro*e1*srcCHA*dstCHA;
                            if ( kerN <= n )
                            {
                                ker = pKerN + (kerN-1)*ro*e1*srcCHA*dstCHA;
                            }

                            const T* pX = pXN + n*ro*e1*srcCHA;
                            T* pBuf =pBufN + n*ro*e1*srcCHA;

                            Gadgetron::math::multiply(ro*e1*srcCHA, pX, ker+dCha*ro*e1*srcCHA, pBuf);
                        //}

                        //for ( n=0; n<N; n++  )
                        //{
                            T* pY = pYN + n*ro*e1*dstCHA;
                            //T* pBuf =pBufN + n*ro*e1*srcCHA;

                            memcpy(pY+dCha*ro*e1, pBuf, sizeof(T)*ro*e1);
                            for ( long long sCha=1; sCha<srcCHA; sCha++ )
                            {
                                Gadgetron::math::add(ro*e1, pY+dCha*ro*e1, pBuf+sCha*ro*e1, pY+dCha*ro*e1);
                            }
                        }
                    }
                }
                else if ( typeid(T)==typeid(GT_Complex16) )
                {
                    for ( n=0; n<N; n++ )
                    {
                        const T* ker = kernel.begin() + n*ro*e1*srcCHA*dstCHA;
                        if ( kerN <= n )
                        {
                            ker = kernel.begin() + (kerN-1)*ro*e1*srcCHA*dstCHA;
                        }

                        const T* pX = x.begin() + n*ro*e1*srcCHA;
                        T* pY = y.begin() + n*ro*e1*dstCHA;

                        for ( long long dCha=0; dCha<dstCHA; dCha++ )
                        {
                            Gadgetron::math::multiply(ro*e1*srcCHA, pX, ker+dCha*ro*e1*srcCHA, pBuf);

                            memcpy(pY+dCha*ro*e1, pBuf, sizeof(T)*ro*e1);
                            for ( long long sCha=1; sCha<srcCHA; sCha++ )
                            {
                                Gadgetron::math::add(ro*e1, pY+dCha*ro*e1, pBuf+sCha*ro*e1, pY+dCha*ro*e1);
                            }
                        }
                    }
                }
            }
        }
        catch (...)
        {
            GADGET_ERROR_MSG("Errors in imageDomainUnwrapping2DT(const hoNDArray<T>& x, const hoNDArray<T>& ker, hoNDArray<T>& buf, hoNDArray<T>& y) ... ");
            return false;
        }
        return true;
    }

    template<typename CoordType, typename T> 
    bool computePeriodicBoundaryValues(const hoNDArray<CoordType>& x, const hoNDArray<T>& y, CoordType start, CoordType end, hoNDArray<CoordType>& vx, hoNDArray<T>& vy)
    {
        try
        {
            typedef typename realType<T>::Type real_value_type;

            size_t N = x.get_size(0);
            size_t M = y.get_size(1);

            GADGET_CHECK_RETURN_FALSE(y.get_size(0)==N);
            GADGET_CHECK_RETURN_FALSE(start<=x(0));
            GADGET_CHECK_RETURN_FALSE(end>=x(N-1));

            vx.create(N+2);
            vy.create(N+2, M);

            size_t m, n;

            vx(0) = start;
            for ( n=0; n<N; n++ )
            {
                vx(n+1) = x(n);
            }
            vx(N+1) = end;

            CoordType dS = x(0) - start;
            CoordType dE = end - x(N-1);

            // start, end
            CoordType wS;
            if ( dE+dS > FLT_EPSILON )
                wS = dE/(dE+dS);
            else
                wS = dE/FLT_EPSILON;

            for ( m=0; m<M; m++ )
            {
                T a = y(0, m);
                T b = y(N-1, m);

                vy(0, m) = b + (real_value_type)wS * ( a - b );
                vy(N+1, m) = vy(0, m);
            }

            // middle
            for ( n=0; n<N; n++ )
            {
                for ( m=0; m<M; m++ )
                {
                    vy(n+1, m) = y(n, m);
                }
            }
        }
        catch (...)
        {
            GADGET_ERROR_MSG("Errors in computePeriodicBoundaryValues(const hoNDArray<CoordType>& x, const hoNDArray<T>& y, CoordType& start, CoordType& end, hoNDArray<T>& r) ... ");
            return false;
        }
        return true;
    }

    template EXPORTGTPLUS bool sumOverLastDimension(const hoNDArray<float>& x, hoNDArray<float>& r);
    template EXPORTGTPLUS bool sumOverLastDimension(const hoNDArray<double>& x, hoNDArray<double>& r);
    template EXPORTGTPLUS bool sumOverLastDimension(const hoNDArray<GT_Complex8>& x, hoNDArray<GT_Complex8>& r);
    template EXPORTGTPLUS bool sumOverLastDimension(const hoNDArray<GT_Complex16>& x, hoNDArray<GT_Complex16>& r);

    template EXPORTGTPLUS bool sumOverSecondLastDimension(const hoNDArray<float>& x, hoNDArray<float>& r);
    template EXPORTGTPLUS bool sumOverSecondLastDimension(const hoNDArray<double>& x, hoNDArray<double>& r);
    template EXPORTGTPLUS bool sumOverSecondLastDimension(const hoNDArray<GT_Complex8>& x, hoNDArray<GT_Complex8>& r);
    template EXPORTGTPLUS bool sumOverSecondLastDimension(const hoNDArray<GT_Complex16>& x, hoNDArray<GT_Complex16>& r);

    template EXPORTGTPLUS bool multiplyOverLastDimension(const hoNDArray<float>& x, const hoNDArray<float>& y, hoNDArray<float>& r);
    template EXPORTGTPLUS bool multiplyOverLastDimension(const hoNDArray<double>& x, const hoNDArray<double>& y, hoNDArray<double>& r);
    template EXPORTGTPLUS bool multiplyOverLastDimension(const hoNDArray<GT_Complex8>& x, const hoNDArray<GT_Complex8>& y, hoNDArray<GT_Complex8>& r);
    template EXPORTGTPLUS bool multiplyOverLastDimension(const hoNDArray<GT_Complex16>& x, const hoNDArray<GT_Complex16>& y, hoNDArray<GT_Complex16>& r);

    template EXPORTGTPLUS bool divideOverLastDimension(const hoNDArray<float>& x, const hoNDArray<float>& y, hoNDArray<float>& r);
    template EXPORTGTPLUS bool divideOverLastDimension(const hoNDArray<double>& x, const hoNDArray<double>& y, hoNDArray<double>& r);
    template EXPORTGTPLUS bool divideOverLastDimension(const hoNDArray<GT_Complex8>& x, const hoNDArray<GT_Complex8>& y, hoNDArray<GT_Complex8>& r);
    template EXPORTGTPLUS bool divideOverLastDimension(const hoNDArray<GT_Complex16>& x, const hoNDArray<GT_Complex16>& y, hoNDArray<GT_Complex16>& r);

    template EXPORTGTPLUS bool sumOver1stDimension(const hoNDArray<float>& x, hoNDArray<float>& r);
    template EXPORTGTPLUS bool sumOver1stDimension(const hoNDArray<double>& x, hoNDArray<double>& r);
    template EXPORTGTPLUS bool sumOver1stDimension(const hoNDArray<GT_Complex8>& x, hoNDArray<GT_Complex8>& r);
    template EXPORTGTPLUS bool sumOver1stDimension(const hoNDArray<GT_Complex16>& x, hoNDArray<GT_Complex16>& r);

    template EXPORTGTPLUS bool sumOver2ndDimension(const hoNDArray<float>& x, hoNDArray<float>& r);
    template EXPORTGTPLUS bool sumOver2ndDimension(const hoNDArray<double>& x, hoNDArray<double>& r);
    template EXPORTGTPLUS bool sumOver2ndDimension(const hoNDArray<GT_Complex8>& x, hoNDArray<GT_Complex8>& r);
    template EXPORTGTPLUS bool sumOver2ndDimension(const hoNDArray<GT_Complex16>& x, hoNDArray<GT_Complex16>& r);

    template EXPORTGTPLUS bool sumOver3rdDimension(const hoNDArray<float>& x, hoNDArray<float>& r);
    template EXPORTGTPLUS bool sumOver3rdDimension(const hoNDArray<double>& x, hoNDArray<double>& r);
    template EXPORTGTPLUS bool sumOver3rdDimension(const hoNDArray<GT_Complex8>& x, hoNDArray<GT_Complex8>& r);
    template EXPORTGTPLUS bool sumOver3rdDimension(const hoNDArray<GT_Complex16>& x, hoNDArray<GT_Complex16>& r);

    template EXPORTGTPLUS bool sumOver4thDimension(const hoNDArray<float>& x, hoNDArray<float>& r);
    template EXPORTGTPLUS bool sumOver4thDimension(const hoNDArray<double>& x, hoNDArray<double>& r);
    template EXPORTGTPLUS bool sumOver4thDimension(const hoNDArray<GT_Complex8>& x, hoNDArray<GT_Complex8>& r);
    template EXPORTGTPLUS bool sumOver4thDimension(const hoNDArray<GT_Complex16>& x, hoNDArray<GT_Complex16>& r);

    template EXPORTGTPLUS bool sumOver5thDimension(const hoNDArray<float>& x, hoNDArray<float>& r);
    template EXPORTGTPLUS bool sumOver5thDimension(const hoNDArray<double>& x, hoNDArray<double>& r);
    template EXPORTGTPLUS bool sumOver5thDimension(const hoNDArray<GT_Complex8>& x, hoNDArray<GT_Complex8>& r);
    template EXPORTGTPLUS bool sumOver5thDimension(const hoNDArray<GT_Complex16>& x, hoNDArray<GT_Complex16>& r);

    template EXPORTGTPLUS bool multiplyOver3rdDimension(const hoNDArray<float>& x3D, const hoNDArray<float>& y4D, hoNDArray<float>& r);
    template EXPORTGTPLUS bool multiplyOver3rdDimension(const hoNDArray<double>& x3D, const hoNDArray<double>& y4D, hoNDArray<double>& r);
    template EXPORTGTPLUS bool multiplyOver3rdDimension(const hoNDArray<GT_Complex8>& x3D, const hoNDArray<GT_Complex8>& y4D, hoNDArray<GT_Complex8>& r);
    template EXPORTGTPLUS bool multiplyOver3rdDimension(const hoNDArray<GT_Complex16>& x3D, const hoNDArray<GT_Complex16>& y4D, hoNDArray<GT_Complex16>& r);

    template EXPORTGTPLUS bool multiplyOver4thDimension(const hoNDArray<float>& x4D, const hoNDArray<float>& y5D, hoNDArray<float>& r);
    template EXPORTGTPLUS bool multiplyOver4thDimension(const hoNDArray<double>& x4D, const hoNDArray<double>& y5D, hoNDArray<double>& r);
    template EXPORTGTPLUS bool multiplyOver4thDimension(const hoNDArray<GT_Complex8>& x4D, const hoNDArray<GT_Complex8>& y5D, hoNDArray<GT_Complex8>& r);
    template EXPORTGTPLUS bool multiplyOver4thDimension(const hoNDArray<GT_Complex16>& x4D, const hoNDArray<GT_Complex16>& y5D, hoNDArray<GT_Complex16>& r);

    template EXPORTGTPLUS bool multiplyOver4thDimensionExcept(const hoNDArray<float>& x4D, const hoNDArray<float>& y5D, size_t n, hoNDArray<float>& r, bool copyY2R);
    template EXPORTGTPLUS bool multiplyOver4thDimensionExcept(const hoNDArray<double>& x4D, const hoNDArray<double>& y5D, size_t n, hoNDArray<double>& r, bool copyY2R);
    template EXPORTGTPLUS bool multiplyOver4thDimensionExcept(const hoNDArray<GT_Complex8>& x4D, const hoNDArray<GT_Complex8>& y5D, size_t n, hoNDArray<GT_Complex8>& r, bool copyY2R);
    template EXPORTGTPLUS bool multiplyOver4thDimensionExcept(const hoNDArray<GT_Complex16>& x4D, const hoNDArray<GT_Complex16>& y5D, size_t n, hoNDArray<GT_Complex16>& r, bool copyY2R);

    template EXPORTGTPLUS bool multiplyOver5thDimension(const hoNDArray<float>& x, const hoNDArray<float>& y, hoNDArray<float>& r);
    template EXPORTGTPLUS bool multiplyOver5thDimension(const hoNDArray<double>& x, const hoNDArray<double>& y, hoNDArray<double>& r);
    template EXPORTGTPLUS bool multiplyOver5thDimension(const hoNDArray<GT_Complex8>& x, const hoNDArray<GT_Complex8>& y, hoNDArray<GT_Complex8>& r);
    template EXPORTGTPLUS bool multiplyOver5thDimension(const hoNDArray<GT_Complex16>& x, const hoNDArray<GT_Complex16>& y, hoNDArray<GT_Complex16>& r);

    template EXPORTGTPLUS bool multiplyOver5thDimensionExcept(const hoNDArray<float>& x, const hoNDArray<float>& y, size_t n, hoNDArray<float>& r, bool copyY2R);
    template EXPORTGTPLUS bool multiplyOver5thDimensionExcept(const hoNDArray<double>& x, const hoNDArray<double>& y, size_t n, hoNDArray<double>& r, bool copyY2R);
    template EXPORTGTPLUS bool multiplyOver5thDimensionExcept(const hoNDArray<GT_Complex8>& x, const hoNDArray<GT_Complex8>& y, size_t n, hoNDArray<GT_Complex8>& r, bool copyY2R);
    template EXPORTGTPLUS bool multiplyOver5thDimensionExcept(const hoNDArray<GT_Complex16>& x, const hoNDArray<GT_Complex16>& y, size_t n, hoNDArray<GT_Complex16>& r, bool copyY2R);

    template EXPORTGTPLUS bool multipleAdd(const hoNDArray<float>& x, const hoNDArray<float>& y, hoNDArray<float>& r);
    template EXPORTGTPLUS bool multipleAdd(const hoNDArray<double>& x, const hoNDArray<double>& y, hoNDArray<double>& r);
    template EXPORTGTPLUS bool multipleAdd(const hoNDArray<GT_Complex8>& x, const hoNDArray<GT_Complex8>& y, hoNDArray<GT_Complex8>& r);
    template EXPORTGTPLUS bool multipleAdd(const hoNDArray<GT_Complex16>& x, const hoNDArray<GT_Complex16>& y, hoNDArray<GT_Complex16>& r);

    template EXPORTGTPLUS bool multipleMultiply(const hoNDArray<float>& x, const hoNDArray<float>& y, hoNDArray<float>& r);
    template EXPORTGTPLUS bool multipleMultiply(const hoNDArray<double>& x, const hoNDArray<double>& y, hoNDArray<double>& r);
    template EXPORTGTPLUS bool multipleMultiply(const hoNDArray<GT_Complex8>& x, const hoNDArray<GT_Complex8>& y, hoNDArray<GT_Complex8>& r);
    template EXPORTGTPLUS bool multipleMultiply(const hoNDArray<GT_Complex16>& x, const hoNDArray<GT_Complex16>& y, hoNDArray<GT_Complex16>& r);

    template EXPORTGTPLUS bool cropUpTo11DArray(const hoNDArray<short>& x, hoNDArray<short>& r, const std::vector<size_t>& start, std::vector<size_t>& size);
    template EXPORTGTPLUS bool cropUpTo11DArray(const hoNDArray<unsigned short>& x, hoNDArray<unsigned short>& r, const std::vector<size_t>& start, std::vector<size_t>& size);
    template EXPORTGTPLUS bool cropUpTo11DArray(const hoNDArray<float>& x, hoNDArray<float>& r, const std::vector<size_t>& start, std::vector<size_t>& size);
    template EXPORTGTPLUS bool cropUpTo11DArray(const hoNDArray<double>& x, hoNDArray<double>& r, const std::vector<size_t>& start, std::vector<size_t>& size);
    template EXPORTGTPLUS bool cropUpTo11DArray(const hoNDArray<GT_Complex8>& x, hoNDArray<GT_Complex8>& r, const std::vector<size_t>& start, std::vector<size_t>& size);
    template EXPORTGTPLUS bool cropUpTo11DArray(const hoNDArray<GT_Complex16>& x, hoNDArray<GT_Complex16>& r, const std::vector<size_t>& start, std::vector<size_t>& size);

    template EXPORTGTPLUS bool setSubArrayUpTo11DArray(const hoNDArray<short>& x, hoNDArray<short>& r, const std::vector<size_t>& start, std::vector<size_t>& size);
    template EXPORTGTPLUS bool setSubArrayUpTo11DArray(const hoNDArray<unsigned short>& x, hoNDArray<unsigned short>& r, const std::vector<size_t>& start, std::vector<size_t>& size);
    template EXPORTGTPLUS bool setSubArrayUpTo11DArray(const hoNDArray<float>& x, hoNDArray<float>& r, const std::vector<size_t>& start, std::vector<size_t>& size);
    template EXPORTGTPLUS bool setSubArrayUpTo11DArray(const hoNDArray<double>& x, hoNDArray<double>& r, const std::vector<size_t>& start, std::vector<size_t>& size);
    template EXPORTGTPLUS bool setSubArrayUpTo11DArray(const hoNDArray<GT_Complex8>& x, hoNDArray<GT_Complex8>& r, const std::vector<size_t>& start, std::vector<size_t>& size);
    template EXPORTGTPLUS bool setSubArrayUpTo11DArray(const hoNDArray<GT_Complex16>& x, hoNDArray<GT_Complex16>& r, const std::vector<size_t>& start, std::vector<size_t>& size);

    template EXPORTGTPLUS bool extractSampledLinesUpTo11DArray(const hoNDArray<short>& x, hoNDArray<short>& r, const hoNDArray<float>& timeStamp, double acceFactorE1, double acceFactorE2);
    template EXPORTGTPLUS bool extractSampledLinesUpTo11DArray(const hoNDArray<unsigned short>& x, hoNDArray<unsigned short>& r, const hoNDArray<float>& timeStamp, double acceFactorE1, double acceFactorE2);
    template EXPORTGTPLUS bool extractSampledLinesUpTo11DArray(const hoNDArray<float>& x, hoNDArray<float>& r, const hoNDArray<float>& timeStamp, double acceFactorE1, double acceFactorE2);
    template EXPORTGTPLUS bool extractSampledLinesUpTo11DArray(const hoNDArray<double>& x, hoNDArray<double>& r, const hoNDArray<float>& timeStamp, double acceFactorE1, double acceFactorE2);
    template EXPORTGTPLUS bool extractSampledLinesUpTo11DArray(const hoNDArray<GT_Complex8>& x, hoNDArray<GT_Complex8>& r, const hoNDArray<float>& timeStamp, double acceFactorE1, double acceFactorE2);
    template EXPORTGTPLUS bool extractSampledLinesUpTo11DArray(const hoNDArray<GT_Complex16>& x, hoNDArray<GT_Complex16>& r, const hoNDArray<float>& timeStamp, double acceFactorE1, double acceFactorE2);

    template EXPORTGTPLUS bool fillSampledLinesUpTo11DArray(const hoNDArray<short>& x, hoNDArray<short>& r, const hoNDArray<float>& timeStamp);
    template EXPORTGTPLUS bool fillSampledLinesUpTo11DArray(const hoNDArray<unsigned short>& x, hoNDArray<unsigned short>& r, const hoNDArray<float>& timeStamp);
    template EXPORTGTPLUS bool fillSampledLinesUpTo11DArray(const hoNDArray<float>& x, hoNDArray<float>& r, const hoNDArray<float>& timeStamp);
    template EXPORTGTPLUS bool fillSampledLinesUpTo11DArray(const hoNDArray<double>& x, hoNDArray<double>& r, const hoNDArray<float>& timeStamp);
    template EXPORTGTPLUS bool fillSampledLinesUpTo11DArray(const hoNDArray<GT_Complex8>& x, hoNDArray<GT_Complex8>& r, const hoNDArray<float>& timeStamp);
    template EXPORTGTPLUS bool fillSampledLinesUpTo11DArray(const hoNDArray<GT_Complex16>& x, hoNDArray<GT_Complex16>& r, const hoNDArray<float>& timeStamp);

    template EXPORTGTPLUS bool cropOver3rdDimension(const hoNDArray<short>& x, hoNDArray<short>& r, size_t start, size_t end);
    template EXPORTGTPLUS bool cropOver3rdDimension(const hoNDArray<unsigned short>& x, hoNDArray<unsigned short>& r, size_t start, size_t end);
    template EXPORTGTPLUS bool cropOver3rdDimension(const hoNDArray<float>& x, hoNDArray<float>& r, size_t start, size_t end);
    template EXPORTGTPLUS bool cropOver3rdDimension(const hoNDArray<double>& x, hoNDArray<double>& r, size_t start, size_t end);
    template EXPORTGTPLUS bool cropOver3rdDimension(const hoNDArray<GT_Complex8>& x, hoNDArray<GT_Complex8>& r, size_t start, size_t end);
    template EXPORTGTPLUS bool cropOver3rdDimension(const hoNDArray<GT_Complex16>& x, hoNDArray<GT_Complex16>& r, size_t start, size_t end);

    template EXPORTGTPLUS bool setSubArrayOver3rdDimension(const hoNDArray<short>& x, hoNDArray<short>& r, size_t start, size_t end);
    template EXPORTGTPLUS bool setSubArrayOver3rdDimension(const hoNDArray<unsigned short>& x, hoNDArray<unsigned short>& r, size_t start, size_t end);
    template EXPORTGTPLUS bool setSubArrayOver3rdDimension(const hoNDArray<float>& x, hoNDArray<float>& r, size_t start, size_t end);
    template EXPORTGTPLUS bool setSubArrayOver3rdDimension(const hoNDArray<double>& x, hoNDArray<double>& r, size_t start, size_t end);
    template EXPORTGTPLUS bool setSubArrayOver3rdDimension(const hoNDArray<GT_Complex8>& x, hoNDArray<GT_Complex8>& r, size_t start, size_t end);
    template EXPORTGTPLUS bool setSubArrayOver3rdDimension(const hoNDArray<GT_Complex16>& x, hoNDArray<GT_Complex16>& r, size_t start, size_t end);

    template EXPORTGTPLUS bool stdOver3rdDimension(const hoNDArray<float>& x, hoNDArray<float>& std, bool NMinusOne);
    template EXPORTGTPLUS bool stdOver3rdDimension(const hoNDArray<double>& x, hoNDArray<double>& std, bool NMinusOne);
    template EXPORTGTPLUS bool stdOver3rdDimension(const hoNDArray<GT_Complex8>& x, hoNDArray<GT_Complex8>& std, bool NMinusOne);
    template EXPORTGTPLUS bool stdOver3rdDimension(const hoNDArray<GT_Complex16>& x, hoNDArray<GT_Complex16>& std, bool NMinusOne);

    template EXPORTGTPLUS bool permuteE2To3rdDimension(const hoNDArray<float>& x, hoNDArray<float>& r);
    template EXPORTGTPLUS bool permuteE2To3rdDimension(const hoNDArray<double>& x, hoNDArray<double>& r);
    template EXPORTGTPLUS bool permuteE2To3rdDimension(const hoNDArray<GT_Complex8>& x, hoNDArray<GT_Complex8>& r);
    template EXPORTGTPLUS bool permuteE2To3rdDimension(const hoNDArray<GT_Complex16>& x, hoNDArray<GT_Complex16>& r);

    template EXPORTGTPLUS bool permuteE2To5thDimension(const hoNDArray<float>& x, hoNDArray<float>& r);
    template EXPORTGTPLUS bool permuteE2To5thDimension(const hoNDArray<double>& x, hoNDArray<double>& r);
    template EXPORTGTPLUS bool permuteE2To5thDimension(const hoNDArray<GT_Complex8>& x, hoNDArray<GT_Complex8>& r);
    template EXPORTGTPLUS bool permuteE2To5thDimension(const hoNDArray<GT_Complex16>& x, hoNDArray<GT_Complex16>& r);

    template EXPORTGTPLUS bool permuteROTo3rdDimensionFor3DRecon(const hoNDArray<float>& x, hoNDArray<float>& r);
    template EXPORTGTPLUS bool permuteROTo3rdDimensionFor3DRecon(const hoNDArray<double>& x, hoNDArray<double>& r);
    template EXPORTGTPLUS bool permuteROTo3rdDimensionFor3DRecon(const hoNDArray<GT_Complex8>& x, hoNDArray<GT_Complex8>& r);
    template EXPORTGTPLUS bool permuteROTo3rdDimensionFor3DRecon(const hoNDArray<GT_Complex16>& x, hoNDArray<GT_Complex16>& r);

    template EXPORTGTPLUS bool permuteROTo4thDimensionFor3DRecon(const hoNDArray<float>& x, hoNDArray<float>& r);
    template EXPORTGTPLUS bool permuteROTo4thDimensionFor3DRecon(const hoNDArray<double>& x, hoNDArray<double>& r);
    template EXPORTGTPLUS bool permuteROTo4thDimensionFor3DRecon(const hoNDArray<GT_Complex8>& x, hoNDArray<GT_Complex8>& r);
    template EXPORTGTPLUS bool permuteROTo4thDimensionFor3DRecon(const hoNDArray<GT_Complex16>& x, hoNDArray<GT_Complex16>& r);

    template EXPORTGTPLUS bool permuteROTo1stDimensionFor3DRecon(const hoNDArray<float>& x, hoNDArray<float>& r);
    template EXPORTGTPLUS bool permuteROTo1stDimensionFor3DRecon(const hoNDArray<double>& x, hoNDArray<double>& r);
    template EXPORTGTPLUS bool permuteROTo1stDimensionFor3DRecon(const hoNDArray<GT_Complex8>& x, hoNDArray<GT_Complex8>& r);
    template EXPORTGTPLUS bool permuteROTo1stDimensionFor3DRecon(const hoNDArray<GT_Complex16>& x, hoNDArray<GT_Complex16>& r);

    template EXPORTGTPLUS bool permute3rdDimensionTo1stDimension(const hoNDArray<float>& x, hoNDArray<float>& r);
    template EXPORTGTPLUS bool permute3rdDimensionTo1stDimension(const hoNDArray<double>& x, hoNDArray<double>& r);
    template EXPORTGTPLUS bool permute3rdDimensionTo1stDimension(const hoNDArray<GT_Complex8>& x, hoNDArray<GT_Complex8>& r);
    template EXPORTGTPLUS bool permute3rdDimensionTo1stDimension(const hoNDArray<GT_Complex16>& x, hoNDArray<GT_Complex16>& r);

    template EXPORTGTPLUS bool permuteROTo5thDimensionFor3DRecon(const hoNDArray<float>& x, hoNDArray<float>& r);
    template EXPORTGTPLUS bool permuteROTo5thDimensionFor3DRecon(const hoNDArray<double>& x, hoNDArray<double>& r);
    template EXPORTGTPLUS bool permuteROTo5thDimensionFor3DRecon(const hoNDArray<GT_Complex8>& x, hoNDArray<GT_Complex8>& r);
    template EXPORTGTPLUS bool permuteROTo5thDimensionFor3DRecon(const hoNDArray<GT_Complex16>& x, hoNDArray<GT_Complex16>& r);

    template EXPORTGTPLUS bool imageDomainUnwrapping2D(const hoNDArray<GT_Complex8>& x, const hoNDArray<GT_Complex8>& ker, hoNDArray<GT_Complex8>& buf, hoNDArray<GT_Complex8>& y);
    template EXPORTGTPLUS bool imageDomainUnwrapping2D(const hoNDArray<GT_Complex16>& x, const hoNDArray<GT_Complex16>& ker, hoNDArray<GT_Complex16>& buf, hoNDArray<GT_Complex16>& y);

    template EXPORTGTPLUS bool imageDomainUnwrapping2DT(const hoNDArray<GT_Complex8>& x, const hoNDArray<GT_Complex8>& ker, hoNDArray<GT_Complex8>& buf, hoNDArray<GT_Complex8>& y);
    template EXPORTGTPLUS bool imageDomainUnwrapping2DT(const hoNDArray<GT_Complex16>& x, const hoNDArray<GT_Complex16>& ker, hoNDArray<GT_Complex16>& buf, hoNDArray<GT_Complex16>& y);

    template EXPORTGTPLUS bool computePeriodicBoundaryValues(const hoNDArray<float>& x, const hoNDArray<float>& y, float start, float end, hoNDArray<float>& vx, hoNDArray<float>& vy);
    template EXPORTGTPLUS bool computePeriodicBoundaryValues(const hoNDArray<float>& x, const hoNDArray<double>& y, float start, float end, hoNDArray<float>& vx, hoNDArray<double>& vy);
    template EXPORTGTPLUS bool computePeriodicBoundaryValues(const hoNDArray<float>& x, const hoNDArray<GT_Complex8>& y, float start, float end, hoNDArray<float>& vx, hoNDArray<GT_Complex8>& vy);
    template EXPORTGTPLUS bool computePeriodicBoundaryValues(const hoNDArray<float>& x, const hoNDArray<GT_Complex16>& y, float start, float end, hoNDArray<float>& vx, hoNDArray<GT_Complex16>& vy);

    template EXPORTGTPLUS bool computePeriodicBoundaryValues(const hoNDArray<double>& x, const hoNDArray<double>& y, double start, double end, hoNDArray<double>& vx, hoNDArray<double>& vy);
    template EXPORTGTPLUS bool computePeriodicBoundaryValues(const hoNDArray<double>& x, const hoNDArray<float>& y, double start, double end, hoNDArray<double>& vx, hoNDArray<float>& vy);
    template EXPORTGTPLUS bool computePeriodicBoundaryValues(const hoNDArray<double>& x, const hoNDArray<GT_Complex8>& y, double start, double end, hoNDArray<double>& vx, hoNDArray<GT_Complex8>& vy);
    template EXPORTGTPLUS bool computePeriodicBoundaryValues(const hoNDArray<double>& x, const hoNDArray<GT_Complex16>& y, double start, double end, hoNDArray<double>& vx, hoNDArray<GT_Complex16>& vy);
}
