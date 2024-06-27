#include "b1_map.h"
#include "cuNDArray_operators.h"
#include "cuNDArray_elemwise.h"
#include "vector_td_utilities.h"
#include "real_utilities.h"
#include "real_utilities_device.h"
#include "complext.h"
#include "check_CUDA.h"
#include "cudaDeviceManager.h"
#include "setup_grid.h"
#include "hoNDArray_fileio.h"
#include "GPUTimer.h"

#include "CUBLASContextProvider.h"
#include <cublas_v2.h>

#include <iostream>
#include <cmath>

using namespace std;

namespace Gadgetron{

    template <class T> int write_cuNDArray_to_disk(cuNDArray<T>* a, const char* filename)
    {
        boost::shared_ptr< hoNDArray<T> > host = a->to_host();
        write_nd_array<complext<float> >(host.get(), filename);
        return 0;
    }

    extern __shared__ char _shared_mem[];

    //
    // Main method
    //

    template<class REAL> EXPORTGPUPMRI bool
    estimate_b1_map_2D_NIH_Souheil( cuNDArray<complext<REAL> >* data, cuNDArray<complext<REAL> >* csm, size_t ks, size_t power, 
                                    cuNDArray<complext<REAL> >& D, cuNDArray<complext<REAL> >& DH_D, 
                                    cuNDArray<complext<REAL> >& V1, cuNDArray<complext<REAL> >& U1)
    {
        if( data->get_number_of_dimensions() < 2 )
        {
            GINFO_STREAM(endl << "estimate_b1_map_2D_NIH_Souheil:: dimensionality mismatch." << endl);
            return false;
        }

        if ( !csm->dimensions_equal(data) )
        {
            csm->create(*data->get_dimensions());
        }

        size_t kss = ks*ks;

        {
            assemble_D( data, &D, ks );
        }

        //{
        //    std::string dstDir = "D:/software/Gadgetron/20130114/gadgetron/toolboxes/gtplus/ut/result/";
        //    std::string filename = dstDir + "D.cplx";
        //    write_cuNDArray_to_disk(&D, filename.c_str());
        //}

        {
            computeDH_D( data, &D, &DH_D, kss );
        }

        //{
        //    std::string dstDir = "D:/software/Gadgetron/20130114/gadgetron/toolboxes/gtplus/ut/result/";
        //    std::string filename = dstDir + "DH_D.cplx";
        //    write_cuNDArray_to_disk(&DH_D, filename.c_str());
        //}

        {
            computeV1( data, &D, &DH_D, &V1, csm, power, kss);
        }

        //{
        //    std::string dstDir = "D:/software/Gadgetron/20130114/gadgetron/toolboxes/gtplus/ut/result/";
        //    std::string filename = dstDir + "V1.cplx";
        //    write_cuNDArray_to_disk(&V1, filename.c_str());
        //}

        {
            computeU1( data, &D, &V1, &U1, kss);
        }

        //{
        //    std::string dstDir = "D:/software/Gadgetron/20130114/gadgetron/toolboxes/gtplus/ut/result/";
        //    std::string filename = dstDir + "U1.cplx";
        //    write_cuNDArray_to_disk(&U1, filename.c_str());
        //}

        {
            extract_csm( data, &V1, &U1, csm, kss);
        }

        //{
        //    std::string dstDir = "D:/software/Gadgetron/20130114/gadgetron/toolboxes/gtplus/ut/result/";
        //    std::string filename = dstDir + "csm.cplx";
        //    write_cuNDArray_to_disk(csm, filename.c_str());
        //}

        return true;
    }

    // assemble_D
    template<class T> __global__ void
    assemble_D_kernel( const T* __restrict__ pData, T* __restrict__ pD, int RO, int E1, int N, int CHA, int kss, int halfKs )
    {
        typedef typename realType<T>::Type REAL;

        const unsigned int cha = threadIdx.y;

        unsigned int n = (blockIdx.x*blockDim.x + threadIdx.x)/(RO*E1);
        unsigned int e1 = (blockIdx.x*blockDim.x + threadIdx.x - n*RO*E1)/RO;
        unsigned int ro = (blockIdx.x*blockDim.x + threadIdx.x - n*RO*E1)%RO;

        // printf("ro=%d, e1=%d, cha=%d, n=%d\n", ro, e1, cha, n);

        if( ro<RO && e1<E1 && n<N )
        {
            // printf("ro=%d, e1=%d\n", ro, e1);

            unsigned int idx2D = cha*RO*E1*kss*N + n*RO*E1 + ro + e1*RO;

            int kro, ke1, de1, dro;

            if ( e1>=halfKs && e1<E1-halfKs && ro>=halfKs && ro<RO-halfKs )
            {
                // printf("e1>=halfKs && e1<E1-halfKs && ro>=halfKs && ro<RO-halfKs\n");

                const T* pDataCurr = pData + n*RO*E1 + cha*RO*E1*N;

                int ind=0;
                for ( ke1=-halfKs; ke1<=halfKs; ke1++ )
                {
                    de1 = e1 + ke1;
                    for ( kro=-halfKs; kro<=halfKs; kro++ )
                    {
                        pD[ind*RO*E1*N + idx2D] = pDataCurr[de1*RO+ro+kro];
                        //printf("pD[idxD]=%f\n", pD[idxD].real());
                        ind++;
                    }
                }
            }
            else
            {
                // printf("boundary\n");
                const T* pDataCurr = pData + n*RO*E1 + cha*RO*E1*N;
                int ind=0;
                for ( ke1=-halfKs; ke1<=halfKs; ke1++ )
                {
                    de1 = e1 + ke1;
                    if ( de1 < 0 ) de1 += E1;
                    if ( de1 >= E1 ) de1 -= E1;

                    for ( kro=-halfKs; kro<=halfKs; kro++ )
                    {
                        dro = ro + kro;
                        if ( dro < 0 ) dro += RO;
                        if ( dro >= RO ) dro -= RO;

                        pD[ind*RO*E1*N+ idx2D] = pDataCurr[de1*RO+dro];
                        ind++;
                    }
                }
            }
        }
    }

    template<class T>
    void assemble_D( cuNDArray<T>* data, cuNDArray<T>* D, size_t ks )
    {
        size_t RO = data->get_size(0);
        size_t E1 = data->get_size(1);
        size_t N(1), CHA;

        size_t NDim = data->get_number_of_dimensions();

        if ( NDim == 3 )
        {
            CHA = data->get_size(2);
        }

        if ( NDim == 4 )
        {
            N = data->get_size(2);
            CHA = data->get_size(3);
        }

        if ( ks%2 != 1 )
        {
            ks++;
        }

        size_t halfKs = ks/2;

        // Setup block/grid dimensions
        int cur_device = cudaDeviceManager::Instance()->getCurrentDevice();
        int warp_size = cudaDeviceManager::Instance()->warp_size(cur_device);
        int max_blockdim = cudaDeviceManager::Instance()->max_blockdim(cur_device);
        dim3 blockDim(((max_blockdim/CHA)/warp_size)*warp_size, CHA);

        if( blockDim.x == 0 )
        {
            blockDim.x = warp_size;
            while ( blockDim.x*CHA*CHA > max_blockdim && blockDim.x>1 )
            {
                blockDim.x /= 2;
            }

            if ( blockDim.x <= 1 )
            {
                blockDim.x = 1;
            }
        }

        dim3 gridDim((RO*E1*N+blockDim.x-1)/blockDim.x);

        // Invoke kernel
        assemble_D_kernel<T><<< gridDim, blockDim >>>( data->get_data_ptr(), D->get_data_ptr(), RO, E1, N, CHA, ks*ks, halfKs );

        CHECK_FOR_CUDA_ERROR();
    }

    // compute DH_D
    template<class T> __global__ void
    computeDH_D_kernel( const T* __restrict__ pD, T* __restrict__ pDH_D, int RO, int E1, int N, int CHA, int kss )
    {
        typedef typename realType<T>::Type REAL;

        // DH_D, [RO E1 CHA CHA_Prime]
        const unsigned int cha = threadIdx.y;
        const unsigned int cha_prime = threadIdx.z;

        unsigned int n = (blockIdx.x*blockDim.x + threadIdx.x)/(RO*E1);
        unsigned int e1 = (blockIdx.x*blockDim.x + threadIdx.x - n*RO*E1)/RO;
        unsigned int ro = (blockIdx.x*blockDim.x + threadIdx.x - n*RO*E1)%RO;

        if( ro<RO && e1<E1 && n<N )
        {
            unsigned int idx = ro + e1*RO + n*RO*E1;

            // every thread compute an element of DH_D for a pixel
            int k;
            T v;
            v = 0;
            for ( k=0; k<kss; k++ )
            {
                v += conj(pD[cha*RO*E1*N*kss + k*RO*E1*N + idx])*pD[cha_prime*RO*E1*N*kss + k*RO*E1*N + idx];
            }

            pDH_D[cha_prime*RO*E1*N*CHA + cha*RO*E1*N + idx] = v;
        }
    }

    // use the shared memory
    template<class T> __global__ void
    computeDH_D_kernel3( const T*  __restrict__ pD, T* __restrict__ pDH_D, int RO, int E1, int N, int CHA, int kss, int ks, int num )
    {
        typedef typename realType<T>::Type REAL;

        // DH_D, [RO E1 CHA CHA_Prime]
        const unsigned int cha = threadIdx.y;
        const unsigned int cha_prime = threadIdx.z;

        unsigned int n = (blockIdx.x*blockDim.x + threadIdx.x)/(RO*E1);
        unsigned int e1 = (blockIdx.x*blockDim.x + threadIdx.x - n*RO*E1)/RO;
        unsigned int ro = (blockIdx.x*blockDim.x + threadIdx.x - n*RO*E1)%RO;

        if( ro<RO && e1<E1 && n<N )
        {
            unsigned int idx = ro + e1*RO + n*RO*E1;
            unsigned int idxD = idx + cha*RO*E1*N*kss;
            unsigned int idxShared = threadIdx.x*kss*CHA;

            T *shared_mem = (T*) _shared_mem;

            int k;

            if ( cha_prime == 0 )
            {
                for ( k=0; k<kss; k++ )
                {
                    shared_mem[idxShared + k + cha*kss ] = pD[idxD + k*RO*E1*N ];
                }
            }

            __syncthreads();

            T v = conj(shared_mem[idxShared + cha*kss])*shared_mem[idxShared + cha_prime*kss];
            for ( k=1; k<kss; k++ )
            {
                v += conj(shared_mem[idxShared + cha*kss + k])*shared_mem[idxShared + cha_prime*kss + k];
            }

            pDH_D[cha_prime*RO*E1*N*CHA + cha*RO*E1*N + idx] = v;
        }
    }

    template<class T>
    void computeDH_D( cuNDArray<T>* data, cuNDArray<T>* D, cuNDArray<T>* DH_D, size_t kss )
    {
        size_t RO = data->get_size(0);
        size_t E1 = data->get_size(1);
        size_t N(1), CHA;

        size_t NDim = data->get_number_of_dimensions();

        if ( NDim == 3 )
        {
            CHA = data->get_size(2);
        }

        if ( NDim == 4 )
        {
            N = data->get_size(2);
            CHA = data->get_size(3);
        }

        // Setup block/grid dimensions
        int cur_device = cudaDeviceManager::Instance()->getCurrentDevice();
        int warp_size = cudaDeviceManager::Instance()->warp_size(cur_device);
        int max_blockdim = cudaDeviceManager::Instance()->max_blockdim(cur_device);
        size_t shared_mem_per_block = cudaDeviceManager::Instance()->shared_mem_per_block(cur_device);

        // estimate how many pixels a block can process
        size_t ks = (size_t)std::sqrt((double)kss);

        // size_t numOfPixels = shared_mem_per_block/4/(sizeof(T)*(kss+ks)*CHA);
        size_t numOfPixels = shared_mem_per_block/4/(sizeof(T)*kss*CHA);

        while ( numOfPixels*ks*CHA>max_blockdim && numOfPixels>0 )
        {
            numOfPixels--;
        }

        if ( numOfPixels > 0 )
        {
            dim3 blockDim(numOfPixels, CHA, CHA);

            dim3 gridDim((RO*E1*N+blockDim.x-1)/blockDim.x);

            computeDH_D_kernel3<T><<< gridDim, blockDim, numOfPixels*sizeof(T)*kss*CHA >>>( D->get_data_ptr(), DH_D->get_data_ptr(), RO, E1, N, CHA, kss, ks, numOfPixels );
        }
        else
        {
            dim3 blockDim(((max_blockdim/(CHA*CHA))/warp_size)*warp_size, CHA, CHA);

            if( blockDim.x == 0 )
            {
                blockDim.x = warp_size;
                while ( blockDim.x*CHA*CHA > max_blockdim && blockDim.x>1 )
                {
                    blockDim.x /= 2;
                }

                if ( blockDim.x <= 1 )
                {
                    blockDim.x = 1;
                }
            }

            dim3 gridDim((RO*E1*N+blockDim.x-1)/blockDim.x);

            // Invoke kernel
            computeDH_D_kernel<T><<< gridDim, blockDim >>>( D->get_data_ptr(), DH_D->get_data_ptr(), RO, E1, N, CHA, kss );
        }

        CHECK_FOR_CUDA_ERROR();
    }

    // compute V1
    template<class T> __global__ void
    computeV1_kernel( const T* __restrict__ pD, T* __restrict__ pV1, int RO, int E1, int N, int CHA, int kss )
    {
        typedef typename realType<T>::Type REAL;

        const unsigned int cha = threadIdx.y;
        unsigned int n = (blockIdx.x*blockDim.x + threadIdx.x)/(RO*E1);
        unsigned int e1 = (blockIdx.x*blockDim.x + threadIdx.x - n*RO*E1)/RO;
        unsigned int ro = (blockIdx.x*blockDim.x + threadIdx.x - n*RO*E1)%RO;

        if( ro<RO && e1<E1 && n<N )
        {
            unsigned int idx = ro + e1*RO + n*RO*E1;
            unsigned int idxD = cha*RO*E1*N*kss + idx;

            T v = 0;
            for ( int ii=0; ii<kss; ii++ )
            {
                v += pD[idxD + ii*RO*E1*N];
            }
            pV1[cha*RO*E1*N + idx] = v;
        }
    }

    template<class T> __global__ void
    power_method_kernel( const T* __restrict__ pDH_D, T* __restrict__ pV1,  T* __restrict__ pV, unsigned int RO, unsigned int E1, unsigned int N, unsigned int CHA, unsigned int kss, unsigned int power )
    {
        typedef typename realType<T>::Type REAL;

        const unsigned int ro = blockIdx.x*blockDim.x+threadIdx.x;
        const unsigned int e1 = blockIdx.y*blockDim.y+threadIdx.y;
        unsigned int n = blockIdx.z;

        if( ro<RO && e1<E1 && n<N )
        {
            unsigned int cha;

            unsigned int idx2D = ro + e1*RO + n*RO*E1;

            unsigned int N3D = RO*E1*N;

            REAL v1Norm(0);
            for ( cha=0; cha<CHA; cha++ )
            {
                v1Norm += norm(pV1[cha*N3D + idx2D]);
            }
            v1Norm = ::sqrt(v1Norm);

            for ( cha=0; cha<CHA; cha++ )
            {
                pV1[cha*N3D + idx2D] /= v1Norm;
            }

            unsigned int po;
            for ( po=0; po<power; po++ )
            {
                for( unsigned j=0; j<CHA; j++)
                {
                    T v = 0;
                    for( unsigned int k=0; k<CHA; k++)
                    {
                        v += pDH_D[k*CHA*N3D+j*N3D+idx2D]*pV1[k*N3D+idx2D];
                    }
                    pV[j*N3D+idx2D] = v;
                }

                for ( cha=0; cha<CHA; cha++ )
                {
                    pV1[cha*N3D + idx2D] = pV[cha*N3D + idx2D];
                }

                v1Norm = 0;
                for ( cha=0; cha<CHA; cha++ )
                {
                    v1Norm += norm(pV1[cha*N3D + idx2D]);
                }
                v1Norm = 1/std::sqrt(v1Norm);

                for ( cha=0; cha<CHA; cha++ )
                {
                    pV1[cha*N3D + idx2D] *= v1Norm;
                }
            }
        }
    }

    template<class T>
    void computeV1( cuNDArray<T>* data, cuNDArray<T>* D, cuNDArray<T>* DH_D, cuNDArray<T>* V1, cuNDArray<T>* V, int power, int kss)
    {
        size_t RO = data->get_size(0);
        size_t E1 = data->get_size(1);
        size_t N(1), CHA;

        size_t NDim = data->get_number_of_dimensions();

        if ( NDim == 3 )
        {
            CHA = data->get_size(2);
        }

        if ( NDim == 4 )
        {
            N = data->get_size(2);
            CHA = data->get_size(3);
        }

        // Setup block/grid dimensions
        int cur_device = cudaDeviceManager::Instance()->getCurrentDevice();
        int warp_size = cudaDeviceManager::Instance()->warp_size(cur_device);
        int max_blockdim = cudaDeviceManager::Instance()->max_blockdim(cur_device);
        dim3 blockDim(((max_blockdim/CHA)/warp_size)*warp_size, CHA);

        if( blockDim.x == 0 )
        {
            GERROR_STREAM("blockDim.x == 0");
            throw std::runtime_error("computeDH_D: dimension exceeds device capacity.");
        }

        dim3 gridDim((RO*E1*N+blockDim.x-1)/blockDim.x);

        // Invoke kernel
        computeV1_kernel<T><<< gridDim, blockDim >>>( D->get_data_ptr(), V1->get_data_ptr(), RO, E1, N, CHA, kss );

        // power method
        dim3 blockDim2(16, 16);
        dim3 gridDim2((RO+blockDim2.x-1)/blockDim2.x, (E1+blockDim2.y-1)/blockDim2.y, N);

        power_method_kernel<T><<< gridDim2, blockDim2 >>>( DH_D->get_data_ptr(), V1->get_data_ptr(), V->get_data_ptr(), RO, E1, N, CHA, kss, power );

        CHECK_FOR_CUDA_ERROR();
    }

    // compute U1
    template<class T> __global__ void
    computeU1_kernel(const  T* __restrict__ pD, const T* __restrict__ pV1, T* __restrict__ pU1, int RO, int E1, int N, int CHA, int kss )
    {
        typedef typename realType<T>::Type REAL;

        const unsigned int k = threadIdx.y;
        unsigned int n = (blockIdx.x*blockDim.x + threadIdx.x)/(RO*E1);
        unsigned int e1 = (blockIdx.x*blockDim.x + threadIdx.x - n*RO*E1)/RO;
        unsigned int ro = (blockIdx.x*blockDim.x + threadIdx.x - n*RO*E1)%RO;

        if( ro<RO && e1<E1 && n<N )
        {
            unsigned int idx = ro + e1*RO + n*RO*E1;
            unsigned int idxD = k*RO*E1*N + idx;

            T v = 0;
            for ( int ii=0; ii<CHA; ii++ )
            {
                v += pD[idxD + ii*kss*RO*E1*N] * pV1[ii*RO*E1*N+idx];
            }
            pU1[k*RO*E1*N + idx] = v;
        }
    }

    template<class T>
    void computeU1( cuNDArray<T>* data, cuNDArray<T>* D, cuNDArray<T>* V1, cuNDArray<T>* U1, int kss)
    {
        size_t RO = data->get_size(0);
        size_t E1 = data->get_size(1);
        size_t N(1), CHA;

        size_t NDim = data->get_number_of_dimensions();

        if ( NDim == 3 )
        {
            CHA = data->get_size(2);
        }

        if ( NDim == 4 )
        {
            N = data->get_size(2);
            CHA = data->get_size(3);
        }

        // Setup block/grid dimensions
        int cur_device = cudaDeviceManager::Instance()->getCurrentDevice();
        int warp_size = cudaDeviceManager::Instance()->warp_size(cur_device);
        int max_blockdim = cudaDeviceManager::Instance()->max_blockdim(cur_device);
        dim3 blockDim(((max_blockdim/kss)/warp_size)*warp_size, kss);

        if( blockDim.x == 0 )
        {
            // GERROR_STREAM("blockDim.x == 0");
            blockDim.x = warp_size;
            while ( blockDim.x*kss > max_blockdim && blockDim.x>1 )
            {
                blockDim.x /= 2;
            }

            if ( blockDim.x <= 1 )
            {
                blockDim.x = 1;
            }
        }

        dim3 gridDim((RO*E1*N+blockDim.x-1)/blockDim.x);

        // Invoke kernel
        computeU1_kernel<T><<< gridDim, blockDim >>>( D->get_data_ptr(), V1->get_data_ptr(), U1->get_data_ptr(), RO, E1, N, CHA, kss );

        CHECK_FOR_CUDA_ERROR();
    }

    // extract the csm
    template<class T> __global__ void
    extract_csm_kernel( const T* __restrict__ pV1, const T* __restrict__ pU1, T* __restrict__ pCSM, unsigned int RO, unsigned int E1, unsigned int N, unsigned int CHA, unsigned int kss )
    {
        typedef typename realType<T>::Type REAL;

        const unsigned int ro = blockIdx.x*blockDim.x+threadIdx.x;
        const unsigned int e1 = blockIdx.y*blockDim.y+threadIdx.y;
        unsigned int n = blockIdx.z;

        if( ro<RO && e1<E1 && n<N )
        {
            unsigned int cha;
            unsigned int idx = ro + e1*RO + n*RO*E1;

            T phaseU1 = pU1[idx];
            for ( int po=1; po<kss; po++ )
            {
                phaseU1 += pU1[idx + po*RO*E1*N];
            }
            phaseU1 /= abs(phaseU1);

            // put the mean object phase to coil map
            for ( cha=0; cha<CHA; cha++ )
            {
                pCSM[cha*RO*E1*N+idx] = phaseU1 * conj(pV1[cha*RO*E1*N+idx]);
            }
        }
    }

    template<class T>
    void extract_csm( cuNDArray<T>* data, cuNDArray<T>* V1, cuNDArray<T>* U1, cuNDArray<T>* csm, int kss)
    {
        size_t RO = data->get_size(0);
        size_t E1 = data->get_size(1);
        size_t N(1), CHA;

        size_t NDim = data->get_number_of_dimensions();

        if ( NDim == 3 )
        {
            CHA = data->get_size(2);
        }

        if ( NDim == 4 )
        {
            N = data->get_size(2);
            CHA = data->get_size(3);
        }

        // Setup block/grid dimensions
        dim3 blockDim(16, 16);
        dim3 gridDim((RO+blockDim.x-1)/blockDim.x, (E1+blockDim.y-1)/blockDim.y, N);

        extract_csm_kernel<T><<< gridDim, blockDim >>>( V1->get_data_ptr(), U1->get_data_ptr(), csm->get_data_ptr(), RO, E1, N, CHA, kss );

        CHECK_FOR_CUDA_ERROR();
    }

    //
    // Template instantiation
    //
    template EXPORTGPUPMRI bool estimate_b1_map_2D_NIH_Souheil<float>( cuNDArray<complext<float> >* data, cuNDArray<complext<float> >* csm, size_t ks, size_t power,
                                    cuNDArray<complext<float> >& D, cuNDArray<complext<float> >& DH_D, 
                                    cuNDArray<complext<float> >& V1, cuNDArray<complext<float> >& U1 );
}
