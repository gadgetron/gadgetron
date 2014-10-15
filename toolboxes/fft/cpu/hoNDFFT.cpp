#include "hoNDFFT.h"
#include "hoMatrix.h"
#include "hoNDArray_math_util.h"
#include "hoNDArray_math.h"

/// uncomment this to disable MKL FFT calls
/// #undef USE_MKL

namespace Gadgetron{

    template<typename T> hoNDFFT<T>* hoNDFFT<T>::instance()
    {
        if (!instance_) instance_ = new hoNDFFT<T>();
        return instance_;
    }

    template<class T> hoNDFFT<T>* hoNDFFT<T>::instance_ = NULL;

    template<class T> void hoNDFFT<T>::fft_int(hoNDArray< ComplexType >* input, size_t dim_to_transform, int sign)
    {
        if (sign != -1 && sign != 1) return;
        if (dim_to_transform >= input->get_number_of_dimensions()) return;

        int stride     = 1;           //Distance between points in transform
        int dist       = 1;           //Distance between vectors
        int trafos     = 1;           //Transformations per chunk
        int chunks     = 1;           //Number of chunks
        int chunk_size = 1;           //Points per chunk
        int length     = 1;           //Length of each transform
        int total_dist = 1;

        T scale = 0.0;

        void* fft_plan        = 0;
        T*    fft_storage     = 0;

        T* fft_buffer = 0;
        T* data_ptr = 0;

        //Set sizes
        length = (int)input->get_size(dim_to_transform);

        if (sign == 1)
        {
            scale = (T)(1.0/length);
        }
        else
        {
            scale = (T)1.0;
        }

        if (dim_to_transform != 0)
        {
            for (size_t i = 0; i < dim_to_transform; i++)
            {
                chunk_size *= (int)input->get_size(i);
            }
            stride = chunk_size;
            trafos = chunk_size;
            chunk_size *= length;

            for (size_t i = dim_to_transform+1; i < input->get_number_of_dimensions(); i++)
            {
                chunks *= (int)input->get_size(i);
            }
        }
        else
        {
            for (size_t i = 1; i < input->get_number_of_dimensions(); i++)
            {
                trafos *= (int)input->get_size(i);
            }
            chunk_size = trafos*length;

            dist = length;
        }

        //*2 real and imag
        chunk_size *= 2;
        dist *= 2;
        total_dist = trafos*dist;


        //Allocate storage and make plan
        {
            mutex_.lock();
            fft_storage = (T*)fftw_malloc_ptr_(sizeof(T)*length*2);
            if (fft_storage == 0)
            {
                std::cout << "Failed to allocate buffer for FFT" << std::endl;
                return;
            }
            fft_buffer = (T*)fft_storage;

            unsigned planner_flags = FFTW_MEASURE | FFTW_DESTROY_INPUT;

            fft_plan = fftw_plan_dft_1d_ptr_(length, fft_storage, fft_storage, sign, planner_flags);

            if (fft_plan == 0)
            {
                fftw_free_ptr_(fft_storage);
                std::cout << "Failed to create plan for FFT" << std::endl;
                return;
            }
            mutex_.unlock();
        }

        //Grab address of data
        data_ptr = reinterpret_cast<T*>(input->get_data_ptr());

        register int idx1_max = chunks*chunk_size;
        register int idx1, idx2;       //Index variables
        register int idx2_limit;
        register int middle_point = ((length+1)>>1)<<1;
        register int length2 = length<<1;
        register int stride2 = stride<<1;

        for (idx1 = 0; idx1 < idx1_max; idx1+=chunk_size) //Loop over all chunks
        {
            idx2_limit = idx1+total_dist;
            for (idx2 = idx1; idx2 < idx2_limit; idx2+=dist) //Loop over all transformations
            {
                ///Copy data to buffer.
                {
                    register int j, idx3 = idx2;
                    for (j = middle_point; j < length2; idx3+=stride2)
                    {
                        fft_buffer[j++] = data_ptr[idx3  ];
                        fft_buffer[j++] = data_ptr[idx3+1];
                    }
                    for (j = 0; j < middle_point; idx3+=stride2)
                    {
                        fft_buffer[j++] = data_ptr[idx3  ];
                        fft_buffer[j++] = data_ptr[idx3+1];
                    }
                }

                fftw_execute_ptr_(fft_plan);

                {
                    register int j, idx3 = idx2;

                    for (j = middle_point; j < length2; idx3+=stride2)
                    {
                        data_ptr[idx3  ] = fft_buffer[j++]*scale;
                        data_ptr[idx3+1] = fft_buffer[j++]*scale;
                    }
                    for (j = 0; j < middle_point; idx3+=stride2)
                    {
                        data_ptr[idx3  ] = fft_buffer[j++]*scale;
                        data_ptr[idx3+1] = fft_buffer[j++]*scale;
                    }
                }

            } //Loop over transformations
        } //Loop over chunks

        //clean up
        {
            mutex_.lock();
            if (fft_plan != 0)
            {
                fftw_destroy_plan_ptr_(fft_plan);
            }

            if (fft_storage != 0)
            {
                fftw_free_ptr_(fft_storage);
            }
            mutex_.unlock();
        }
    }

    template<> void hoNDFFT<float>::set_function_pointers()
    {
        fftw_import_wisdom_from_file_ptr_ = &fftwf_import_wisdom_from_file;
        fftw_export_wisdom_to_file_ptr_ = &fftwf_export_wisdom_to_file;
        fftw_cleanup_ptr_ = &fftwf_cleanup;
        fftw_malloc_ptr_ = &fftwf_malloc;
        fftw_free_ptr_ = &fftwf_free;
        fftw_execute_ptr_ = (void (*)(void*))(&fftwf_execute);
        fftw_plan_dft_1d_ptr_ = (void* (*)(int, void*, void*, int, unsigned))(&fftwf_plan_dft_1d);
        fftw_destroy_plan_ptr_ = (void (*)(void*))(&fftwf_destroy_plan);
    }

    template<> void hoNDFFT<double>::set_function_pointers()
    {
        fftw_import_wisdom_from_file_ptr_ = &fftw_import_wisdom_from_file;
        fftw_export_wisdom_to_file_ptr_ = &fftw_export_wisdom_to_file;
        fftw_cleanup_ptr_ = &fftw_cleanup;
        fftw_malloc_ptr_ = &fftw_malloc;
        fftw_free_ptr_ = &fftw_free;
        fftw_execute_ptr_ = (void (*)(void*))(&fftw_execute);
        fftw_plan_dft_1d_ptr_ = (void* (*)(int, void*, void*, int, unsigned))(&fftw_plan_dft_1d);
        fftw_destroy_plan_ptr_ = (void (*)(void*))(&fftw_destroy_plan);
    }

    template<typename T> 
    inline size_t hoNDFFT<T>::fftshiftPivot(size_t x)
    {
        return (size_t)(ceil(x*0.5));
    }

    template<typename T> 
    inline size_t hoNDFFT<T>::ifftshiftPivot(size_t x)
    {
        return (size_t)(floor(x*0.5));
    }

    template<typename T> 
    inline bool hoNDFFT<T>::fftshift1D(const ComplexType* a, ComplexType* r, size_t x, size_t pivot)
    {
        try
        {
            memcpy(r, a+pivot, sizeof(ComplexType)*(x-pivot));
            memcpy(r+x-pivot, a, sizeof(ComplexType)*pivot);
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Errors in hoNDFFT<T>::fftshift1D(const ComplexType* a, ComplexType* r, size_t x, size_t pivot) ...");
            return false;
        }

        return true;
    }

    template<typename T> 
    inline bool hoNDFFT<T>::ifftshift1D(const ComplexType* a, ComplexType* r, size_t x, size_t pivot)
    {
        return fftshift1D(a, r, x, pivot);
    }

    template<typename T> 
    bool hoNDFFT<T>::fftshiftPivot1D(ComplexType* a, size_t x, size_t n, size_t pivot)
    {
        try
        {
            long long counter;

            #pragma omp parallel private(counter) shared(n, x, pivot, a) if ( n > 256 )
            {
                hoNDArray< ComplexType > aTmp(x);

                #pragma omp for
                for ( counter=0; counter<(long long)n; counter++ )
                {
                    fftshift1D(a+counter*x, aTmp.begin(), x, pivot);
                    memcpy(a+counter*x, aTmp.begin(), sizeof(ComplexType)*x);
                }
            }
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Errors in hoNDFFT<T>::fftshiftPivot1D(ComplexType* a, size_t x, size_t n, size_t pivot) ...");
            return false;
        }

        return true;
    }

    template<typename T> 
    bool hoNDFFT<T>::fftshiftPivot1D(const ComplexType* a, ComplexType* r, size_t x, size_t n, size_t pivot)
    {
        try
        {
            long long counter;

            #pragma omp parallel for private(counter) shared(n, x, pivot, a, r) if ( n > 256 )
            for ( counter=0; counter<(long long)n; counter++ )
            {
                fftshift1D(a+counter*x, r+counter*x, x, pivot);
            }
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Errors in hoNDFFT<T>::fftshiftPivot1D(const ComplexType* a, ComplexType* r, size_t x, size_t n, size_t pivot) ...");
            return false;
        }

        return true;
    }

    template<typename T> 
    bool hoNDFFT<T>::fftshift1D(hoNDArray< ComplexType >& a)
    {
        try
        {
            size_t x = a.get_size(0);
            size_t pivot = fftshiftPivot(x);
            size_t numOfShifts = a.get_number_of_elements()/x;

            GADGET_CHECK_RETURN_FALSE(fftshiftPivot1D(a.begin(), x, numOfShifts, pivot));
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Errors in hoNDFFT<T>::fftshift1D(hoNDArray< ComplexType >& a) ...");
            return false;
        }

        return true;
    }

    template<typename T> 
    bool hoNDFFT<T>::fftshift1D(const hoNDArray< ComplexType >& a, hoNDArray< ComplexType >& r)
    {
        try
        {
            if ( !r.dimensions_equal(&a) )
            {
                r = a;
            }

            size_t x = a.get_size(0);
            size_t pivot = fftshiftPivot(x);
            size_t numOfShifts = a.get_number_of_elements()/x;

            GADGET_CHECK_RETURN_FALSE(fftshiftPivot1D(a.begin(), r.begin(), x, numOfShifts, pivot));
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Errors in hoNDFFT<T>::fftshift1D(const hoNDArray< ComplexType >& a, hoNDArray< ComplexType >& r) ...");
            return false;
        }

        return true;
    }

    template<typename T> 
    bool hoNDFFT<T>::ifftshift1D(hoNDArray< ComplexType >& a)
    {
        try
        {
            size_t x = a.get_size(0);
            size_t pivot = ifftshiftPivot(x);
            size_t numOfShifts = a.get_number_of_elements()/x;

            GADGET_CHECK_RETURN_FALSE(fftshiftPivot1D(a.begin(), x, numOfShifts, pivot));
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Errors in hoNDFFT<T>::ifftshift1D(hoNDArray< ComplexType >& a) ...");
            return false;
        }

        return true;
    }

    template<typename T> 
    bool hoNDFFT<T>::ifftshift1D(const hoNDArray< ComplexType >& a, hoNDArray< ComplexType >& r)
    {
        try
        {
            if ( !r.dimensions_equal(&a) )
            {
                r = a;
            }

            size_t x = a.get_size(0);
            size_t pivot = ifftshiftPivot(x);
            size_t numOfShifts = a.get_number_of_elements()/x;

            GADGET_CHECK_RETURN_FALSE(fftshiftPivot1D(a.begin(), r.begin(), x, numOfShifts, pivot));
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Errors in hoNDFFT<T>::ifftshift1D(const hoNDArray< ComplexType >& a, hoNDArray< ComplexType >& r) ...");
            return false;
        }

        return true;
    }

    template<typename T> 
    bool hoNDFFT<T>::fftshiftPivot2D(const ComplexType* a, ComplexType* r, size_t x, size_t y, size_t n, size_t pivotx, size_t pivoty)
    {
        try
        {
            GADGET_CHECK_RETURN_FALSE(a!=NULL);
            GADGET_CHECK_RETURN_FALSE(r!=NULL);

            long long tt;

            #pragma omp parallel for private(tt) shared(a, r, x, y, n, pivotx, pivoty) if (n>1)
            for ( tt=0; tt<(long long)n; tt++ )
            {
                const ComplexType* ac = a + tt*x*y;
                ComplexType* rc = r + tt*x*y;

                size_t ay, ry;

                for ( ay=pivoty; ay<y; ay++ )
                {
                    ry = ay - pivoty;
                    memcpy(rc+ry*x, ac+ay*x+pivotx, sizeof(ComplexType)*(x-pivotx));
                    memcpy(rc+ry*x+x-pivotx, ac+ay*x, sizeof(ComplexType)*pivotx);
                }

                for ( ay=0; ay<pivoty; ay++ )
                {
                    ry = ay + y - pivoty;
                    memcpy(rc+ry*x, ac+ay*x+pivotx, sizeof(ComplexType)*(x-pivotx));
                    memcpy(rc+ry*x+x-pivotx, ac+ay*x, sizeof(ComplexType)*pivotx);
                }
            }
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Errors in hoNDFFT<T>::fftshiftPivot2D(const ComplexType* a, ComplexType* r, size_t x, size_t y, size_t n, size_t pivotx, size_t pivoty) ...");
            return false;
        }

        return true;
    }

    template<typename T> 
    bool hoNDFFT<T>::fftshiftPivot2D(ComplexType* a, size_t x, size_t y, size_t n, size_t pivotx, size_t pivoty)
    {
        try
        {
            GADGET_CHECK_RETURN_FALSE(a!=NULL);

            long long tt;

            #pragma omp parallel private(tt) shared(a, x, y, n, pivotx, pivoty) if (n>1)
            {
                hoNDArray< ComplexType > aTmp(x*y);
                ComplexType* rc = aTmp.begin();

                #pragma omp for
                for ( tt=0; tt<(long long)n; tt++ )
                {
                    ComplexType* ac = a + tt*x*y;

                    size_t ay, ry;

                    for ( ay=pivoty; ay<y; ay++ )
                    {
                        ry = ay - pivoty;
                        memcpy(rc+ry*x, ac+ay*x+pivotx, sizeof(ComplexType)*(x-pivotx));
                        memcpy(rc+ry*x+x-pivotx, ac+ay*x, sizeof(ComplexType)*pivotx);
                    }

                    for ( ay=0; ay<pivoty; ay++ )
                    {
                        ry = ay + y - pivoty;
                        memcpy(rc+ry*x, ac+ay*x+pivotx, sizeof(ComplexType)*(x-pivotx));
                        memcpy(rc+ry*x+x-pivotx, ac+ay*x, sizeof(ComplexType)*pivotx);
                    }

                    memcpy(ac, rc, sizeof(ComplexType)*x*y);
                }
            }
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Errors in hoNDFFT<T>::fftshiftPivot2D(ComplexType* a, size_t x, size_t y, size_t n, size_t pivotx, size_t pivoty) ...");
            return false;
        }

        return true;
    }

    template<typename T> 
    inline bool hoNDFFT<T>::fftshift2D(const ComplexType* a, ComplexType* r, size_t x, size_t y, size_t n)
    {
        try
        {
            GADGET_CHECK_RETURN_FALSE(a!=NULL);
            GADGET_CHECK_RETURN_FALSE(r!=NULL);

            size_t pivotx = fftshiftPivot(x);
            size_t pivoty = fftshiftPivot(y);

            GADGET_CHECK_RETURN_FALSE(fftshiftPivot2D(a, r, x, y, n, pivotx, pivoty));
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Errors in hoNDFFT<T>::fftshift2D(const ComplexType* a, ComplexType* r, size_t x, size_t y, size_t n) ...");
            return false;
        }

        return true;
    }

    template<typename T> 
    inline bool hoNDFFT<T>::ifftshift2D(const ComplexType* a, ComplexType* r, size_t x, size_t y, size_t n)
    {
        try
        {
            GADGET_CHECK_RETURN_FALSE(a!=NULL);
            GADGET_CHECK_RETURN_FALSE(r!=NULL);

            size_t pivotx = ifftshiftPivot(x);
            size_t pivoty = ifftshiftPivot(y);

            GADGET_CHECK_RETURN_FALSE(fftshiftPivot2D(a, r, x, y, n, pivotx, pivoty));
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Errors in hoNDFFT<T>::ifftshift2D(const ComplexType* a, ComplexType* r, size_t x, size_t y, size_t n) ...");
            return false;
        }

        return true;
    }

    template<typename T> 
    inline bool hoNDFFT<T>::fftshift2D(ComplexType* a, size_t x, size_t y, size_t n)
    {
        try
        {
            GADGET_CHECK_RETURN_FALSE(a!=NULL);

            size_t pivotx = fftshiftPivot(x);
            size_t pivoty = fftshiftPivot(y);

            GADGET_CHECK_RETURN_FALSE(fftshiftPivot2D(a, x, y, n, pivotx, pivoty));
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Errors in hoNDFFT<T>::fftshift2D(ComplexType* a, size_t x, size_t y, size_t n) ...");
            return false;
        }

        return true;
    }

    template<typename T> 
    inline bool hoNDFFT<T>::ifftshift2D(ComplexType* a, size_t x, size_t y, size_t n)
    {
        try
        {
            GADGET_CHECK_RETURN_FALSE(a!=NULL);

            size_t pivotx = ifftshiftPivot(x);
            size_t pivoty = ifftshiftPivot(y);

            GADGET_CHECK_RETURN_FALSE(fftshiftPivot2D(a, x, y, n, pivotx, pivoty));
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Errors in hoNDFFT<T>::ifftshift2D(ComplexType* a, size_t x, size_t y, size_t n) ...");
            return false;
        }

        return true;
    }

    template<typename T> 
    inline bool hoNDFFT<T>::fftshift2D(hoNDArray< ComplexType >& a)
    {
        size_t n = a.get_number_of_elements()/(a.get_size(0)*a.get_size(1));
        return fftshift2D(a.begin(), a.get_size(0), a.get_size(1), n);
    }

    template<typename T> 
    inline bool hoNDFFT<T>::fftshift2D(const hoNDArray< ComplexType >& a, hoNDArray< ComplexType >& r)
    {
        if ( !r.dimensions_equal(&a) )
        {
            r = a;
        }

        size_t n = a.get_number_of_elements()/(a.get_size(0)*a.get_size(1));
        return fftshift2D(a.begin(), r.begin(), a.get_size(0), a.get_size(1), n);
    }

    template<typename T> 
    inline bool hoNDFFT<T>::ifftshift2D(hoNDArray< ComplexType >& a)
    {
        size_t n = a.get_number_of_elements()/(a.get_size(0)*a.get_size(1));
        return ifftshift2D(a.begin(), a.get_size(0), a.get_size(1), n);
    }

    template<typename T> 
    inline bool hoNDFFT<T>::ifftshift2D(const hoNDArray< ComplexType >& a, hoNDArray< ComplexType >& r)
    {
        if ( !r.dimensions_equal(&a) )
        {
            r = a;
        }

        size_t n = a.get_number_of_elements()/(a.get_size(0)*a.get_size(1));
        return ifftshift2D(a.begin(), r.begin(), a.get_size(0), a.get_size(1), n);
    }

    template<typename T> 
    bool hoNDFFT<T>::fftshiftPivot3D(const ComplexType* a, ComplexType* r, size_t x, size_t y, size_t z, size_t n, size_t pivotx, size_t pivoty,  size_t pivotz)
    {
        try
        {
            GADGET_CHECK_RETURN_FALSE(a!=NULL);
            GADGET_CHECK_RETURN_FALSE(r!=NULL);

            long long tt;

#pragma omp parallel for private(tt) shared(a, r, x, y, z, n, pivotx, pivoty, pivotz) if (n>1)
            for ( tt=0; tt<(long long)n; tt++ )
            {
                size_t ay, ry, az, rz;

                for ( az=pivotz; az<z; az++ )
                {
                    rz = az - pivotz;

                    const ComplexType* ac = a + tt*x*y*z + az*x*y;
                    ComplexType* rc = r + tt*x*y*z + rz*x*y;

                    for ( ay=pivoty; ay<y; ay++ )
                    {
                        ry = ay - pivoty;
                        memcpy(rc+ry*x, ac+ay*x+pivotx, sizeof(ComplexType)*(x-pivotx));
                        memcpy(rc+ry*x+x-pivotx, ac+ay*x, sizeof(ComplexType)*pivotx);
                    }

                    for ( ay=0; ay<pivoty; ay++ )
                    {
                        ry = ay + y - pivoty;
                        memcpy(rc+ry*x, ac+ay*x+pivotx, sizeof(ComplexType)*(x-pivotx));
                        memcpy(rc+ry*x+x-pivotx, ac+ay*x, sizeof(ComplexType)*pivotx);
                    }
                }

                for ( az=0; az<pivotz; az++ )
                {
                    rz = az + z - pivotz;

                    const ComplexType* ac = a + tt*x*y*z + az*x*y;
                    ComplexType* rc = r + tt*x*y*z + rz*x*y;

                    for ( ay=pivoty; ay<y; ay++ )
                    {
                        ry = ay - pivoty;
                        memcpy(rc+ry*x, ac+ay*x+pivotx, sizeof(ComplexType)*(x-pivotx));
                        memcpy(rc+ry*x+x-pivotx, ac+ay*x, sizeof(ComplexType)*pivotx);
                    }

                    for ( ay=0; ay<pivoty; ay++ )
                    {
                        ry = ay + y - pivoty;
                        memcpy(rc+ry*x, ac+ay*x+pivotx, sizeof(ComplexType)*(x-pivotx));
                        memcpy(rc+ry*x+x-pivotx, ac+ay*x, sizeof(ComplexType)*pivotx);
                    }
                }
            }
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Errors in hoNDFFT<T>::fftshiftPivot3D(const ComplexType* a, ComplexType* r, size_t x, size_t y, size_t z, size_t n, unsigned pivotx, unsigned pivoty,  unsigned pivotz) ...");
            return false;
        }

        return true;
    }

    template<typename T> 
    bool hoNDFFT<T>::fftshiftPivot3D(ComplexType* a, size_t x, size_t y, size_t z, size_t n, size_t pivotx, size_t pivoty,  size_t pivotz)
    {
        try
        {
            GADGET_CHECK_RETURN_FALSE(a!=NULL);

            long long tt;

#pragma omp parallel private(tt) shared(a, x, y, z, n, pivotx, pivoty, pivotz) if (n>1)
            {
                hoNDArray< ComplexType > aTmp(x*y*z);

#pragma omp for
                for ( tt=0; tt<(long long)n; tt++ )
                {
                    size_t ay, ry, az, rz;

                    for ( az=pivotz; az<z; az++ )
                    {
                        rz = az - pivotz;

                        const ComplexType* ac = a + tt*x*y*z + az*x*y;
                        ComplexType* rc = aTmp.begin() + rz*x*y;

                        for ( ay=pivoty; ay<y; ay++ )
                        {
                            ry = ay - pivoty;
                            memcpy(rc+ry*x, ac+ay*x+pivotx, sizeof(ComplexType)*(x-pivotx));
                            memcpy(rc+ry*x+x-pivotx, ac+ay*x, sizeof(ComplexType)*pivotx);
                        }

                        for ( ay=0; ay<pivoty; ay++ )
                        {
                            ry = ay + y - pivoty;
                            memcpy(rc+ry*x, ac+ay*x+pivotx, sizeof(ComplexType)*(x-pivotx));
                            memcpy(rc+ry*x+x-pivotx, ac+ay*x, sizeof(ComplexType)*pivotx);
                        }
                    }

                    for ( az=0; az<pivotz; az++ )
                    {
                        rz = az + z - pivotz;

                        const ComplexType* ac = a + tt*x*y*z + az*x*y;
                        ComplexType* rc = aTmp.begin() + rz*x*y;

                        for ( ay=pivoty; ay<y; ay++ )
                        {
                            ry = ay - pivoty;
                            memcpy(rc+ry*x, ac+ay*x+pivotx, sizeof(ComplexType)*(x-pivotx));
                            memcpy(rc+ry*x+x-pivotx, ac+ay*x, sizeof(ComplexType)*pivotx);
                        }

                        for ( ay=0; ay<pivoty; ay++ )
                        {
                            ry = ay + y - pivoty;
                            memcpy(rc+ry*x, ac+ay*x+pivotx, sizeof(ComplexType)*(x-pivotx));
                            memcpy(rc+ry*x+x-pivotx, ac+ay*x, sizeof(ComplexType)*pivotx);
                        }
                    }

                    memcpy(a+tt*x*y*z, aTmp.begin(), sizeof(ComplexType)*x*y*z);
                }
            }
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Errors in hoNDFFT<T>::fftshiftPivot3D(ComplexType* a, size_t x, size_t y, size_t z, size_t n, unsigned pivotx, unsigned pivoty,  unsigned pivotz) ...");
            return false;
        }

        return true;
    }

    template<typename T> 
    inline bool hoNDFFT<T>::fftshift3D(const ComplexType* a, ComplexType* r, size_t x, size_t y, size_t z, size_t n)
    {
        try
        {
            GADGET_CHECK_RETURN_FALSE(a!=NULL);
            GADGET_CHECK_RETURN_FALSE(r!=NULL);

            size_t pivotx = fftshiftPivot(x);
            size_t pivoty = fftshiftPivot(y);
            size_t pivotz = fftshiftPivot(z);

            GADGET_CHECK_RETURN_FALSE(fftshiftPivot3D(a, r, x, y, z, n, pivotx, pivoty, pivotz));
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Errors in hoNDFFT<T>::fftshift3D(const ComplexType* a, ComplexType* r, size_t x, size_t y, size_t z, size_t n) ...");
            return false;
        }

        return true;
    }

    template<typename T> 
    inline bool hoNDFFT<T>::ifftshift3D(const ComplexType* a, ComplexType* r, size_t x, size_t y, size_t z, size_t n)
    {
        try
        {
            GADGET_CHECK_RETURN_FALSE(a!=NULL);
            GADGET_CHECK_RETURN_FALSE(r!=NULL);

            size_t pivotx = ifftshiftPivot(x);
            size_t pivoty = ifftshiftPivot(y);
            size_t pivotz = ifftshiftPivot(z);

            GADGET_CHECK_RETURN_FALSE(fftshiftPivot3D(a, r, x, y, z, n, pivotx, pivoty, pivotz));
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Errors in hoNDFFT<T>::ifftshift3D(const ComplexType* a, ComplexType* r, size_t x, size_t y, size_t z, size_t n) ...");
            return false;
        }

        return true;
    }

    template<typename T> 
    inline bool hoNDFFT<T>::fftshift3D(ComplexType* a, size_t x, size_t y, size_t z, size_t n)
    {
        try
        {
            GADGET_CHECK_RETURN_FALSE(a!=NULL);

            size_t pivotx = fftshiftPivot(x);
            size_t pivoty = fftshiftPivot(y);
            size_t pivotz = fftshiftPivot(z);

            GADGET_CHECK_RETURN_FALSE(fftshiftPivot3D(a, x, y, z, n, pivotx, pivoty, pivotz));
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Errors in hoNDFFT<T>::fftshift3D(ComplexType* a, size_t x, size_t y, size_t z, size_t n) ...");
            return false;
        }

        return true;
    }

    template<typename T> 
    inline bool hoNDFFT<T>::ifftshift3D(ComplexType* a, size_t x, size_t y, size_t z, size_t n)
    {
        try
        {
            GADGET_CHECK_RETURN_FALSE(a!=NULL);

            size_t pivotx = ifftshiftPivot(x);
            size_t pivoty = ifftshiftPivot(y);
            size_t pivotz = ifftshiftPivot(z);

            GADGET_CHECK_RETURN_FALSE(fftshiftPivot3D(a, x, y, z, n, pivotx, pivoty, pivotz));
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Errors in hoNDFFT<T>::ifftshift3D(ComplexType* a, size_t x, size_t y, size_t z, size_t n) ...");
            return false;
        }

        return true;
    }

    template<typename T> 
    inline bool hoNDFFT<T>::fftshift3D(hoNDArray< ComplexType >& a)
    {
        size_t n = a.get_number_of_elements()/(a.get_size(0)*a.get_size(1)*a.get_size(2));
        return fftshift3D(a.begin(), a.get_size(0), a.get_size(1), a.get_size(2), n);
    }

    template<typename T> 
    inline bool hoNDFFT<T>::fftshift3D(const hoNDArray< ComplexType >& a, hoNDArray< ComplexType >& r)
    {
        if ( !r.dimensions_equal(&a) )
        {
            r = a;
        }

        size_t n = a.get_number_of_elements()/(a.get_size(0)*a.get_size(1)*a.get_size(2));
        return fftshift3D(a.begin(), r.begin(), a.get_size(0), a.get_size(1), a.get_size(2), n);
    }

    template<typename T> 
    inline bool hoNDFFT<T>::ifftshift3D(hoNDArray< ComplexType >& a)
    {
        size_t n = a.get_number_of_elements()/(a.get_size(0)*a.get_size(1)*a.get_size(2));
        return ifftshift3D(a.begin(), a.get_size(0), a.get_size(1), a.get_size(2), n);
    }

    template<typename T> 
    inline bool hoNDFFT<T>::ifftshift3D(const hoNDArray< ComplexType >& a, hoNDArray< ComplexType >& r)
    {
        if ( !r.dimensions_equal(&a) )
        {
            r = a;
        }

        size_t n = a.get_number_of_elements()/(a.get_size(0)*a.get_size(1)*a.get_size(2));
        return ifftshift3D(a.begin(), r.begin(), a.get_size(0), a.get_size(1), a.get_size(2), n);
    }

    // -----------------------------------------------------------------------------------------

    template<typename T> 
    inline bool hoNDFFT<T>::fft1(hoNDArray< ComplexType >& a)
    {
        return fft1(a, true);
    }

    template<typename T> 
    inline bool hoNDFFT<T>::ifft1(hoNDArray< ComplexType >& a)
    {
        return fft1(a, false);
    }

    template<typename T> 
    inline bool hoNDFFT<T>::fft1(const hoNDArray< ComplexType >& a, hoNDArray< ComplexType >& r)
    {
        if ( !r.dimensions_equal(&a) )
        {
            r.create(a.get_dimensions());
        }

        return fft1(const_cast<hoNDArray< ComplexType >&>(a), r, true);
    }

    template<typename T> 
    inline bool hoNDFFT<T>::ifft1(const hoNDArray< ComplexType >& a, hoNDArray< ComplexType >& r)
    {
        if ( !r.dimensions_equal(&a) )
        {
            r.create(a.get_dimensions());
        }

        return fft1(const_cast<hoNDArray< ComplexType >&>(a), r, false);
    }

    template<typename T> 
    inline bool hoNDFFT<T>::fft1c(hoNDArray< ComplexType >& a)
    {
        GADGET_CHECK_RETURN_FALSE(ifftshift1D(a));
        GADGET_CHECK_RETURN_FALSE(fft1(a));
        GADGET_CHECK_RETURN_FALSE(fftshift1D(a));
        return true;
    }

    template<typename T> 
    inline bool hoNDFFT<T>::ifft1c(hoNDArray< ComplexType >& a)
    {
        GADGET_CHECK_RETURN_FALSE(ifftshift1D(a));
        GADGET_CHECK_RETURN_FALSE(ifft1(a));
        GADGET_CHECK_RETURN_FALSE(fftshift1D(a));
        return true;
    }

    template<typename T> 
    inline bool hoNDFFT<T>::fft1c(const hoNDArray< ComplexType >& a, hoNDArray< ComplexType >& r)
    {
        GADGET_CHECK_RETURN_FALSE(ifftshift1D(a, r));
        GADGET_CHECK_RETURN_FALSE(fft1(r));
        GADGET_CHECK_RETURN_FALSE(fftshift1D(r));
        return true;
    }

    template<typename T> 
    inline bool hoNDFFT<T>::ifft1c(const hoNDArray< ComplexType >& a, hoNDArray< ComplexType >& r)
    {
        GADGET_CHECK_RETURN_FALSE(ifftshift1D(a, r));
        GADGET_CHECK_RETURN_FALSE(ifft1(r));
        GADGET_CHECK_RETURN_FALSE(fftshift1D(r));
        return true;
    }

    template<typename T> 
    inline bool hoNDFFT<T>::fft1c(const hoNDArray< ComplexType >& a, hoNDArray< ComplexType >& r, hoNDArray< ComplexType >& buf)
    {
        GADGET_CHECK_RETURN_FALSE(ifftshift1D(a, r));
        GADGET_CHECK_RETURN_FALSE(fft1(r, buf));
        GADGET_CHECK_RETURN_FALSE(fftshift1D(buf, r));
        return true;
    }

    template<typename T> 
    inline bool hoNDFFT<T>::ifft1c(const hoNDArray< ComplexType >& a, hoNDArray< ComplexType >& r, hoNDArray< ComplexType >& buf)
    {
        GADGET_CHECK_RETURN_FALSE(ifftshift1D(a, r));
        GADGET_CHECK_RETURN_FALSE(ifft1(r, buf));
        GADGET_CHECK_RETURN_FALSE(fftshift1D(buf, r));
        return true;
    }

    // -----------------------------------------------------------------------------------------

    template<typename T> 
    inline bool hoNDFFT<T>::fft2(hoNDArray< ComplexType >& a)
    {
        return fft2(a, true);
    }

    template<typename T> 
    inline bool hoNDFFT<T>::ifft2(hoNDArray< ComplexType >& a)
    {
        return fft2(a, false);
    }

    template<typename T> 
    inline bool hoNDFFT<T>::fft2(const hoNDArray< ComplexType >& a, hoNDArray< ComplexType >& r)
    {
        //r = a;
        //return fft2(r);
        if ( !r.dimensions_equal(&a) )
        {
            r.create(a.get_dimensions());
        }

        return fft2(const_cast<hoNDArray< ComplexType >&>(a), r, true);
    }

    template<typename T> 
    inline bool hoNDFFT<T>::ifft2(const hoNDArray< ComplexType >& a, hoNDArray< ComplexType >& r)
    {
        /*r = a;
        return ifft2(r);*/

        if ( !r.dimensions_equal(&a) )
        {
            r.create(a.get_dimensions());
        }

        return fft2(const_cast<hoNDArray< ComplexType >&>(a), r, false);
    }

    template<typename T> 
    inline bool hoNDFFT<T>::fft2c(hoNDArray< ComplexType >& a)
    {
        GADGET_CHECK_RETURN_FALSE(ifftshift2D(a));
        GADGET_CHECK_RETURN_FALSE(fft2(a));
        GADGET_CHECK_RETURN_FALSE(fftshift2D(a));
        return true;
    }

    template<typename T> 
    inline bool hoNDFFT<T>::ifft2c(hoNDArray< ComplexType >& a)
    {
        GADGET_CHECK_RETURN_FALSE(ifftshift2D(a));
        GADGET_CHECK_RETURN_FALSE(ifft2(a));
        GADGET_CHECK_RETURN_FALSE(fftshift2D(a));
        return true;
    }

    template<typename T> 
    inline bool hoNDFFT<T>::fft2c(const hoNDArray< ComplexType >& a, hoNDArray< ComplexType >& r)
    {
        GADGET_CHECK_RETURN_FALSE(ifftshift2D(a, r));
        GADGET_CHECK_RETURN_FALSE(fft2(r));
        GADGET_CHECK_RETURN_FALSE(fftshift2D(r));
        return true;
    }

    template<typename T> 
    inline bool hoNDFFT<T>::ifft2c(const hoNDArray< ComplexType >& a, hoNDArray< ComplexType >& r)
    {
        GADGET_CHECK_RETURN_FALSE(ifftshift2D(a, r));
        GADGET_CHECK_RETURN_FALSE(ifft2(r));
        GADGET_CHECK_RETURN_FALSE(fftshift2D(r));
        return true;
    }

    template<typename T> 
    inline bool hoNDFFT<T>::fft2c(const hoNDArray< ComplexType >& a, hoNDArray< ComplexType >& r, hoNDArray< ComplexType >& buf)
    {
        GADGET_CHECK_RETURN_FALSE(ifftshift2D(a, r));
        GADGET_CHECK_RETURN_FALSE(fft2(r, buf));
        GADGET_CHECK_RETURN_FALSE(fftshift2D(buf, r));
        return true;
    }

    template<typename T> 
    inline bool hoNDFFT<T>::ifft2c(const hoNDArray< ComplexType >& a, hoNDArray< ComplexType >& r, hoNDArray< ComplexType >& buf)
    {
        GADGET_CHECK_RETURN_FALSE(ifftshift2D(a, r));
        GADGET_CHECK_RETURN_FALSE(ifft2(r, buf));
        GADGET_CHECK_RETURN_FALSE(fftshift2D(buf, r));
        return true;
    }

    // -----------------------------------------------------------------------------------------

    template<typename T> 
    inline bool hoNDFFT<T>::fft3(hoNDArray< ComplexType >& a)
    {
        return fft3(a, true);
    }

    template<typename T> 
    inline bool hoNDFFT<T>::ifft3(hoNDArray< ComplexType >& a)
    {
        return fft3(a, false);
    }

    template<typename T> 
    inline bool hoNDFFT<T>::fft3(const hoNDArray< ComplexType >& a, hoNDArray< ComplexType >& r)
    {
        /*r = a;
        return fft3(r);*/
        if ( !r.dimensions_equal(&a) )
        {
            r.create(a.get_dimensions());
        }

        return fft3(const_cast<hoNDArray< ComplexType >&>(a), r, true);
    }

    template<typename T> 
    inline bool hoNDFFT<T>::ifft3(const hoNDArray< ComplexType >& a, hoNDArray< ComplexType >& r)
    {
        /*r = a;
        return ifft3(r);*/
        if ( !r.dimensions_equal(&a) )
        {
            r.create(a.get_dimensions());
        }

        return fft3(const_cast<hoNDArray< ComplexType >&>(a), r, false);
    }

    template<typename T> 
    inline bool hoNDFFT<T>::fft3c(hoNDArray< ComplexType >& a)
    {
        GADGET_CHECK_RETURN_FALSE(ifftshift3D(a));
        GADGET_CHECK_RETURN_FALSE(fft3(a));
        GADGET_CHECK_RETURN_FALSE(fftshift3D(a));
        return true;
    }

    template<typename T> 
    inline bool hoNDFFT<T>::ifft3c(hoNDArray< ComplexType >& a)
    {
        GADGET_CHECK_RETURN_FALSE(ifftshift3D(a));
        GADGET_CHECK_RETURN_FALSE(ifft3(a));
        GADGET_CHECK_RETURN_FALSE(fftshift3D(a));
        return true;
    }

    template<typename T> 
    inline bool hoNDFFT<T>::fft3c(const hoNDArray< ComplexType >& a, hoNDArray< ComplexType >& r)
    {
        GADGET_CHECK_RETURN_FALSE(ifftshift3D(a, r));
        GADGET_CHECK_RETURN_FALSE(fft3(r));
        GADGET_CHECK_RETURN_FALSE(fftshift3D(r));
        return true;
    }

    template<typename T> 
    inline bool hoNDFFT<T>::ifft3c(const hoNDArray< ComplexType >& a, hoNDArray< ComplexType >& r)
    {
        GADGET_CHECK_RETURN_FALSE(ifftshift3D(a, r));
        GADGET_CHECK_RETURN_FALSE(ifft3(r));
        GADGET_CHECK_RETURN_FALSE(fftshift3D(r));
        return true;
    }

    template<typename T> 
    inline bool hoNDFFT<T>::fft3c(const hoNDArray< ComplexType >& a, hoNDArray< ComplexType >& r, hoNDArray< ComplexType >& buf)
    {
        GADGET_CHECK_RETURN_FALSE(ifftshift3D(a, r));
        GADGET_CHECK_RETURN_FALSE(fft3(r, buf));
        GADGET_CHECK_RETURN_FALSE(fftshift3D(buf, r));
        return true;
    }

    template<typename T> 
    inline bool hoNDFFT<T>::ifft3c(const hoNDArray< ComplexType >& a, hoNDArray< ComplexType >& r, hoNDArray< ComplexType >& buf)
    {
        GADGET_CHECK_RETURN_FALSE(ifftshift3D(a, r));
        GADGET_CHECK_RETURN_FALSE(ifft3(r, buf));
        GADGET_CHECK_RETURN_FALSE(fftshift3D(buf, r));
        return true;
    }

    template<typename T> 
    bool hoNDFFT<T>::fft1(hoNDArray< ComplexType >& a, bool forward)
    {
#ifdef USE_MKL
        return fft1_mkl(a, forward);
#else
        hoNDArray< ComplexType > res(a);
        if ( !fft1(res, a, forward) )
        {
            return false;
        }

        return true;
#endif //USE_MKL
    }

    template<typename T> 
    bool hoNDFFT<T>::fft2(hoNDArray< ComplexType >& a, bool forward)
    {
#ifdef USE_MKL
        return fft2_mkl(a, forward);
#else
        hoNDArray< ComplexType > res(a);
        if ( !fft2(res, a, forward) )
        {
            return false;
        }

        return true;
#endif //USE_MKL
    }

    template<typename T> 
    bool hoNDFFT<T>::fft3(hoNDArray< ComplexType >& a, bool forward)
    {
#ifdef USE_MKL
        return fft3_mkl(a, forward);
#else
        hoNDArray< ComplexType > res(a);
        if ( !fft3(res, a, forward) )
        {
            return false;
        }

        return true;
#endif //USE_MKL
    }

    template<typename T> 
    bool hoNDFFT<T>::fft1(hoNDArray< ComplexType >& a, hoNDArray< ComplexType >& r, bool forward)
    {
#ifdef USE_MKL
        return fft1_mkl(a, r, forward);
#else
        r = a;

        int n0 = a.get_size(0);
        T fftRatio = 1.0/std::sqrt( T(n0) );

        size_t num = a.get_number_of_elements()/n0;
        long long n;

        if ( typeid(T) == typeid(float) )
        {
            fftwf_plan p;

            {
                mutex_.lock();
                if ( forward )
                {
                    p = fftwf_plan_dft_1d(n0, 
                            reinterpret_cast<fftwf_complex*>(a.begin()), 
                            reinterpret_cast<fftwf_complex*>(r.begin()),
                            FFTW_FORWARD, FFTW_ESTIMATE);
                }
                else
                {
                    p = fftwf_plan_dft_1d(n0, 
                            reinterpret_cast<fftwf_complex*>(a.begin()), 
                            reinterpret_cast<fftwf_complex*>(r.begin()),
                            FFTW_BACKWARD, FFTW_ESTIMATE);
                }
                mutex_.unlock();
            }

            #pragma omp parallel for private(n) shared(num, p, a, n0, r)
            for ( n=0; n<num; n++ )
            {
                fftwf_execute_dft(p, reinterpret_cast<fftwf_complex*>(a.begin()+n*n0), 
                    reinterpret_cast<fftwf_complex*>(r.begin()+n*n0));
            }

            {
                mutex_.lock();
                fftwf_destroy_plan(p);
                mutex_.unlock();
            }
        }
        else if ( typeid(T) == typeid(double) )
        {
            fftw_plan p;

            {
                mutex_.lock();
                if ( forward )
                {
                    p = fftw_plan_dft_1d(n0, 
                            reinterpret_cast<fftw_complex*>(a.begin()), 
                            reinterpret_cast<fftw_complex*>(r.begin()),
                            FFTW_FORWARD, FFTW_ESTIMATE);
                }
                else
                {
                    p = fftw_plan_dft_1d(n0, 
                            reinterpret_cast<fftw_complex*>(a.begin()), 
                            reinterpret_cast<fftw_complex*>(r.begin()),
                            FFTW_BACKWARD, FFTW_ESTIMATE);
                }
                mutex_.unlock();
            }

            #pragma omp parallel for private(n) shared(num, p, a, n0, r)
            for ( n=0; n<num; n++ )
            {
                fftw_execute_dft(p, reinterpret_cast<fftw_complex*>(a.begin()+n*n0), 
                    reinterpret_cast<fftw_complex*>(r.begin()+n*n0));
            }

            {
                mutex_.lock();
                fftw_destroy_plan(p);
                mutex_.unlock();
            }
        }

        Gadgetron::scal(fftRatio, r);

        return true;
#endif //USE_MKL
    }

    template<typename T> 
    bool hoNDFFT<T>::fft2(hoNDArray< ComplexType >& a, hoNDArray< ComplexType >& r, bool forward)
    {
#ifdef USE_MKL
        return fft2_mkl(a, r, forward);
#else
        r = a;

        int n0 = a.get_size(1);
        int n1 = a.get_size(0);

        T fftRatio = 1.0/std::sqrt( T(n0*n1) );

        size_t num = a.get_number_of_elements()/(n0*n1);
        long long n;

        if ( typeid(T) == typeid(float) )
        {
            fftwf_plan p;

            {
                mutex_.lock();
                if ( forward )
                {
                    p = fftwf_plan_dft_2d(n0, n1,
                            reinterpret_cast<fftwf_complex*>(a.begin()), 
                            reinterpret_cast<fftwf_complex*>(r.begin()),
                            FFTW_FORWARD, FFTW_ESTIMATE);
                }
                else
                {
                    p = fftwf_plan_dft_2d(n0, n1,
                            reinterpret_cast<fftwf_complex*>(a.begin()), 
                            reinterpret_cast<fftwf_complex*>(r.begin()),
                            FFTW_BACKWARD, FFTW_ESTIMATE);
                }
                mutex_.unlock();
            }

            #pragma omp parallel for private(n) shared(num, p, a, n0, n1, r)
            for ( n=0; n<num; n++ )
            {
                fftwf_execute_dft(p, reinterpret_cast<fftwf_complex*>(a.begin()+n*n0*n1), 
                    reinterpret_cast<fftwf_complex*>(r.begin()+n*n0*n1));
            }

            {
                mutex_.lock();
                fftwf_destroy_plan(p);
                mutex_.unlock();
            }
        }
        else if ( typeid(T) == typeid(double) )
        {
            fftw_plan p;

            {
                mutex_.lock();
                if ( forward )
                {
                    p = fftw_plan_dft_2d(n0, n1,
                            reinterpret_cast<fftw_complex*>(a.begin()), 
                            reinterpret_cast<fftw_complex*>(r.begin()),
                            FFTW_FORWARD, FFTW_ESTIMATE);
                }
                else
                {
                    p = fftw_plan_dft_2d(n0, n1,
                            reinterpret_cast<fftw_complex*>(a.begin()), 
                            reinterpret_cast<fftw_complex*>(r.begin()),
                            FFTW_BACKWARD, FFTW_ESTIMATE);
                }
                mutex_.unlock();
            }

            #pragma omp parallel for private(n) shared(num, p, a, n0, n1, r)
            for ( n=0; n<num; n++ )
            {
                fftw_execute_dft(p, reinterpret_cast<fftw_complex*>(a.begin()+n*n0*n1), 
                    reinterpret_cast<fftw_complex*>(r.begin()+n*n0*n1));
            }

            {
                mutex_.lock();
                fftw_destroy_plan(p);
                mutex_.unlock();
            }
        }

        Gadgetron::scal(fftRatio, r);

        return true;
#endif //USE_MKL
    }

    template<typename T> 
    bool hoNDFFT<T>::fft3(hoNDArray< ComplexType >& a, hoNDArray< ComplexType >& r, bool forward)
    {
#ifdef USE_MKL
        return fft3_mkl(a, r, forward);
#else
        r = a;

        int n2 = a.get_size(0);
        int n1 = a.get_size(1);
        int n0 = a.get_size(2);

        T fftRatio = 1.0/std::sqrt( T(n0*n1*n2) );

        size_t num = a.get_number_of_elements()/(n0*n1*n2);
        long long n;

        if ( typeid(T) == typeid(float) )
        {
            fftwf_plan p;

            {
                mutex_.lock();
                if ( forward )
                {
                    p = fftwf_plan_dft_3d(n0, n1, n2, 
                            reinterpret_cast<fftwf_complex*>(a.begin()), 
                            reinterpret_cast<fftwf_complex*>(r.begin()),
                            FFTW_FORWARD, FFTW_ESTIMATE);
                }
                else
                {
                    p = fftwf_plan_dft_3d(n0, n1, n2, 
                            reinterpret_cast<fftwf_complex*>(a.begin()), 
                            reinterpret_cast<fftwf_complex*>(r.begin()),
                            FFTW_BACKWARD, FFTW_ESTIMATE);
                }
                mutex_.unlock();
            }

            #pragma omp parallel for private(n) shared(num, p, a, n0, n1, n2, r) if (num > 8)
            for ( n=0; n<num; n++ )
            {
                fftwf_execute_dft(p, reinterpret_cast<fftwf_complex*>(a.begin()+n*n0*n1*n2), 
                    reinterpret_cast<fftwf_complex*>(r.begin()+n*n0*n1*n2));
            }

            {
                mutex_.lock();
                fftwf_destroy_plan(p);
                mutex_.unlock();
            }
        }
        else if ( typeid(T) == typeid(double) )
        {
            fftw_plan p;

            {
                mutex_.lock();
                if ( forward )
                {
                    p = fftw_plan_dft_3d(n0, n1, n2, 
                            reinterpret_cast<fftw_complex*>(a.begin()), 
                            reinterpret_cast<fftw_complex*>(r.begin()),
                            FFTW_FORWARD, FFTW_ESTIMATE);
                }
                else
                {
                    p = fftw_plan_dft_3d(n0, n1, n2, 
                            reinterpret_cast<fftw_complex*>(a.begin()), 
                            reinterpret_cast<fftw_complex*>(r.begin()),
                            FFTW_BACKWARD, FFTW_ESTIMATE);
                }
                mutex_.unlock();
            }

            #pragma omp parallel for private(n) shared(num, p, a, n0, n1, n2, r) if (num > 8)
            for ( n=0; n<num; n++ )
            {
                fftw_execute_dft(p, reinterpret_cast<fftw_complex*>(a.begin()+n*n0*n1*n2), 
                    reinterpret_cast<fftw_complex*>(r.begin()+n*n0*n1*n2));
            }

            {
                mutex_.lock();
                fftw_destroy_plan(p);
                mutex_.unlock();
            }
        }

        Gadgetron::scal(fftRatio, r);

        return true;
#endif //USE_MKL
    }

    // -----------------------------------------------------------------------------------------

    // MKL related

#ifdef USE_MKL

    template<typename T> 
    bool hoNDFFT<T>::configureFFTHandle(long long NDim, MKL_LONG* dim, DFTI_CONFIG_VALUE fftPresion, size_t n, DFTI_DESCRIPTOR_HANDLE& handle)
    {
        long long ii;

        MKL_LONG res;

        if ( NDim == 1 )
        {
            if ( (res=DftiCreateDescriptor( &handle, fftPresion, DFTI_COMPLEX, NDim, dim[0])) != 0 )
            {
                GADGET_ERROR_MSG( DftiErrorMessage(res) );
                return false;
            }
        }
        else
        {
            if ( (res=DftiCreateDescriptor( &handle, fftPresion, DFTI_COMPLEX, NDim, dim)) != 0 )
            {
                GADGET_ERROR_MSG( DftiErrorMessage(res) );
                return false;
            }
        }

        double fftScaling = 1.0;
        for ( ii=0; ii<NDim; ii++ )
        {
            fftScaling *= dim[ii];
        }

        if ( (res=DftiSetValue( handle, DFTI_FORWARD_SCALE, 1.0/std::sqrt(fftScaling))) != 0 )
        {
            GADGET_ERROR_MSG( DftiErrorMessage(res) );
            return false;
        }

        if ( (res=DftiSetValue( handle, DFTI_BACKWARD_SCALE, 1.0/std::sqrt(fftScaling))) != 0 )
        {
            GADGET_ERROR_MSG( DftiErrorMessage(res) );
            return false;
        }

        if ( (res=DftiSetValue( handle, DFTI_PLACEMENT, DFTI_INPLACE)) != 0 )
        {
            GADGET_ERROR_MSG( DftiErrorMessage(res) );
            return false;
        }

        if ( n > 1 )
        {
            if ( (res=DftiSetValue( handle, DFTI_NUMBER_OF_TRANSFORMS, n)) != 0 )
            {
                GADGET_ERROR_MSG( DftiErrorMessage(res) );
                return false;
            }

            if ( (res=DftiSetValue( handle, DFTI_INPUT_DISTANCE, (MKL_INT)fftScaling)) != 0 )
            {
                GADGET_ERROR_MSG( DftiErrorMessage(res) );
                return false;
            }

            if ( (res=DftiSetValue( handle, DFTI_OUTPUT_DISTANCE, (MKL_INT)fftScaling)) != 0 )
            {
                GADGET_ERROR_MSG( DftiErrorMessage(res) );
                return false;
            }
        }

        if ( (res=DftiCommitDescriptor( handle)) != 0 )
        {
            GADGET_ERROR_MSG( DftiErrorMessage(res) );
            return false;
        }

        return true;
    }

    template<typename T> 
    bool hoNDFFT<T>::configureFFTHandleOutOfPlace(long long NDim, MKL_LONG* dim, DFTI_CONFIG_VALUE fftPresion, size_t n, DFTI_DESCRIPTOR_HANDLE& handle)
    {
        long long ii;

        MKL_LONG res;

        if ( NDim == 1 )
        {
            if ( (res=DftiCreateDescriptor( &handle, fftPresion, DFTI_COMPLEX, NDim, dim[0])) != 0 )
            {
                GADGET_ERROR_MSG( DftiErrorMessage(res) );
                return false;
            }
        }
        else
        {
            if ( (res=DftiCreateDescriptor( &handle, fftPresion, DFTI_COMPLEX, NDim, dim)) != 0 )
            {
                GADGET_ERROR_MSG( DftiErrorMessage(res) );
                return false;
            }
        }

        double fftScaling = 1.0;
        for ( ii=0; ii<NDim; ii++ )
        {
            fftScaling *= dim[ii];
        }

        if ( (res=DftiSetValue( handle, DFTI_FORWARD_SCALE, 1.0/std::sqrt(fftScaling))) != 0 )
        {
            GADGET_ERROR_MSG( DftiErrorMessage(res) );
            return false;
        }

        if ( (res=DftiSetValue( handle, DFTI_BACKWARD_SCALE, 1.0/std::sqrt(fftScaling))) != 0 )
        {
            GADGET_ERROR_MSG( DftiErrorMessage(res) );
            return false;
        }

        if ( (res=DftiSetValue( handle, DFTI_PLACEMENT, DFTI_NOT_INPLACE)) != 0 )
        {
            GADGET_ERROR_MSG( DftiErrorMessage(res) );
            return false;
        }

        if ( n > 1 )
        {
            if ( (res=DftiSetValue( handle, DFTI_NUMBER_OF_TRANSFORMS, n)) != 0 )
            {
                GADGET_ERROR_MSG( DftiErrorMessage(res) );
                return false;
            }

            if ( (res=DftiSetValue( handle, DFTI_INPUT_DISTANCE, (MKL_INT)fftScaling)) != 0 )
            {
                GADGET_ERROR_MSG( DftiErrorMessage(res) );
                return false;
            }

            if ( (res=DftiSetValue( handle, DFTI_OUTPUT_DISTANCE, (MKL_INT)fftScaling)) != 0 )
            {
                GADGET_ERROR_MSG( DftiErrorMessage(res) );
                return false;
            }
        }

        if ( (res=DftiCommitDescriptor( handle)) != 0 )
        {
            GADGET_ERROR_MSG( DftiErrorMessage(res) );
            return false;
        }

        return true;
    }

    template<typename T> 
    bool hoNDFFT<T>::fft1_mkl(hoNDArray< ComplexType >& a, bool forward)
    {
        size_t n = a.get_number_of_elements()/a.get_size(0);
        MKL_LONG dim = a.get_size(0);

        DFTI_DESCRIPTOR_HANDLE handle;

        if ( typeid(T) == typeid(float) )
        {
            GADGET_CHECK_RETURN_FALSE(configureFFTHandle(1, &dim, DFTI_SINGLE, n, handle));
        }
        else if ( typeid(T) == typeid(double) )
        {
            GADGET_CHECK_RETURN_FALSE(configureFFTHandle(1, &dim, DFTI_DOUBLE, n, handle));
        }
        else
        {
            GADGET_ERROR_MSG("hoNDFFT<T>::fft1_mkl(hoNDArray< ComplexType >& a), only float and double are supported ... ");
            return false;
        }

        MKL_LONG res;

        if ( forward )
        {
            if ( ( res=DftiComputeForward(handle, reinterpret_cast<T*>(a.begin())) ) != 0 ) 
            { 
                GADGET_ERROR_MSG( DftiErrorMessage(res) ); 
                return false; 
            }
        }
        else
        {
            if ( ( res=DftiComputeBackward(handle, reinterpret_cast<T*>(a.begin())) ) != 0 ) 
            { 
                GADGET_ERROR_MSG( DftiErrorMessage(res) ); 
                return false; 
            }
        }

        if ( ( res=DftiFreeDescriptor(&handle) ) != 0 ) 
        { 
            GADGET_ERROR_MSG( DftiErrorMessage(res) ); 
            return false; 
        }

        return true;
    }

    template<typename T> 
    bool hoNDFFT<T>::fft1_mkl(hoNDArray< ComplexType >& a, hoNDArray< ComplexType >& r, bool forward)
    {
        size_t n = a.get_number_of_elements()/a.get_size(0);
        MKL_LONG dim = a.get_size(0);

        DFTI_DESCRIPTOR_HANDLE handle;

        if ( typeid(T) == typeid(float) )
        {
            GADGET_CHECK_RETURN_FALSE(configureFFTHandleOutOfPlace(1, &dim, DFTI_SINGLE, n, handle));
        }
        else if ( typeid(T) == typeid(double) )
        {
            GADGET_CHECK_RETURN_FALSE(configureFFTHandleOutOfPlace(1, &dim, DFTI_DOUBLE, n, handle));
        }
        else
        {
            GADGET_ERROR_MSG("hoNDFFT<T>::fft1_mkl(hoNDArray< ComplexType >& a, hoNDArray< ComplexType >& r), only float and double are supported ... ");
            return false;
        }

        MKL_LONG res;

        if ( forward )
        {
            if ( ( res=DftiComputeForward( handle, reinterpret_cast<T*>(a.begin()), reinterpret_cast<T*>(r.begin()) ) ) != 0 ) 
            { 
                GADGET_ERROR_MSG( DftiErrorMessage(res) ); 
                return false; 
            }
        }
        else
        {
            if ( ( res=DftiComputeBackward( handle, reinterpret_cast<T*>(a.begin()), reinterpret_cast<T*>(r.begin()) ) ) != 0 ) 
            { 
                GADGET_ERROR_MSG( DftiErrorMessage(res) ); 
                return false; 
            }
        }

        if ( ( res=DftiFreeDescriptor(&handle) ) != 0 ) 
        { 
            GADGET_ERROR_MSG( DftiErrorMessage(res) ); 
            return false; 
        }

        return true;
    }

    template<typename T> 
    bool hoNDFFT<T>::fft2_mkl(hoNDArray< ComplexType >& a, bool forward)
    {
        size_t n = a.get_number_of_elements()/(a.get_size(0)*a.get_size(1));
        MKL_LONG dim[2];
        dim[0] = a.get_size(1);
        dim[1] = a.get_size(0);

        DFTI_DESCRIPTOR_HANDLE handle;

        if ( typeid(T) == typeid(float) )
        {
            GADGET_CHECK_RETURN_FALSE(configureFFTHandle(2, dim, DFTI_SINGLE, n, handle));
        }
        else if ( typeid(T) == typeid(double) )
        {
            GADGET_CHECK_RETURN_FALSE(configureFFTHandle(2, dim, DFTI_DOUBLE, n, handle));
        }
        else
        {
            GADGET_ERROR_MSG("hoNDFFT<T>::fft2_mkl(hoNDArray< ComplexType >& a), only float and double are supported ... ");
            return false;
        }

        MKL_LONG res;
        if ( forward )
        {
            if ( ( res=DftiComputeForward(handle, reinterpret_cast<T*>(a.begin())) ) != 0 ) 
            { 
                GADGET_ERROR_MSG( DftiErrorMessage(res) ); 
                return false; 
            }
        }
        else
        {
            if ( ( res=DftiComputeBackward(handle, reinterpret_cast<T*>(a.begin())) ) != 0 ) 
            { 
                GADGET_ERROR_MSG( DftiErrorMessage(res) ); 
                return false; 
            }
        }

        if ( ( res=DftiFreeDescriptor(&handle) ) != 0 ) 
        { 
            GADGET_ERROR_MSG( DftiErrorMessage(res) ); 
            return false; 
        }

        return true;
    }

    template<typename T> 
    bool hoNDFFT<T>::fft2_mkl(hoNDArray< ComplexType >& a, hoNDArray< ComplexType >& r, bool forward)
    {
        size_t n = a.get_number_of_elements()/(a.get_size(0)*a.get_size(1));
        MKL_LONG dim[2];
        dim[0] = a.get_size(1);
        dim[1] = a.get_size(0);

        DFTI_DESCRIPTOR_HANDLE handle;

        if ( typeid(T) == typeid(float) )
        {
            GADGET_CHECK_RETURN_FALSE(configureFFTHandleOutOfPlace(2, dim, DFTI_SINGLE, n, handle));
        }
        else if ( typeid(T) == typeid(double) )
        {
            GADGET_CHECK_RETURN_FALSE(configureFFTHandleOutOfPlace(2, dim, DFTI_DOUBLE, n, handle));
        }
        else
        {
            GADGET_ERROR_MSG("hoNDFFT<T>::fft2_mkl(hoNDArray< ComplexType >& a, hoNDArray< ComplexType >& r), only float and double are supported ... ");
            return false;
        }

        MKL_LONG res;
        if ( forward )
        {
            if ( ( res=DftiComputeForward(handle, reinterpret_cast<T*>(a.begin()), reinterpret_cast<T*>(r.begin())) ) != 0 ) 
            { 
                GADGET_ERROR_MSG( DftiErrorMessage(res) ); 
                return false; 
            }
        }
        else
        {
            if ( ( res=DftiComputeBackward(handle, reinterpret_cast<T*>(a.begin()), reinterpret_cast<T*>(r.begin())) ) != 0 ) 
            { 
                GADGET_ERROR_MSG( DftiErrorMessage(res) ); 
                return false; 
            }
        }

        if ( ( res=DftiFreeDescriptor(&handle) ) != 0 ) 
        { 
            GADGET_ERROR_MSG( DftiErrorMessage(res) ); 
            return false; 
        }

        return true;
    }

    template<typename T> 
    bool hoNDFFT<T>::fft3_mkl(hoNDArray< ComplexType >& a, bool forward)
    {
        size_t n = a.get_number_of_elements()/(a.get_size(0)*a.get_size(1)*a.get_size(2));

        MKL_LONG dim[3];
        dim[0] = a.get_size(2);
        dim[1] = a.get_size(1);
        dim[2] = a.get_size(0);

        DFTI_DESCRIPTOR_HANDLE handle;

        if ( typeid(T) == typeid(float) )
        {
            GADGET_CHECK_RETURN_FALSE(configureFFTHandle(3, dim, DFTI_SINGLE, n, handle));
        }
        else if ( typeid(T) == typeid(double) )
        {
            GADGET_CHECK_RETURN_FALSE(configureFFTHandle(3, dim, DFTI_DOUBLE, n, handle));
        }
        else
        {
            GADGET_ERROR_MSG("hoNDFFT<T>::fft3_mkl(hoNDArray< ComplexType >& a), only float and double are supported ... ");
            return false;
        }

        MKL_LONG res;
        if ( forward )
        {
            if ( ( res=DftiComputeForward(handle, reinterpret_cast<T*>(a.begin())) ) != 0 ) 
            { 
                GADGET_ERROR_MSG( DftiErrorMessage(res) ); 
                return false; 
            }
        }
        else
        {
            if ( ( res=DftiComputeBackward(handle, reinterpret_cast<T*>(a.begin())) ) != 0 ) 
            { 
                GADGET_ERROR_MSG( DftiErrorMessage(res) ); 
                return false; 
            }
        }

        if ( ( res=DftiFreeDescriptor(&handle) ) != 0 ) 
        { 
            GADGET_ERROR_MSG( DftiErrorMessage(res) ); 
            return false; 
        }

        return true;
    }

    template<typename T> 
    bool hoNDFFT<T>::fft3_mkl(hoNDArray< ComplexType >& a, hoNDArray< ComplexType >& r, bool forward)
    {
        size_t n = a.get_number_of_elements()/(a.get_size(0)*a.get_size(1)*a.get_size(2));

        MKL_LONG dim[3];
        dim[0] = a.get_size(2);
        dim[1] = a.get_size(1);
        dim[2] = a.get_size(0);

        DFTI_DESCRIPTOR_HANDLE handle;

        if ( typeid(T) == typeid(float) )
        {
            GADGET_CHECK_RETURN_FALSE(configureFFTHandleOutOfPlace(3, dim, DFTI_SINGLE, n, handle));
        }
        else if ( typeid(T) == typeid(double) )
        {
            GADGET_CHECK_RETURN_FALSE(configureFFTHandleOutOfPlace(3, dim, DFTI_DOUBLE, n, handle));
        }
        else
        {
            GADGET_ERROR_MSG("hoNDFFT<T>::fft3_mkl(hoNDArray< ComplexType >& a, hoNDArray< ComplexType >& r), only float and double are supported ... ");
            return false;
        }

        MKL_LONG res;
        if ( forward )
        {
            if ( ( res=DftiComputeForward(handle, reinterpret_cast<T*>(a.begin()), reinterpret_cast<T*>(r.begin())) ) != 0 ) 
            { 
                GADGET_ERROR_MSG( DftiErrorMessage(res) ); 
                return false; 
            }
        }
        else
        {
            if ( ( res=DftiComputeBackward(handle, reinterpret_cast<T*>(a.begin()), reinterpret_cast<T*>(r.begin())) ) != 0 ) 
            { 
                GADGET_ERROR_MSG( DftiErrorMessage(res) ); 
                return false; 
            }
        }

        if ( ( res=DftiFreeDescriptor(&handle) ) != 0 ) 
        { 
            GADGET_ERROR_MSG( DftiErrorMessage(res) ); 
            return false; 
        }

        return true;
    }

#endif // USE_MKL

    // 
    // Instantiation
    //

    template class EXPORTCPUFFT hoNDFFT<float>;
    template class EXPORTCPUFFT hoNDFFT<double>;
}
