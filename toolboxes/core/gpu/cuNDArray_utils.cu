#include "cuNDArray_utils.h"
#include "vector_td_utilities.h"
#include "cudaDeviceManager.h"

#include <cmath>

namespace Gadgetron {

    template <class T> 
    __global__ void cuNDArray_permute_kernel(
        T* in, T* out, 
        unsigned int ndim,
        unsigned int* dims,
        unsigned int* strides_out,
        unsigned long int elements,
        int shift_mode)
    {
        unsigned long idx_in = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x+threadIdx.x;
        unsigned long idx_out = 0;

        unsigned long idx_in_tmp = idx_in;
        if (idx_in < elements) {

            unsigned int cur_index;
            for (unsigned int i = 0; i < ndim; i++) {
                unsigned long idx_in_remainder = idx_in_tmp / dims[i];
                cur_index = idx_in_tmp-(idx_in_remainder*dims[i]); //cur_index = idx_in_tmp%dims[i];
                if (shift_mode < 0) { //IFFTSHIFT
                    idx_out += ((cur_index+(dims[i]>>1))%dims[i])*strides_out[i];
                } else if (shift_mode > 0) { //FFTSHIFT
                    idx_out += ((cur_index+((dims[i]+1)>>1))%dims[i])*strides_out[i];
                } else {
                    idx_out += cur_index*strides_out[i];
                }
                idx_in_tmp = idx_in_remainder;
            }
            out[idx_in] = in[idx_out];
        }
    }

    template <class T> void cuNDArray_permute(cuNDArray<T>* in,
        cuNDArray<T>* out,
        std::vector<unsigned int> *order,
        int shift_mode)
    {

        if( out == 0x0 ){
            BOOST_THROW_EXCEPTION(cuda_error("cuNDArray_permute(internal): 0x0 output"));
        }

        cudaError_t err;

        T* in_ptr = in->get_data_ptr();
        T* out_ptr = 0;

        if (out) {
            out_ptr = out->get_data_ptr();
        } else {
            if (cudaMalloc((void**) &out_ptr, in->get_number_of_elements()*sizeof(T)) != cudaSuccess) {
                BOOST_THROW_EXCEPTION(cuda_error("cuNDArray_permute : Error allocating CUDA memory"));
            }
        }

        unsigned int* dims        = new unsigned int[in->get_number_of_dimensions()];
        unsigned int* strides_out = new unsigned int[in->get_number_of_dimensions()];
        if (!dims || !strides_out) {
            BOOST_THROW_EXCEPTION(cuda_error("cuNDArray_permute: failed to allocate temporary storage for arrays"));
        }

        for (unsigned int i = 0; i < in->get_number_of_dimensions(); i++) {
            dims[i] = (*in->get_dimensions())[(*order)[i]];
            strides_out[i] = 1;    
            for (unsigned int j = 0; j < (*order)[i]; j++) {
                strides_out[i] *= (*in->get_dimensions())[j];
            }
        }

        unsigned int* dims_dev        = 0;
        unsigned int* strides_out_dev = 0;

        if (cudaMalloc((void**) &dims_dev, in->get_number_of_dimensions()*sizeof(unsigned int)) != cudaSuccess) {
            BOOST_THROW_EXCEPTION(cuda_error("cuNDArray_permute : Error allocating CUDA dims memory"));      
        }

        if (cudaMalloc((void**) &strides_out_dev, in->get_number_of_dimensions()*sizeof(unsigned int)) != cudaSuccess) {
            BOOST_THROW_EXCEPTION(cuda_error("cuNDArray_permute : Error allocating CUDA strides_out memory"));
        }

        if (cudaMemcpy(dims_dev, dims, in->get_number_of_dimensions()*sizeof(unsigned int), cudaMemcpyHostToDevice) != cudaSuccess) {
            err = cudaGetLastError();
            std::stringstream ss;
            ss << "cuNDArray_permute : Error uploading dimensions to device, " << cudaGetErrorString(err);
            BOOST_THROW_EXCEPTION(cuda_error(ss.str()));
        }

        if (cudaMemcpy(strides_out_dev, strides_out, in->get_number_of_dimensions()*sizeof(unsigned int), cudaMemcpyHostToDevice) != cudaSuccess) {
            BOOST_THROW_EXCEPTION(cuda_error("cuNDArray_permute : Error uploading strides to device"));      
        }

        dim3 blockDim(512,1,1);
        dim3 gridDim;
        if( in->get_number_of_dimensions() > 2 ){
            gridDim = dim3((unsigned int) std::ceil((double)in->get_size(0)*in->get_size(1)/blockDim.x), 1, 1 );
            for( unsigned int d=2; d<in->get_number_of_dimensions(); d++ )
                gridDim.y *= in->get_size(d);
        }
        else
            gridDim = dim3((unsigned int) std::ceil((double)in->get_number_of_elements()/blockDim.x), 1, 1 );

        cuNDArray_permute_kernel<<< gridDim, blockDim >>>( in_ptr, out_ptr, in->get_number_of_dimensions(), 
            dims_dev, strides_out_dev, in->get_number_of_elements(), shift_mode);

        err = cudaGetLastError();
        if( err != cudaSuccess ){
            std::stringstream ss;
            ss <<"cuNDArray_permute : Error during kernel call: " << cudaGetErrorString(err);
            BOOST_THROW_EXCEPTION(cuda_error(ss.str()));      
        }

        if (cudaFree(dims_dev) != cudaSuccess) {
            err = cudaGetLastError();
            std::stringstream ss;
            ss << "cuNDArray_permute: failed to delete device memory (dims_dev) " << cudaGetErrorString(err);
            BOOST_THROW_EXCEPTION(cuda_error(ss.str()));
        }

        if (cudaFree(strides_out_dev) != cudaSuccess) {
            err = cudaGetLastError();
            std::stringstream ss;
            ss << "cuNDArray_permute: failed to delete device memory (strides_out_dev) "<< cudaGetErrorString(err);
            BOOST_THROW_EXCEPTION(cuda_error(ss.str()));
        }    
        delete [] dims;
        delete [] strides_out;    
    }  

    template <class T> boost::shared_ptr< cuNDArray<T> >
    permute( cuNDArray<T> *in, std::vector<unsigned int> *dim_order, int shift_mode )
    {
        if( in == 0x0 || dim_order == 0x0 ) {
            BOOST_THROW_EXCEPTION(runtime_error("permute(): invalid pointer provided"));
        }    

        std::vector<unsigned int> dims;
        for (unsigned int i = 0; i < dim_order->size(); i++)
            dims.push_back(in->get_dimensions()->at(dim_order->at(i)));
        boost::shared_ptr< cuNDArray<T> > out( new cuNDArray<T>() );    
        out->create(&dims);
        permute( in, out.get(), dim_order, shift_mode );
        return out;
    }

    template <class T> void
        permute( cuNDArray<T> *in, cuNDArray<T> *out, std::vector<unsigned int> *dim_order, int shift_mode )
    {
        if( in == 0x0 || out == 0x0 || dim_order == 0x0 ) {
            BOOST_THROW_EXCEPTION(runtime_error("permute(): invalid pointer provided"));
        }    

        //Check ordering array
        if (dim_order->size() > in->get_number_of_dimensions()) {
            BOOST_THROW_EXCEPTION(runtime_error("permute(): invalid length of dimension ordering array"));
        }

        std::vector<unsigned int> dim_count(in->get_number_of_dimensions(),0);
        for (unsigned int i = 0; i < dim_order->size(); i++) {
            if ((*dim_order)[i] >= in->get_number_of_dimensions()) {
                BOOST_THROW_EXCEPTION(runtime_error("permute(): invalid dimension order array"));
            }
            dim_count[(*dim_order)[i]]++;
        }

        //Create an internal array to store the dimensions
        std::vector<unsigned int> dim_order_int;

        //Check that there are no duplicate dimensions
        for (unsigned int i = 0; i < dim_order->size(); i++) {
            if (dim_count[(*dim_order)[i]] != 1) {
                BOOST_THROW_EXCEPTION(runtime_error("permute(): invalid dimension order array (duplicates)"));
            }
            dim_order_int.push_back((*dim_order)[i]);
        }

        for (unsigned int i = 0; i < dim_order_int.size(); i++) {
            if ((*in->get_dimensions())[dim_order_int[i]] != out->get_size(i)) {
                BOOST_THROW_EXCEPTION(runtime_error("permute(): dimensions of output array do not match the input array"));
            }
        }

        //Pad dimension order array with dimension not mentioned in order array
        if (dim_order_int.size() < in->get_number_of_dimensions()) {
            for (unsigned int i = 0; i < dim_count.size(); i++) {
                if (dim_count[i] == 0) {
                    dim_order_int.push_back(i);
                }
            }
        }    
        cuNDArray_permute(in, out, &dim_order_int, shift_mode);
    }

    template<class T> boost::shared_ptr< cuNDArray<T> >
    shift_dim( cuNDArray<T> *in, int shift )
    {
        if( in == 0x0 ) {
            BOOST_THROW_EXCEPTION(runtime_error("shift_dim(): invalid input pointer provided"));
        }    

        std::vector<unsigned int> order;
        for (int i = 0; i < in->get_number_of_dimensions(); i++) {
            order.push_back(static_cast<unsigned int>((i+shift)%in->get_number_of_dimensions()));
        }
        return permute(in,&order);
    }

    template<class T> 
    void shift_dim( cuNDArray<T> *in, cuNDArray<T> *out, int shift )
    {
        if( in == 0x0 || out == 0x0 ) {
            BOOST_THROW_EXCEPTION(runtime_error("shift_dim(): invalid pointer provided"));
        }    

        std::vector<unsigned int> order;
        for (int i = 0; i < in->get_number_of_dimensions(); i++) {
            order.push_back(static_cast<unsigned int>((i+shift)%in->get_number_of_dimensions()));
        }
        permute(in,out,&order);
    }

    template<class T> static void find_stride( cuNDArray<T> *in, unsigned int dim, unsigned int *stride, std::vector<unsigned int> *dims )
    {
        *stride = 1;
        for( unsigned int i=0; i<in->get_number_of_dimensions(); i++ ){
            if( i != dim )
                dims->push_back(in->get_size(i));
            if( i < dim )
                *stride *= in->get_size(i);
        }
    }

    void setup_grid( unsigned int number_of_elements, dim3 *blockDim, dim3* gridDim, unsigned int num_batches = 1 )
    {    
        int cur_device = cudaDeviceManager::Instance()->getCurrentDevice();
        int maxGridDim = cudaDeviceManager::Instance()->max_griddim(cur_device);

        // For small arrays we keep the block dimension fairly small
        *blockDim = dim3(256);
        *gridDim = dim3((number_of_elements+blockDim->x-1)/blockDim->x, num_batches);

        // Extend block/grid dimensions for large arrays
        if( gridDim->x > maxGridDim){
            blockDim->x = maxGridDim;
            gridDim->x = (number_of_elements+blockDim->x-1)/blockDim->x;
        }

        if( gridDim->x > maxGridDim ){
            gridDim->x = ((unsigned int)std::sqrt((float)number_of_elements)+blockDim->x-1)/blockDim->x;
            gridDim->y *= ((number_of_elements+blockDim->x*gridDim->x-1)/(blockDim->x*gridDim->x));
        }

        if( gridDim->x >maxGridDim || gridDim->y >maxGridDim){      
            BOOST_THROW_EXCEPTION(cuda_error("Grid dimension larger than supported by device"));
        }
    }

    // Expand
    //
    template<class T> 
    __global__ void expand_kernel( 
        T *in, T *out, 
        unsigned int number_of_elements_in, unsigned int number_of_elements_out, unsigned int new_dim_size )
    {
        const unsigned int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x+threadIdx.x;    
        if( idx < number_of_elements_out ){
            out[idx] = in[idx%number_of_elements_in];
        }
    }

    // Expand
    //
    template<class T> boost::shared_ptr< cuNDArray<T> > 
    expand( cuNDArray<T> *in, unsigned int new_dim_size )
    {
        unsigned int number_of_elements_out = in->get_number_of_elements()*new_dim_size;

        // Setup block/grid dimensions
        dim3 blockDim; dim3 gridDim;
        setup_grid( number_of_elements_out, &blockDim, &gridDim );

        // Find element stride
        std::vector<unsigned int> dims = *in->get_dimensions();
        dims.push_back(new_dim_size);

        // Invoke kernel
        boost::shared_ptr< cuNDArray<T> > out( new cuNDArray<T>());
        out->create(&dims);

        expand_kernel<T><<< gridDim, blockDim >>>( in->get_data_ptr(), out->get_data_ptr(), 
            in->get_number_of_elements(), number_of_elements_out, new_dim_size );

        CHECK_FOR_CUDA_ERROR();    
        return out;
    }

    // Sum
    //
    template<class T> 
    __global__ void sum_kernel( 
        T *in, T *out, 
        unsigned int stride, unsigned int number_of_batches, unsigned int number_of_elements )
    {
        const unsigned int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x+threadIdx.x;

        if( idx < number_of_elements ){

            unsigned int in_idx = (idx/stride)*stride*number_of_batches+(idx%stride);

            T val = in[in_idx];

            for( unsigned int i=1; i<number_of_batches; i++ ) 
                val += in[i*stride+in_idx];

            out[idx] = val; 
        }
    }

    // Sum
    //
    template<class T>  boost::shared_ptr< cuNDArray<T> > sum( cuNDArray<T> *in, unsigned int dim )
    {
        // Some validity checks
        if( !(in->get_number_of_dimensions()>1) ){
            BOOST_THROW_EXCEPTION(runtime_error("sum: underdimensioned."));
        }

        if( dim > in->get_number_of_dimensions()-1 ){
            BOOST_THROW_EXCEPTION(runtime_error( "sum: dimension out of range."));
        }

        unsigned int number_of_batches = in->get_size(dim);
        unsigned int number_of_elements = in->get_number_of_elements()/number_of_batches;

        // Setup block/grid dimensions
        dim3 blockDim; dim3 gridDim;
        setup_grid( number_of_elements, &blockDim, &gridDim );

        // Find element stride
        unsigned int stride; std::vector<unsigned int> dims;
        find_stride<T>( in, dim, &stride, &dims );

        // Invoke kernel
        boost::shared_ptr< cuNDArray<T> > out(new cuNDArray<T>());
        out->create(&dims);

        sum_kernel<T><<< gridDim, blockDim >>>( in->get_data_ptr(), out->get_data_ptr(), stride, number_of_batches, number_of_elements );

        CHECK_FOR_CUDA_ERROR();
        return out;
    }

    // Crop
    template<class T, unsigned int D> __global__ void crop_kernel( 
        vector_td<unsigned int,D> offset, vector_td<unsigned int,D> matrix_size_in, vector_td<unsigned int,D> matrix_size_out,
        T *in, T *out, unsigned int num_batches, unsigned int num_elements )
    {
        typedef vector_td<unsigned int,D> uintd;
        const unsigned int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x+threadIdx.x;
        const unsigned int frame_offset = idx/num_elements;

        if( idx < num_elements*num_batches ){
            const uintd co = idx_to_co<D>( idx-frame_offset*num_elements, matrix_size_out );
            const uintd co_os = offset + co;
            const unsigned int in_idx = co_to_idx<D>(co_os, matrix_size_in)+frame_offset*prod(matrix_size_in);
            out[idx] = in[in_idx];
        }
    }

    // Crop
    template<class T, unsigned int D>
    void crop( typename uintd<D>::Type offset, cuNDArray<T> *in, cuNDArray<T> *out )
    {
        if( in == 0x0 || out == 0x0 ){
            BOOST_THROW_EXCEPTION(runtime_error("crop: 0x0 ndarray provided"));
        }

        if( in->get_number_of_dimensions() != out->get_number_of_dimensions() ){
            BOOST_THROW_EXCEPTION(runtime_error("crop: image dimensions mismatch"));
        }

        if( in->get_number_of_dimensions() < D ){
            std::stringstream ss;
            ss << "crop: number of image dimensions should be at least " << D;
            BOOST_THROW_EXCEPTION(runtime_error(ss.str()));
        }

        typename uintd<D>::Type matrix_size_in = from_std_vector<unsigned int,D>( *in->get_dimensions() );
        typename uintd<D>::Type matrix_size_out = from_std_vector<unsigned int,D>( *out->get_dimensions() );

        unsigned int number_of_batches = 1;
        for( unsigned int d=D; d<in->get_number_of_dimensions(); d++ ){
            number_of_batches *= in->get_size(d);
        }

        if( weak_greater(offset+matrix_size_out, matrix_size_in) ){
            BOOST_THROW_EXCEPTION(runtime_error( "crop: cropping size mismatch"));
        }

        // Setup block/grid dimensions
        dim3 blockDim; dim3 gridDim;
        setup_grid( prod(matrix_size_out), &blockDim, &gridDim, number_of_batches );

        // Invoke kernel
        crop_kernel<T,D><<< gridDim, blockDim >>>
            ( offset, matrix_size_in, matrix_size_out, in->get_data_ptr(), out->get_data_ptr(), number_of_batches, prod(matrix_size_out) );

        CHECK_FOR_CUDA_ERROR();
    }

    template<class T, unsigned int D> boost::shared_ptr< cuNDArray<T> > 
    crop( typename uintd<D>::Type offset, typename uintd<D>::Type size, cuNDArray<T> *in )
    {
        if( in == 0x0 ){
            BOOST_THROW_EXCEPTION(runtime_error("crop: 0x0 array provided"));
        }
        std::vector<unsigned int> dims = to_std_vector(size);
        for( unsigned int d=D; d<in->get_number_of_dimensions(); d++ ){
            dims.push_back(in->get_size(d));
        }
        boost::shared_ptr< cuNDArray<T> > result( new cuNDArray<T>(&dims) );
        crop<T,D>(offset, in, result.get());
        return result;
    }  

    // Expand and zero fill
    template<class T, unsigned int D> 
    __global__ void pad_kernel( 
        vector_td<unsigned int,D> matrix_size_in, vector_td<unsigned int,D> matrix_size_out,
        T *in, T *out, unsigned int number_of_batches, unsigned int num_elements, T val )
    {
        typedef vector_td<unsigned int,D> uintd;
        const unsigned int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x+threadIdx.x;
        const unsigned int frame_offset = idx/num_elements;

        if( idx < num_elements*number_of_batches ){

            const uintd co_out = idx_to_co<D>( idx-frame_offset*num_elements, matrix_size_out );
            const uintd offset = (matrix_size_out-matrix_size_in)>>1;
            T _out;
            bool inside = (co_out>=offset) && (co_out<(matrix_size_in+offset));

            if( inside )
                _out = in[co_to_idx<D>(co_out-offset, matrix_size_in)+frame_offset*prod(matrix_size_in)];
            else{      
                _out = val;
            }

            out[idx] = _out;
        }
    }

    template<class T, unsigned int D> 
    void pad( cuNDArray<T> *in, cuNDArray<T> *out, T val )
    { 
        if( in == 0x0 || out == 0x0 ){
            BOOST_THROW_EXCEPTION(runtime_error("pad: 0x0 ndarray provided"));
        }

        if( in->get_number_of_dimensions() != out->get_number_of_dimensions() ){
            BOOST_THROW_EXCEPTION(runtime_error("pad: image dimensions mismatch"));
        }

        if( in->get_number_of_dimensions() < D ){
            std::stringstream ss;
            ss << "pad: number of image dimensions should be at least " << D;
            BOOST_THROW_EXCEPTION(runtime_error(ss.str()));
        }

        typename uintd<D>::Type matrix_size_in = from_std_vector<unsigned int,D>( *in->get_dimensions() );
        typename uintd<D>::Type matrix_size_out = from_std_vector<unsigned int,D>( *out->get_dimensions() );

        unsigned int number_of_batches = 1;
        for( unsigned int d=D; d<in->get_number_of_dimensions(); d++ ){
            number_of_batches *= in->get_size(d);
        }

        if( weak_greater(matrix_size_in,matrix_size_out) ){
            std::runtime_error("pad: size mismatch, cannot expand");
        }

        // Setup block/grid dimensions
        dim3 blockDim; dim3 gridDim;
        setup_grid( prod(matrix_size_out), &blockDim, &gridDim, number_of_batches );

        // Invoke kernel
        pad_kernel<T,D><<< gridDim, blockDim >>> 
            ( matrix_size_in, matrix_size_out, in->get_data_ptr(), out->get_data_ptr(), number_of_batches, prod(matrix_size_out), val );

        CHECK_FOR_CUDA_ERROR();
    }

    template<class T, unsigned int D> boost::shared_ptr< cuNDArray<T> >
    pad( typename uintd<D>::Type size, cuNDArray<T> *in, T val )
    {
        if( in == 0x0 ){
            BOOST_THROW_EXCEPTION(runtime_error("pad: 0x0 array provided"));
        }
        std::vector<unsigned int> dims = to_std_vector(size);
        for( unsigned int d=D; d<in->get_number_of_dimensions(); d++ ){
            dims.push_back(in->get_size(d));
        }
        boost::shared_ptr< cuNDArray<T> > result( new cuNDArray<T>(&dims) );
        pad<T,D>(in, result.get(), val);
        return result;
    }

    template<class T, unsigned int D> 
    __global__ void fill_border_kernel( 
        vector_td<unsigned int,D> matrix_size_in, vector_td<unsigned int,D> matrix_size_out,
        T *image, unsigned int number_of_batches, unsigned int number_of_elements, T val )
    {
        const unsigned int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x+threadIdx.x;

        if( idx < number_of_elements ){
            const vector_td<unsigned int,D> co_out = idx_to_co<D>( idx, matrix_size_out );
            const vector_td<unsigned int,D> offset = (matrix_size_out-matrix_size_in)>>1;
            if( weak_less( co_out, offset ) || weak_greater_equal( co_out, matrix_size_in+offset ) ){
                for( unsigned int batch=0; batch<number_of_batches; batch++ ){
                    image[idx+batch*number_of_elements] = val;
                }
            }
            else
                ; // do nothing
        }
    }

    // Zero fill border (rectangular)
    template<class T, unsigned int D> 
    void fill_border( typename uintd<D>::Type matrix_size_in, cuNDArray<T> *in_out, T val )
    { 
        typename uintd<D>::Type matrix_size_out = from_std_vector<unsigned int,D>( *in_out->get_dimensions() );

        if( weak_greater(matrix_size_in, matrix_size_out) ){
            BOOST_THROW_EXCEPTION(runtime_error("fill_border: size mismatch, cannot zero fill"));
        }

        unsigned int number_of_batches = 1;
        for( unsigned int d=D; d<in_out->get_number_of_dimensions(); d++ ){
            number_of_batches *= in_out->get_size(d);
        }

        // Setup block/grid dimensions
        dim3 blockDim; dim3 gridDim;
        setup_grid( prod(matrix_size_out), &blockDim, &gridDim );

        // Invoke kernel
        fill_border_kernel<T,D><<< gridDim, blockDim >>>
            ( matrix_size_in, matrix_size_out, in_out->get_data_ptr(), number_of_batches, prod(matrix_size_out), val );

        CHECK_FOR_CUDA_ERROR();
    }

    //
    // Instantiation
    //

    template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> > permute( cuNDArray<float>*, std::vector<unsigned int>*, int );
    template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> > permute( cuNDArray<double>*, std::vector<unsigned int>*, int );
    template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float_complext> > permute( cuNDArray<float_complext>*, std::vector<unsigned int>*, int );
    template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double_complext> > permute( cuNDArray<double_complext>*, std::vector<unsigned int>*, int );  

    template EXPORTGPUCORE void permute( cuNDArray<float>*, cuNDArray<float>*, std::vector<unsigned int>*, int);
    template EXPORTGPUCORE void permute( cuNDArray<double>*, cuNDArray<double>*, std::vector<unsigned int>*, int);
    template EXPORTGPUCORE void permute( cuNDArray<float_complext>*, cuNDArray<float_complext>*, std::vector<unsigned int>*, int);
    template EXPORTGPUCORE void permute( cuNDArray<double_complext>*, cuNDArray<double_complext>*, std::vector<unsigned int>*, int);

    template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> > shift_dim( cuNDArray<float>*, int );
    template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> > shift_dim( cuNDArray<double>*, int );
    template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float_complext> > shift_dim( cuNDArray<float_complext>*, int );
    template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double_complext> > shift_dim( cuNDArray<double_complext>*, int );

    template EXPORTGPUCORE void shift_dim( cuNDArray<float>*, cuNDArray<float>*, int shift );
    template EXPORTGPUCORE void shift_dim( cuNDArray<double>*, cuNDArray<double>*, int shift );
    template EXPORTGPUCORE void shift_dim( cuNDArray<float_complext>*, cuNDArray<float_complext>*, int shift );
    template EXPORTGPUCORE void shift_dim( cuNDArray<double_complext>*, cuNDArray<double_complext>*, int shift );

    template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> > expand<float>( cuNDArray<float>*, unsigned int);
    template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> > expand<double>( cuNDArray<double>*, unsigned int);
    template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float_complext> > expand<float_complext>( cuNDArray<float_complext>*, unsigned int);
    template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double_complext> > expand<double_complext>( cuNDArray<double_complext>*, unsigned int);

    template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> > sum<float>( cuNDArray<float>*, unsigned int);
    template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> > sum<double>( cuNDArray<double>*, unsigned int);
    template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float_complext> > sum<float_complext>( cuNDArray<float_complext>*, unsigned int);
    template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double_complext> > sum<double_complext>( cuNDArray<double_complext>*, unsigned int);  

    template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> > crop<float,1>( typename uintd<1>::Type, typename uintd<1>::Type, cuNDArray<float>*);
    template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> > crop<float,2>( typename uintd<2>::Type, typename uintd<2>::Type, cuNDArray<float>*);
    template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> > crop<float,3>( typename uintd<3>::Type, typename uintd<3>::Type, cuNDArray<float>*);
    template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> > crop<float,4>( typename uintd<4>::Type, typename uintd<4>::Type, cuNDArray<float>*);

    template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float_complext> > crop<float_complext,1>( typename uintd<1>::Type, typename uintd<1>::Type, cuNDArray<float_complext>*);
    template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float_complext> > crop<float_complext,2>( typename uintd<2>::Type, typename uintd<2>::Type, cuNDArray<float_complext>*);
    template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float_complext> > crop<float_complext,3>( typename uintd<3>::Type, typename uintd<3>::Type, cuNDArray<float_complext>*);
    template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float_complext> > crop<float_complext,4>( typename uintd<4>::Type, typename uintd<4>::Type, cuNDArray<float_complext>*);

    template EXPORTGPUCORE void crop<float,1>( uintd1, cuNDArray<float>*, cuNDArray<float>*);
    template EXPORTGPUCORE void crop<float,2>( uintd2, cuNDArray<float>*, cuNDArray<float>*);
    template EXPORTGPUCORE void crop<float,3>( uintd3, cuNDArray<float>*, cuNDArray<float>*);
    template EXPORTGPUCORE void crop<float,4>( uintd4, cuNDArray<float>*, cuNDArray<float>*);

    template EXPORTGPUCORE void crop<complext<float>,1>( uintd1, cuNDArray<complext<float> >*, cuNDArray< complext<float> >*);
    template EXPORTGPUCORE void crop<complext<float>,2>( uintd2, cuNDArray<complext<float> >*, cuNDArray< complext<float> >*);
    template EXPORTGPUCORE void crop<complext<float>,3>( uintd3, cuNDArray<complext<float> >*, cuNDArray< complext<float> >*);
    template EXPORTGPUCORE void crop<complext<float>,4>( uintd4, cuNDArray<complext<float> >*, cuNDArray< complext<float> >*);

    template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> > pad<float,1>( typename uintd<1>::Type, cuNDArray<float>*, float );
    template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> > pad<float,2>( typename uintd<2>::Type, cuNDArray<float>*, float );
    template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> > pad<float,3>( typename uintd<3>::Type, cuNDArray<float>*, float );
    template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> > pad<float,4>( typename uintd<4>::Type, cuNDArray<float>*, float );

    template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float_complext> > pad<float_complext,1>( typename uintd<1>::Type, cuNDArray<float_complext>*, float_complext );
    template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float_complext> > pad<float_complext,2>( typename uintd<2>::Type, cuNDArray<float_complext>*, float_complext );
    template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float_complext> > pad<float_complext,3>( typename uintd<3>::Type, cuNDArray<float_complext>*, float_complext );
    template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float_complext> > pad<float_complext,4>( typename uintd<4>::Type, cuNDArray<float_complext>*, float_complext );

    template EXPORTGPUCORE void pad<float,1>( cuNDArray<float>*, cuNDArray<float>*, float);
    template EXPORTGPUCORE void pad<float,2>( cuNDArray<float>*, cuNDArray<float>*, float);
    template EXPORTGPUCORE void pad<float,3>( cuNDArray<float>*, cuNDArray<float>*, float);
    template EXPORTGPUCORE void pad<float,4>( cuNDArray<float>*, cuNDArray<float>*, float);

    template EXPORTGPUCORE void pad<float_complext,1>( cuNDArray<float_complext>*, cuNDArray<float_complext>*, float_complext);
    template EXPORTGPUCORE void pad<float_complext,2>( cuNDArray<float_complext>*, cuNDArray<float_complext>*, float_complext);  
    template EXPORTGPUCORE void pad<float_complext,3>( cuNDArray<float_complext>*, cuNDArray<float_complext>*, float_complext);
    template EXPORTGPUCORE void pad<float_complext,4>( cuNDArray<float_complext>*, cuNDArray<float_complext>*, float_complext);

    template EXPORTGPUCORE void fill_border<float,1>(uintd1, cuNDArray<float>*,float);
    template EXPORTGPUCORE void fill_border<float,2>(uintd2, cuNDArray<float>*,float);
    template EXPORTGPUCORE void fill_border<float,3>(uintd3, cuNDArray<float>*,float);
    template EXPORTGPUCORE void fill_border<float,4>(uintd4, cuNDArray<float>*,float);

    template EXPORTGPUCORE void fill_border<float_complext,1>(uintd1, cuNDArray<float_complext>*,float_complext);
    template EXPORTGPUCORE void fill_border<float_complext,2>(uintd2, cuNDArray<float_complext>*,float_complext);
    template EXPORTGPUCORE void fill_border<float_complext,3>(uintd3, cuNDArray<float_complext>*,float_complext);
    template EXPORTGPUCORE void fill_border<float_complext,4>(uintd4, cuNDArray<float_complext>*,float_complext);

    template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> > crop<double,1>( typename uintd<1>::Type, typename uintd<1>::Type, cuNDArray<double>*);
    template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> > crop<double,2>( typename uintd<2>::Type, typename uintd<2>::Type, cuNDArray<double>*);
    template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> > crop<double,3>( typename uintd<3>::Type, typename uintd<3>::Type, cuNDArray<double>*);
    template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> > crop<double,4>( typename uintd<4>::Type, typename uintd<4>::Type, cuNDArray<double>*);

    template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double_complext> > crop<double_complext,1>( typename uintd<1>::Type, typename uintd<1>::Type, cuNDArray<double_complext>*);
    template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double_complext> > crop<double_complext,2>( typename uintd<2>::Type, typename uintd<2>::Type, cuNDArray<double_complext>*);
    template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double_complext> > crop<double_complext,3>( typename uintd<3>::Type, typename uintd<3>::Type, cuNDArray<double_complext>*);
    template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double_complext> > crop<double_complext,4>( typename uintd<4>::Type, typename uintd<4>::Type, cuNDArray<double_complext>*);

    template EXPORTGPUCORE void crop<double,1>( uintd1, cuNDArray<double>*, cuNDArray<double>*);
    template EXPORTGPUCORE void crop<double,2>( uintd2, cuNDArray<double>*, cuNDArray<double>*);
    template EXPORTGPUCORE void crop<double,3>( uintd3, cuNDArray<double>*, cuNDArray<double>*);
    template EXPORTGPUCORE void crop<double,4>( uintd4, cuNDArray<double>*, cuNDArray<double>*);

    template EXPORTGPUCORE void crop<complext<double>,1>( uintd1, cuNDArray<complext<double> >*, cuNDArray< complext<double> >*);
    template EXPORTGPUCORE void crop<complext<double>,2>( uintd2, cuNDArray<complext<double> >*, cuNDArray< complext<double> >*);
    template EXPORTGPUCORE void crop<complext<double>,3>( uintd3, cuNDArray<complext<double> >*, cuNDArray< complext<double> >*);
    template EXPORTGPUCORE void crop<complext<double>,4>( uintd4, cuNDArray<complext<double> >*, cuNDArray< complext<double> >*);

    template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> > pad<double,1>( typename uintd<1>::Type, cuNDArray<double>*, double );
    template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> > pad<double,2>( typename uintd<2>::Type, cuNDArray<double>*, double );
    template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> > pad<double,3>( typename uintd<3>::Type, cuNDArray<double>*, double );
    template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> > pad<double,4>( typename uintd<4>::Type, cuNDArray<double>*, double );

    template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double_complext> > pad<double_complext,1>( typename uintd<1>::Type, cuNDArray<double_complext>*, double_complext );
    template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double_complext> > pad<double_complext,2>( typename uintd<2>::Type, cuNDArray<double_complext>*, double_complext );
    template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double_complext> > pad<double_complext,3>( typename uintd<3>::Type, cuNDArray<double_complext>*, double_complext );
    template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double_complext> > pad<double_complext,4>( typename uintd<4>::Type, cuNDArray<double_complext>*, double_complext );

    template EXPORTGPUCORE void pad<double,1>( cuNDArray<double>*, cuNDArray<double>*, double);
    template EXPORTGPUCORE void pad<double,2>( cuNDArray<double>*, cuNDArray<double>*, double);
    template EXPORTGPUCORE void pad<double,3>( cuNDArray<double>*, cuNDArray<double>*, double);
    template EXPORTGPUCORE void pad<double,4>( cuNDArray<double>*, cuNDArray<double>*, double);

    template EXPORTGPUCORE void pad<double_complext,1>( cuNDArray<double_complext>*, cuNDArray<double_complext>*, double_complext);
    template EXPORTGPUCORE void pad<double_complext,2>( cuNDArray<double_complext>*, cuNDArray<double_complext>*, double_complext);  
    template EXPORTGPUCORE void pad<double_complext,3>( cuNDArray<double_complext>*, cuNDArray<double_complext>*, double_complext);
    template EXPORTGPUCORE void pad<double_complext,4>( cuNDArray<double_complext>*, cuNDArray<double_complext>*, double_complext);

    template EXPORTGPUCORE void fill_border<double,1>(uintd1, cuNDArray<double>*,double);
    template EXPORTGPUCORE void fill_border<double,2>(uintd2, cuNDArray<double>*,double);
    template EXPORTGPUCORE void fill_border<double,3>(uintd3, cuNDArray<double>*,double);
    template EXPORTGPUCORE void fill_border<double,4>(uintd4, cuNDArray<double>*,double);

    template EXPORTGPUCORE void fill_border<double_complext,1>(uintd1, cuNDArray<double_complext>*,double_complext);
    template EXPORTGPUCORE void fill_border<double_complext,2>(uintd2, cuNDArray<double_complext>*,double_complext);
    template EXPORTGPUCORE void fill_border<double_complext,3>(uintd3, cuNDArray<double_complext>*,double_complext);
    template EXPORTGPUCORE void fill_border<double_complext,4>(uintd4, cuNDArray<double_complext>*,double_complext);
}
