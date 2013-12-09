#pragma once

#include "hoNDArray.h"
#include "vector_td_utilities.h"

namespace Gadgetron {

    class ArrayIterator
    {
    public:

        ArrayIterator(std::vector<unsigned long long> *dimensions, std::vector<unsigned long long> *order)
        {
            dimensions_  = boost::shared_ptr< std::vector<unsigned long long> >      (new std::vector<unsigned long long>);
            order_       = boost::shared_ptr< std::vector<unsigned long long> >      (new std::vector<unsigned long long>);
            current_     = boost::shared_ptr< std::vector<unsigned long long> >      (new std::vector<unsigned long long>);
            block_sizes_ = boost::shared_ptr< std::vector<unsigned long int> > (new std::vector<unsigned long int>);

            block_sizes_->push_back(1);
            for (unsigned long long i = 0; i < order->size(); i++) {
                dimensions_->push_back((*dimensions)[i]);
                order_->push_back((*order)[i]);
                current_->push_back(0);
                if (i > 0) {
                    block_sizes_->push_back((*block_sizes_)[i-1]*(*dimensions_)[i-1]);
                }
            }
            current_idx_ = 0;
        }

        inline unsigned long int advance()
        {
            unsigned long long order_index = 0;
            (*current_)[(*order_)[order_index]]++;
            while ((*current_)[(*order_)[order_index]] >= (*dimensions_)[(*order_)[order_index]]) {
                (*current_)[(*order_)[order_index]] = 0;
                order_index = (order_index+1)%dimensions_->size();
                (*current_)[(*order_)[order_index]]++;
            }

            current_idx_ = 0;
            for (unsigned long long i = 0; i < dimensions_->size(); i++) {
                current_idx_ += (*current_)[i]*(*block_sizes_)[i];
            }	
            return current_idx_;
        }

        inline unsigned long int get_current_idx() {
            return current_idx_;
        }

        boost::shared_ptr< std::vector<unsigned long long> > get_current_sub() {
            return current_;
        }

    protected:
        boost::shared_ptr< std::vector<unsigned long long> > dimensions_;
        boost::shared_ptr< std::vector<unsigned long long> > order_;
        boost::shared_ptr< std::vector<unsigned long long> > current_;
        boost::shared_ptr< std::vector<unsigned long int> > block_sizes_;
        unsigned long int current_idx_;
    };

    template<class T> boost::shared_ptr< hoNDArray<T> > shift_dim( hoNDArray<T> *in, int shift )  
    {
        if( in == 0x0 ) {
            throw std::runtime_error("shift_dim(): invalid input pointer provided");;
        }    
        std::vector<unsigned long long> order;
        for (unsigned long long i = 0; i < in->get_number_of_dimensions(); i++) {
            order.push_back(static_cast<unsigned long long>((i+shift)%in->get_number_of_dimensions()));
        }
        return permute(in,&order);
    }

    template<class T> void shift_dim( hoNDArray<T> *in, hoNDArray<T> *out, int shift )
    {
        if( in == 0x0 || out == 0x0 ) {
            throw std::runtime_error("shift_dim(): invalid pointer provided");;
        }    
        std::vector<unsigned long long> order;
        for (unsigned long long i = 0; i < in->get_number_of_dimensions(); i++) {
            order.push_back(static_cast<unsigned long long>((i+shift)%in->get_number_of_dimensions()));
        }
        permute(in,out,&order);
    }

    template<class T> boost::shared_ptr< hoNDArray<T> > 
    permute( hoNDArray<T> *in, std::vector<unsigned long long> *dim_order, int shift_mode = 0) 
    {
        if( in == 0x0 || dim_order == 0x0 ) {
            throw std::runtime_error("permute(): invalid pointer provided");;
        }    

        std::vector<unsigned long long> dims;
        for (unsigned long long i = 0; i < dim_order->size(); i++)
            dims.push_back(in->get_dimensions()->at(dim_order->at(i)));
        boost::shared_ptr< hoNDArray<T> > out( new hoNDArray<T>() );    
        out->create(&dims);
        permute( in, out.get(), dim_order, shift_mode );
        return out;
    }

    template<class T> void 
        permute( hoNDArray<T> *in, hoNDArray<T> *out, std::vector<unsigned long long> *dim_order, int shift_mode = 0) 
    {
        if( in == 0x0 || out == 0x0 || dim_order == 0x0 ) {
            throw std::runtime_error("permute(): invalid pointer provided");;
        }    

        // Check ordering array
        if (dim_order->size() > in->get_number_of_dimensions()) {
            throw std::runtime_error("hoNDArray::permute - Invalid length of dimension ordering array");;
        }

        std::vector<unsigned long long> dim_count(in->get_number_of_dimensions(),0);
        for (unsigned long long i = 0; i < dim_order->size(); i++) {
            if ((*dim_order)[i] >= in->get_number_of_dimensions()) {
                throw std::runtime_error("hoNDArray::permute - Invalid dimension order array");;
            }
            dim_count[(*dim_order)[i]]++;
        }

        // Create an internal array to store the dimensions
        std::vector<unsigned long long> dim_order_int;

        // Check that there are no duplicate dimensions
        for (unsigned long long i = 0; i < dim_order->size(); i++) {
            if (dim_count[(*dim_order)[i]] != 1) {
                throw std::runtime_error("hoNDArray::permute - Invalid dimension order array (duplicates)");;

            }
            dim_order_int.push_back((*dim_order)[i]);
        }

        for (unsigned long long i = 0; i < dim_order_int.size(); i++) {
            if ((*in->get_dimensions())[dim_order_int[i]] != out->get_size(i)) {
                throw std::runtime_error("permute(): dimensions of output array do not match the input array");;
            }
        }

        // Pad dimension order array with dimension not mentioned in order array
        if (dim_order_int.size() < in->get_number_of_dimensions()) {
            for (unsigned long long i = 0; i < dim_count.size(); i++) {
                if (dim_count[i] == 0) {
                    dim_order_int.push_back(i);
                }
            }
        }

        T* o = out->get_data_ptr();

        ArrayIterator it(in->get_dimensions().get(),&dim_order_int);
        for (unsigned long int i = 0; i < in->get_number_of_elements(); i++) {
            o[i] = in->get_data_ptr()[it.get_current_idx()];
            it.advance();
        }
    }
   
    // Expand array to new dimension
    template<class T> boost::shared_ptr<hoNDArray<T> > 
    expand(hoNDArray<T> *in, unsigned long long new_dim_size )
    {
      if( in == 0x0 ){
	throw std::runtime_error("expand(): illegal input pointer.");;
      }
      
      const unsigned long long number_of_elements_in = in->get_number_of_elements();    

      std::vector<unsigned long long> dims = *in->get_dimensions(); 
      dims.push_back(new_dim_size);

      boost::shared_ptr< hoNDArray<T> > out(new hoNDArray<T>(&dims));
      
#ifdef USE_OMP
#pragma omp parallel for
#endif
      for( int idx=0; idx<number_of_elements_in*new_dim_size; idx++ ){
	(*out)[idx] = in->at(idx%number_of_elements_in);
      }
      return out;
    }
  
    // Sum over dimension
    template<class T> boost::shared_ptr<hoNDArray<T> > 
    sum(hoNDArray<T> *in, unsigned long long dim )
    {
        if( in == 0x0 ){
            throw std::runtime_error("sum(): illegal input pointer.");;
        }

        if( !(in->get_number_of_dimensions()>1) ){
            throw std::runtime_error("sum(): underdimensioned.");;
        }

        if( dim > in->get_number_of_dimensions()-1 ){
            throw std::runtime_error( "sum(): dimension out of range.");;
        }

        unsigned long long number_of_batches = in->get_size(dim);
        unsigned long long number_of_elements = in->get_number_of_elements()/number_of_batches;
        std::vector<unsigned long long> dims = *in->get_dimensions(); dims.pop_back();

        boost::shared_ptr< hoNDArray<T> > out(new hoNDArray<T>());
        out->create(&dims);

#ifdef USE_OMP
#pragma omp parallel for
#endif
        for( int idx=0; idx<(int)number_of_elements; idx++ ){
            T val(0);
            for( unsigned long long j=0; j<number_of_batches; j++ ){
                unsigned long long in_idx = j*number_of_elements+idx;
                val += in->get_data_ptr()[in_idx];      
            }
            out->get_data_ptr()[idx] = val;       
        }
        return out;
    } 

    template<class T, unsigned long long D> boost::shared_ptr< hoNDArray<T> >
    crop( const vector_td< unsigned long long, D >& crop_offset, const vector_td< unsigned long long, D >& crop_size, hoNDArray<T> *in )
    {
        if( in == 0x0 ){
            throw std::runtime_error("crop: 0x0 array provided");;
        }

        if( in->get_number_of_dimensions() < D ){
            std::stringstream ss;
            ss << "crop: number of image dimensions should be at least " << D;
            throw std::runtime_error(ss.str());;
        }

        std::vector<unsigned long long> dims = to_std_vector(crop_size);
        for( unsigned long long d=D; d<in->get_number_of_dimensions(); d++ ){
            dims.push_back(in->get_size(d));
        }
        boost::shared_ptr< hoNDArray<T> > out( new hoNDArray<T>(&dims) );

        typename uintd<D>::Type matrix_size_in = from_std_vector<unsigned long long,D>( *in->get_dimensions() );
        typename uintd<D>::Type matrix_size_out = from_std_vector<unsigned long long,D>( *out->get_dimensions() );

        unsigned long long num_batches = 1;
        for( unsigned long long d=D; d<in->get_number_of_dimensions(); d++ ){
            num_batches *= in->get_size(d);
        }

        if( weak_greater(crop_offset+matrix_size_out, matrix_size_in) ){
            throw std::runtime_error( "crop: cropping size mismatch");;
        }

        const int num_elements_in = prod(matrix_size_in);
        const int num_elements_out = prod(matrix_size_out);

        T *in_ptr = in->get_data_ptr();
        T *out_ptr = out->get_data_ptr();

        for( unsigned long long frame_offset=0; frame_offset<num_batches; frame_offset++ ){
#ifdef USE_OMP
#pragma omp parallel for
#endif
            for( int idx=0; idx<num_elements_out; idx++ ){
                const typename uintd<D>::Type co = idx_to_co<D>( idx, matrix_size_out );
                const typename uintd<D>::Type co_os = crop_offset + co;
                const unsigned long long in_idx = co_to_idx<D>(co_os, matrix_size_in)+frame_offset*num_elements_in;
                out_ptr[idx+frame_offset*num_elements_out] = in_ptr[in_idx];
            }
        }
        return out;
    }    

    /**
     * @param[in] size Size of the output array
     * @param[in] in Input array
     * @param[in] val Value to use for padding
     * @returns New array of the specified size, containing the original input array in the center and val outside.
     */
    template<class T, unsigned long long D> boost::shared_ptr< hoNDArray<T> >
    pad( const typename uintd<D>::Type& size, hoNDArray<T> *in, T val = T(0) )
    {
        if( in == 0x0 ){
            throw std::runtime_error("pad: 0x0 array provided");;
        }

        if( in->get_number_of_dimensions() < D ){
            std::stringstream ss;
            ss << "pad: number of image dimensions should be at least " << D;
            throw std::runtime_error(ss.str());;
        }

        std::vector<unsigned long long> dims = to_std_vector(size);
        for( unsigned long long d=D; d<in->get_number_of_dimensions(); d++ ){
            dims.push_back(in->get_size(d));
        }
        boost::shared_ptr< hoNDArray<T> > out( new hoNDArray<T>(&dims) );

        typename uintd<D>::Type matrix_size_in = from_std_vector<unsigned long long,D>( *in->get_dimensions() );
        typename uintd<D>::Type matrix_size_out = from_std_vector<unsigned long long,D>( *out->get_dimensions() );

        unsigned long long num_batches = 1;
        for( unsigned long long d=D; d<in->get_number_of_dimensions(); d++ ){
            num_batches *= in->get_size(d);
        }

        if( weak_greater(matrix_size_in,matrix_size_out) ){
            throw std::runtime_error("pad: size mismatch, cannot expand");
        }

        const int num_elements_in = prod(matrix_size_in);
        const int num_elements_out = prod(matrix_size_out);
        const typename uintd<D>::Type offset = (matrix_size_out-matrix_size_in)>>1;

        T *in_ptr = in->get_data_ptr();
        T *out_ptr = out->get_data_ptr();

        for( unsigned long long frame_offset=0; frame_offset<num_batches; frame_offset++ ){
#ifdef USE_OMP
#pragma omp parallel for
#endif
            for( int idx=0; idx<num_elements_out; idx++ ){
                const typename uintd<D>::Type co_out = idx_to_co<D>( idx, matrix_size_out );
                T _out;
                bool inside = (co_out>=offset) && (co_out<(matrix_size_in+offset));

                if( inside )
                    _out = in_ptr[co_to_idx<D>( co_out-offset, matrix_size_in)+frame_offset*num_elements_in];
                else{
                    _out = val;
                }
                out_ptr[idx+frame_offset*num_elements_out] = _out;
            }
        }
        return out;
    }

    template<typename T> 
    bool permuteFirstTwoDimensions(const hoNDArray<T>& x, hoNDArray<T>& r)
    {
        try
        {
            unsigned long long NDim = x.get_number_of_dimensions();
            if ( NDim == 1 )
            {
                r = x;
                return true;
            }

            boost::shared_ptr< std::vector<unsigned long long> > dimX = x.get_dimensions();

            unsigned long long RO = x.get_size(0);
            unsigned long long E1 = x.get_size(1);
            unsigned long long numOfPermute =  x.get_number_of_elements()/(RO*E1);

            std::vector<unsigned long long> dimR(NDim);
            dimR = *dimX;
            dimR[0] = E1;
            dimR[1] = RO;

            if ( r.dimensions_equal(&dimR) )
            {
                r.create(dimR);
            }

            int n;

            #pragma omp parallel for default(none) private(n) shared(RO, E1, numOfPermute, x, r)
            for ( n=0; n<(int)numOfPermute; n++ )
            {
                const T* pX = x.begin() + n*RO*E1;
                T* pR = r.begin() + n*RO*E1;

                for ( unsigned long long e=0; e<E1; e++ )
                {
                    for ( unsigned long long r=0; r<RO; r++ )
                    {
                        pR[e+r*E1] = pX[r+e*RO];
                    }
                }
            }
        }
        catch (...)
        {
            GADGET_ERROR_MSG("Errors in permuteFirstTwoDimensions(const hoNDArray<T>& x, hoNDArray<T>& r) ... ");
            return false;
        }
        return true;
    }

    template<typename T> 
    bool permuteLastTwoDimensions(const hoNDArray<T>& x, hoNDArray<T>& r)
    {
        try
        {
            unsigned long long NDim = x.get_number_of_dimensions();
            if ( NDim == 1 )
            {
                r = x;
                return true;
            }

            boost::shared_ptr< std::vector<unsigned long long> > dimX = x.get_dimensions();

            unsigned long long lastDim = x.get_size(NDim-1);
            unsigned long long secondLastDim = x.get_size(NDim-2);
            unsigned long long N =  x.get_number_of_elements()/(lastDim*secondLastDim);

            std::vector<unsigned long long> dimR(NDim);
            dimR = *dimX;
            dimR[NDim-2] = lastDim;
            dimR[NDim-1] = secondLastDim;

            if ( !r.dimensions_equal(&dimR) )
            {
                r.create(dimR);
            }

            int l;

            #ifdef GCC_OLD_FLAG
                #pragma omp parallel for default(none) private(l) shared(lastDim, secondLastDim, N)
            #else
                #pragma omp parallel for default(none) private(l) shared(lastDim, secondLastDim, x, r, N)
            #endif
            for ( l=0; l<(int)lastDim; l++ )
            {
                for ( unsigned long long sl=0; sl<secondLastDim; sl++ )
                {
                    const T* pX = x.begin() + sl*N + l*N*secondLastDim;
                    T* pR = r.begin() + l*N + sl*N*lastDim;
                    memcpy(pR, pX, sizeof(T)*N);
                }
            }
        }
        catch (...)
        {
            GADGET_ERROR_MSG("Errors in permuteLastTwoDimensions(const hoNDArray<T>& x, hoNDArray<T>& r) ... ");
            return false;
        }
        return true;
    }

    /// copy the sub array x(:, indLastDim) to all other places of the last dimensions
    template<typename T> 
    bool repmatLastDimension(hoNDArray<T>& x, unsigned long long indLastDim)
    {
        try
        {
            unsigned long long NDim = x.get_number_of_dimensions();
            unsigned long long lastDim = x.get_size(NDim-1);
            GADGET_CHECK_RETURN_FALSE( indLastDim < lastDim );

            std::vector<unsigned long long> ind(NDim, 0);
            ind[NDim-1] = indLastDim;
            int offsetIndLastDim = x.calculate_offset(ind);

            unsigned long long N = x.get_number_of_elements() / lastDim;

            int l;
            #ifdef GCC_OLD_FLAG
                #pragma omp parallel for default(none) private(l) shared(lastDim, offsetIndLastDim, ind, indLastDim, N, NDim)
            #else
                #pragma omp parallel for default(none) private(l) shared(lastDim, offsetIndLastDim, x, ind, indLastDim, N, NDim)
            #endif
            for ( l=0; l<(int)lastDim; l++ )
            {
                if ( l==indLastDim ) continue;
                ind[NDim-1] = l;
                int offsetInd = x.calculate_offset(ind);

                memcpy(x.begin()+offsetInd, x.begin()+offsetIndLastDim, sizeof(T)*N);
            }
        }
        catch (...)
        {
            GADGET_ERROR_MSG("Errors in repmatLastDimension(hoNDArray<T>& x, unsigned long long indLastDim) ... ");
            return false;
        }
        return true;
    }

}
