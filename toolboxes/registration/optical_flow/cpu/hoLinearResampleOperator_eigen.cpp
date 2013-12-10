#include "hoLinearResampleOperator_eigen.h"
#include "vector_td_utilities.h"
#include "vector_td_operators.h"

#include "GadgetronTimer.h"

#include <stdio.h>
#include <algorithm>
#include <Eigen/Core>

namespace Gadgetron{

    template <class T, unsigned int D> void
        hoLinearResampleOperator_eigen<T,D>::mult_M( hoNDArray<T> *in, hoNDArray<T> *out, bool accumulate )
    {
        if( !this->preprocessed_ ){
            throw std::runtime_error("hoLinearResampleOperator::mult_M(): displacements not set." );
        }

        if( !in || !in->get_data_ptr() || !out || !out->get_data_ptr() ){
            throw std::runtime_error("hoLinearResampleOperator::mult_M(): illegal input/output array." );
        }

        Eigen::Map< Eigen::Matrix<typename realType<T>::Type, 1, Eigen::Dynamic> > in_vec( in->get_data_ptr(), in->get_number_of_elements() );
        Eigen::Map< Eigen::Matrix<typename realType<T>::Type, 1, Eigen::Dynamic> > out_vec( out->get_data_ptr(), out->get_number_of_elements() );

        out_vec = in_vec * (*R_);
    }

    template <class T, unsigned int D> void
        hoLinearResampleOperator_eigen<T,D>::mult_MH( hoNDArray<T> *in, hoNDArray<T> *out, bool accumulate )
    {
        if( !this->preprocessed_ ){
            throw std::runtime_error("hoLinearResampleOperator::mult_M(): displacements not set." );
        }

        if( !in || !in->get_data_ptr() || !out || !out->get_data_ptr() ){
            throw std::runtime_error("hoLinearResampleOperator::mult_M(): illegal input/output array." );
        }

        Eigen::Map< Eigen::Matrix<typename realType<T>::Type, Eigen::Dynamic, 1> > in_vec( in->get_data_ptr(), in->get_number_of_elements() );
        Eigen::Map< Eigen::Matrix<typename realType<T>::Type, Eigen::Dynamic, 1> > out_vec( out->get_data_ptr(), out->get_number_of_elements() );

        out_vec = (*R_) * in_vec;
    }

    template <class T, unsigned int D> void
        hoLinearResampleOperator_eigen<T,D>::set_displacement_field( boost::shared_ptr< hoNDArray<typename realType<T>::Type> > displacements )
    {
        if( displacements.get() == 0x0 ){
            throw std::runtime_error("hoLinearResampleOperator_eigen_eigen::set_displacement_field : displacements ptr is 0x0." );
        }  

        const int surplus = displacements->get_number_of_dimensions()-D;

        if( !( surplus == 1 || surplus == 2 ) ){
            throw std::runtime_error("hoLinearResampleOperator_eigen::set_displacement_field : unexpected array dimensionality." );
        }  

        // Determine the number of registrations performed
        const size_t extended_dim = (surplus == 1) ? 1 : displacements->get_size(D); 
        temporal_dim_size_ = extended_dim;

        const size_t field_dim = (surplus == 1) ? displacements->get_size(D) : displacements->get_size(D+1);

        if( !(field_dim == D || field_dim == D+1 )){
            throw std::runtime_error("hoLinearResampleOperator_eigen::set_displacement_field : illegal tailing array dim" );
        }

        const typename uint64d<D>::Type matrix_size = from_std_vector<size_t,D>( *(displacements->get_dimensions()));

        const size_t num_elements_mat = prod(matrix_size);
        const size_t num_elements_ext = prod(matrix_size)*extended_dim;

        R_ = boost::shared_ptr< Eigen::SparseMatrix<typename realType<T>::Type> >
            ( new Eigen::SparseMatrix<typename realType<T>::Type>( num_elements_mat, num_elements_ext ) );

        std::vector< Eigen::Triplet<typename realType<T>::Type> > coefficients;

        for( size_t idx=0; idx<num_elements_ext; idx++ ){

            const size_t batch_no = idx/num_elements_mat;
            const size_t idx_in_batch = idx-batch_no*num_elements_mat;

            const typename uint64d<D>::Type co = idx_to_co<D>( idx_in_batch, matrix_size );

            typename reald<typename realType<T>::Type,D>::Type co_disp = to_reald<typename realType<T>::Type,size_t,D>(co);
            for( size_t dim=0; dim<D; dim++ ){
                typename realType<T>::Type tmp = displacements->get_data_ptr()[dim*num_elements_ext+batch_no*num_elements_mat+idx_in_batch];
                co_disp.vec[dim] += tmp;
            } 

            // Determine the number of neighbors
            //

            const typename uint64d<D>::Type twos = to_vector_td<size_t,D>(2);
            const size_t num_neighbors = this->get_num_neighbors();

            // Weights are non-zero only if all neighbors exist
            //

            if( this->is_border_pixel(co_disp, matrix_size) )
                continue;

            // Iterate over all neighbors
            //

            //
            // Eigen asks us to build the matrix column by column 
            // It is more easy then to construct the transpose
            //

            size_t mat_j = idx;
            size_t mat_i;

            for( size_t i=0; i<num_neighbors; i++ ){

                // Determine image coordinate of current neighbor
                //

                const typename uint64d<D>::Type stride = idx_to_co<D>( i, twos );

                if( weak_greater_equal( stride, matrix_size ) ) continue; // For dimensions of size 1

                typename reald<typename realType<T>::Type,D>::Type co_stride;

                for( size_t dim=0; dim<D; dim++ ){
                    if( stride.vec[dim] == 0 ){
                        co_stride.vec[dim] = std::floor(co_disp.vec[dim]);
                    }
                    else{
                        co_stride.vec[dim] = std::ceil(co_disp.vec[dim]);
                        if( co_stride.vec[dim] == co_disp.vec[dim] )
                            co_stride.vec[dim] += typename realType<T>::Type(1.0);
                    }
                }

                // Validate that the coordinate is within the expected range
                //

                typename uint64d<D>::Type ones = to_vector_td<size_t,D>(1);
                typename uint64d<D>::Type co_stride_uint64d = to_uint64d<typename realType<T>::Type,D>(co_stride);

                if( weak_greater( co_stride_uint64d, matrix_size-ones ) ){

                    for( size_t dim=0; dim<D; dim++ ){
                        if( co_stride[dim] < typename realType<T>::Type(0) )
                            co_stride_uint64d[dim] = 0;
                        if( co_stride[dim] > (typename realType<T>::Type(matrix_size[dim])-typename realType<T>::Type(1)) )
                            co_stride_uint64d[dim] = matrix_size[dim]-1;
                    }
                }

                mat_i = co_to_idx<D>(co_stride_uint64d, matrix_size);

                // Determine weight
                //

                typename realType<T>::Type weight = typename realType<T>::Type(1);

                for( size_t dim=0; dim<D; dim++ ){	  
                    if( stride.vec[dim] == 0 ){
                        weight *= (typename realType<T>::Type(1.0)-(co_disp.vec[dim]-co_stride.vec[dim])); }
                    else{
                        weight *= (typename realType<T>::Type(1.0)-(co_stride.vec[dim]-co_disp.vec[dim])); }
                }

                // Insert weight in resampling matrix R_
                //

                //R_->insert( mat_i, mat_j ) =  weight;
                coefficients.push_back(Eigen::Triplet<typename realType<T>::Type>(mat_i, mat_j, weight));
            }
        }  
        //R_->finalize();
        R_->setFromTriplets(coefficients.begin(), coefficients.end());
        this->preprocessed_ = true;
    }

    template <class T, unsigned int D> bool
        hoLinearResampleOperator_eigen<T,D>::is_border_pixel( typename reald<typename realType<T>::Type,D>::Type co, typename uint64d<D>::Type dims )
    {
        for( size_t dim=0; dim<D; dim++ ){
            if( dims[dim] > 1 && ( co[dim] < typename realType<T>::Type(0) || co[dim] >= (typename realType<T>::Type(dims[dim])-typename realType<T>::Type(1)) ) )
                return true;
        }
        return false;
    }

    template <class T, unsigned int D> size_t
        hoLinearResampleOperator_eigen<T,D>::get_num_neighbors()
    {
        return 1 << D;
    }


    template class EXPORTCPUREG hoLinearResampleOperator_eigen<float,1>;
    template class EXPORTCPUREG hoLinearResampleOperator_eigen<float,2>;
    template class EXPORTCPUREG hoLinearResampleOperator_eigen<float,3>;
    template class EXPORTCPUREG hoLinearResampleOperator_eigen<float,4>;

    template class EXPORTCPUREG hoLinearResampleOperator_eigen<double,1>;
    template class EXPORTCPUREG hoLinearResampleOperator_eigen<double,2>;
    template class EXPORTCPUREG hoLinearResampleOperator_eigen<double,3>;
    template class EXPORTCPUREG hoLinearResampleOperator_eigen<double,4>;
}
