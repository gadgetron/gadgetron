#pragma once

#include "cuTV1DOperator.h"
#include "hoNDArray.h"
#include "sparsifyingOperator.h"
#include "check_CUDA.h"

template<class REAL, class T, unsigned int D> class hoTV1DOperator : public sparsifyingOperator< REAL, hoNDArray<T> >{
public:
	hoTV1DOperator() : sparsifyingOperator<REAL, hoNDArray<T> >(), cuTV(){
		limit_ = REAL(1e-8);


		cuTV.set_limit(limit_);
	}
	virtual ~hoTV1DOperator(){}
	void set_limit(REAL limit){
		limit_ = limit;
		cuTV.set_limit(limit);
	}
	virtual void apply(hoNDArray<T>* in ,hoNDArray<T>* out, bool accumulate=false){
		const vector_td<unsigned int,D> dims = from_std_vector<unsigned int, D>(*(in->get_dimensions()));
		int elements = in->get_number_of_elements();

		for (int i =0; i < (elements/prod(dims)); i++){
			hoNDArray<T> tmp_in;
			std::vector<unsigned int> dimensions = to_std_vector(dims);
			tmp_in.create(&dimensions,in->get_data_ptr()+i*prod(dims));
			hoNDArray<T> tmp_out;
			tmp_out.create(&dimensions,out->get_data_ptr()+i*prod(dims));
			cuNDArray<T> cuIn(&tmp_in);
			cuNDArray<T> cuOut(&tmp_out);

			cuTV.apply(&cuIn,&cuOut,accumulate);
			boost::shared_ptr< hoNDArray<T> > tmp = cuOut.to_host();
			tmp_out = *tmp;
		}

	}
	virtual void gradient(hoNDArray<T>* in ,hoNDArray<T>* out, bool accumulate=false){
		/*if (dim_ == 3){
			unsigned int k[] = {0,1,3,2};
			std::vector<unsigned int> p(k,k+4);
			std::vector<unsigned int> in_dims = *in->get_dimensions();
			int n[] = {in_dims[0],in_dims[1],in_dims[3],in_dims[2]};
			std::vector<unsigned int> p_dims(n,n+4);
			hoNDArray<T> p_in;
			p_in.create(&p_dims);
			in->permute(&p,&p_in);
			std::vector<unsigned int> dimensions = p_dims;
			dimensions[3]=1;
			int slice_elems = p_dims[0]*p_dims[1]*p_dims[2];
			hoNDArray<T> p_out;
				p_out.create(&p_dims);
			if (accumulate){
				out->permute(&p,&p_in);
			} else {
				hoNDA_clear(&p_out);
			}
			for (int i =0; i < p_dims[3]; i++){

				hoNDArray<T> tmp_in;
				tmp_in.create(&dimensions,p_in.get_data_ptr()+i*slice_elems);
				hoNDArray<T> tmp_out;
				tmp_out.create(&dimensions,p_out.get_data_ptr()+i*slice_elems);
				cuNDArray<T> cuIn(&tmp_in);
				cuNDArray<T> cuOut(&tmp_out);

				cuTV.gradient(&cuIn,&cuOut,accumulate);
				boost::shared_ptr< hoNDArray<T> > tmp = cuOut.to_host();
				tmp_out = *tmp;
			}
			p_out.permute(&p,out);


		} else {*/
		const vector_td<unsigned int,D> dims = from_std_vector<unsigned int, D>(*(in->get_dimensions()));
		int elements = in->get_number_of_elements();

		for (int i =0; i < (elements/prod(dims)); i++){
			std::vector<unsigned int> dimensions = to_std_vector(dims);
			hoNDArray<T> tmp_in;
			tmp_in.create(&dimensions,in->get_data_ptr()+i*prod(dims));
			hoNDArray<T> tmp_out;
			tmp_out.create(&dimensions,out->get_data_ptr()+i*prod(dims));
			cuNDArray<T> cuIn(&tmp_in);
			cuNDArray<T> cuOut(&tmp_out);

			cuTV.gradient(&cuIn,&cuOut,accumulate);
			boost::shared_ptr< hoNDArray<T> > tmp = cuOut.to_host();
			tmp_out = *tmp;
		}
		//}

	}

	virtual void set_weight(REAL weight){
		this->weight_ = weight;
		cuTV.set_weight(weight);
	}


protected:
	REAL limit_;
	cuTV1DOperator<REAL,T,D> cuTV;

};



template<class REAL, class T> class hoTV1DOperator<REAL,T,4> : public sparsifyingOperator< REAL, hoNDArray<T> >{
public:
	hoTV1DOperator() : sparsifyingOperator<REAL, hoNDArray<T> >(), cuTV(){
		limit_ = REAL(1e-8);


		cuTV.set_limit(limit_);
	}
	virtual ~hoTV1DOperator(){}
	void set_limit(REAL limit){
		limit_ = limit;
		cuTV.set_limit(limit);
	}
	virtual void apply(hoNDArray<T>* in ,hoNDArray<T>* out, bool accumulate=false){
		/*const vector_td<unsigned int,D> dims = from_std_vector<unsigned int, D>(*(in->get_dimensions()));
		int elements = in->get_number_of_elements();

		for (int i =0; i < (elements/prod(dims)); i++){
			hoNDArray<T> tmp_in;
			std::vector<unsigned int> dimensions = to_std_vector(dims);
			tmp_in.create(&dimensions,in->get_data_ptr()+i*prod(dims));
			hoNDArray<T> tmp_out;
			tmp_out.create(&dimensions,out->get_data_ptr()+i*prod(dims));
			cuNDArray<T> cuIn(&tmp_in);
			cuNDArray<T> cuOut(&tmp_out);

			cuTV.apply(&cuIn,&cuOut,accumulate);
			boost::shared_ptr< hoNDArray<T> > tmp = cuOut.to_host();
			tmp_out = *tmp;
		}*/

	}
	virtual void gradient(hoNDArray<T>* in ,hoNDArray<T>* out, bool accumulate=false){
		std::vector<unsigned int> in_dims = *in->get_dimensions();

		int slice_size = in_dims[0]*in_dims[1];
		unsigned int tmp[] = {in_dims[0],in_dims[1],in_dims[3]};
		std::vector<unsigned int> cu_dims(tmp,tmp+3);
		cuNDArray<T> cu_in;
		cu_in.create(&cu_dims);
		T* cu_in_ptr = cu_in.get_data_ptr();
		cuNDArray<T> cu_out;
		cu_out.create(&cu_dims);
		T* cu_out_ptr = cu_out.get_data_ptr();

		T* in_ptr = in->get_data_ptr();
		T* out_ptr = out->get_data_ptr();

		for (int iz = 0; iz < in_dims[2]; iz++){
			for (int it =0; it < in_dims[3]; it++){
				cudaMemcpy(cu_in_ptr+it*slice_size,
						in_ptr+it*slice_size*in_dims[2]+iz*slice_size, slice_size*sizeof(T), cudaMemcpyHostToDevice);
				if (accumulate){
					cudaMemcpy(cu_out_ptr+it*slice_size,
							out_ptr+it*slice_size*in_dims[2]+iz*slice_size,
							slice_size*sizeof(T), cudaMemcpyHostToDevice);
				}
				CHECK_FOR_CUDA_ERROR();
			}
			cuTV.gradient(&cu_in,&cu_out,accumulate);
			for (int it =0; it < in_dims[3]; it++){
				cudaMemcpy(out_ptr+it*slice_size*in_dims[2]+iz*slice_size,
						cu_out_ptr+it*slice_size, slice_size*sizeof(T), cudaMemcpyDeviceToHost);
				CHECK_FOR_CUDA_ERROR();
			}
		}





	}

	virtual void set_weight(REAL weight){
		this->weight_ = weight;
		cuTV.set_weight(weight);
	}


protected:
	REAL limit_;
	cuTV1DOperator<REAL,T,3> cuTV;

};

