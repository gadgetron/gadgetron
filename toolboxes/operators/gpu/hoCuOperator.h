#pragma once
#include "hoCuNDArray_math.h"
#include "linearOperator.h"

#include <boost/shared_ptr.hpp>
namespace Gadgetron{


template<class T > class hoCuOperator : public linearOperator<hoCuNDArray<T> > {

	public:
		hoCuOperator(){};
		hoCuOperator(boost::shared_ptr<linearOperator<hoNDArray<T> > > _op): op(_op) {};
		virtual ~hoCuOperator(){};

		virtual void mult_M(hoCuNDArray<T>* in, hoCuNDArray<T>* out, bool accumulate=false){
			op->mult_M(in,out,accumulate);
		}
		virtual void mult_MH(hoCuNDArray<T>* in, hoCuNDArray<T>* out, bool accumulate=false){
			op->mult_MH(in,out,accumulate);
		}

		virtual void gradient(hoCuNDArray<T>* in, hoCuNDArray<T>* out, bool accumulate=false){
			op->gradient(in,out,accumulate);
		}
		virtual void mult_MH_M(hoCuNDArray<T>* in, hoCuNDArray<T>* out, bool accumulate=false){
			op->mult_MH_M(in,out,accumulate);
		}

		virtual boost::shared_ptr< linearOperator< hoCuNDArray<T> > > clone() {
			return linearOperator< hoCuNDArray<T> >::clone(this);
		}
		virtual boost::shared_ptr< std::vector<unsigned int> > get_codomain_dimensions(){
			return op->get_codomain_dimensions();
		}
		virtual boost::shared_ptr< std::vector<unsigned int> > get_domain_dimensions(){
			return op->get_domain_dimensions();
		}
		 virtual void set_weight( typename realType<T>::Type weight ){ op->set_weight(weight); };
		virtual typename realType<T>::Type get_weight(){ return op->get_weight(); };
		virtual void set_codomain_dimensions( std::vector<unsigned int> *dims ){
			op->set_codomain_dimensions(dims);
		}
		virtual void set_domain_dimensions( std::vector<unsigned int> *dims ){
			op->set_domain_dimensions(dims);
		}
	private:
	 boost::shared_ptr<linearOperator<hoNDArray<T> > > op;
};

template<class T> boost::shared_ptr<linearOperator<hoCuNDArray<T> > > to_hoCu(boost::shared_ptr<linearOperator<hoNDArray<T> > > _op){
	return boost::shared_ptr<hoCuOperator<T> > (new hoCuOperator<T>(_op));
}

}
