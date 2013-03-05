#pragma once

#include "linearSolver.h"
#include "GadgetronException.h"
#include "real_utilities.h"
#include "complext.h"
#include <vector>
#include <iostream>

namespace Gadgetron{
/* Using adaptive step size from Zhou et al, 2006, Computational Optimization and Applications,
 * DOI: 10.1007/s10589-006-6446-0
 */

template <class ARRAY_TYPE> class gpSolver
: public linearSolver<ARRAY_TYPE> {
public:
	typedef typename ARRAY_TYPE::element_type ELEMENT_TYPE;
	typedef typename realType<ELEMENT_TYPE>::type REAL;


	virtual void set_domain_dimensions(std::vector<unsigned int> *dims ){
		for (int i = 0;  i < operators.size(); i++) operators[i]->set_domain_dimensions(dims);
	}
	virtual ~gpSolver(){}

	virtual bool add_regularization_operator(boost::shared_ptr< linearOperator<ARRAY_TYPE> > op ){
		return add_regularization_operator(op,2);
	}
	virtual bool add_regularization_operator(boost::shared_ptr< linearOperator<ARRAY_TYPE> > op, int L_norm ){
	  if (L_norm==1){
		  operators.push_back(boost::shared_ptr<gpRegularizationOperator>(new l1GPRegularizationOperator(op)));
	  }else{
		  operators.push_back(boost::shared_ptr<gpRegularizationOperator>(new l2GPRegularizationOperator(op)));
	  }
	  return true;
	}

	virtual bool add_regularization_operator(boost::shared_ptr< linearOperator<ARRAY_TYPE> > op, boost::shared_ptr<ARRAY_TYPE> prior, int L_norm=2 ){
	  if (L_norm==1){
		  operators.push_back(boost::shared_ptr<gpRegularizationOperator>(new l1GPRegularizationOperator(op,prior)));
	  }else{
		  operators.push_back(boost::shared_ptr<gpRegularizationOperator>(new l2GPRegularizationOperator(op,prior)));
	  }
	  return true;
	}

	virtual void add_regularization_group_operator ( boost::shared_ptr< linearOperator<ARRAY_TYPE> > op )
	{
	 current_group.push_back(op);
	}

	virtual void add_group(int L_norm=1)
	{
	  if(current_group.size()==0){
	  	BOOST_THROW_EXCEPTION(runtime_error( "Error: gpBBSolver::add_group : no regularization group operators added" ));
	  }
	  if (L_norm==2){
		for (int i =0; i < current_group.size(); i++){
			add_regularization_operator(current_group[i]);
		}

	  } else {

		  boost::shared_ptr<gpRegularizationOperator> new_group(new l1GroupGPRegularizationOperator(current_group));
		  operators.push_back(new_group);
	  }
	  current_group = std::vector<boost::shared_ptr<linearOperator<ARRAY_TYPE> > >();
	}

	virtual void add_group(boost::shared_ptr<ARRAY_TYPE> prior, int L_norm=1)
	{
	  if(current_group.size()==0){
		BOOST_THROW_EXCEPTION(runtime_error( "Error: gpBBSolver::add_group : no regularization group operators added" ));

	  }
	  if (L_norm==2){
		for (int i =0; i < current_group.size(); i++){
			add_regularization_operator(current_group[i],prior);
		}

	  } else {

		  boost::shared_ptr<gpRegularizationOperator> new_group(new l1GroupGPRegularizationOperator(current_group,prior));
		  operators.push_back(new_group);
	  }
	  current_group = std::vector<boost::shared_ptr<linearOperator<ARRAY_TYPE> > >();
	}

protected:


	virtual REAL add_gradient(ARRAY_TYPE* x, ARRAY_TYPE* g){
		REAL cost = REAL(0);

		for (int i = 0; i < operators.size(); i++){
			boost::shared_ptr<gpRegularizationOperator> op = operators[i];
			cost += op->add_gradient(x,g);

		}
		return cost;
	}



	class gpRegularizationOperator {
	public:
	  gpRegularizationOperator(){
	  }
	  gpRegularizationOperator(
			  boost::shared_ptr<ARRAY_TYPE> _prior){
		  prior = _prior;


	  }

	  virtual void set_domain_dimensions(std::vector<unsigned int> *dims)=0;

	  virtual REAL add_gradient(ARRAY_TYPE* x, ARRAY_TYPE* g)=0;

	protected:

	  boost::shared_ptr<ARRAY_TYPE> prior;

	};

	class l1GPRegularizationOperator : public gpRegularizationOperator {
	public:
	  l1GPRegularizationOperator(boost::shared_ptr<linearOperator<ARRAY_TYPE> > _op){

		  op = _op;
	  }
	  l1GPRegularizationOperator(
			  boost::shared_ptr<linearOperator<ARRAY_TYPE> > _op,
			  boost::shared_ptr<ARRAY_TYPE> _prior): gpRegularizationOperator(_prior){op = _op;}


	  virtual REAL add_gradient(ARRAY_TYPE* x, ARRAY_TYPE* g){
		  ARRAY_TYPE tmp(op->get_codomain_dimensions());
		  ARRAY_TYPE q(op->get_domain_dimensions());
		  ARRAY_TYPE* x2 = x;
		  if (this->prior.get()){
			  x2 = new ARRAY_TYPE;
			  *x2 = *x;
			  axpy(REAL(-1),this->prior.get(),x2);

		  }
		  op->mult_M(x2,&tmp);

		  REAL cost = op->get_weight()*asum(&tmp);
		  inplace_sgn(&tmp);
		  op->mult_MH(&tmp,&q,false);
		  axpy(op->get_weight(),&q,g);
		  if (this->prior.get()){ delete x2; }
		  return cost;
	  }
	  virtual void set_domain_dimensions(std::vector<unsigned int> *dims){
		  op->set_domain_dimensions(dims);
		  if (op->get_codomain_dimensions()->size() == 0){
			  std::cout << "WARNING: Codomain dimension not set. Setting to domain_dimension" << std::endl;
			  op->set_codomain_dimensions(dims);
		  }
	  }
	  boost::shared_ptr<linearOperator<ARRAY_TYPE> > op;


	};

	class l2GPRegularizationOperator : public gpRegularizationOperator {
	public:
	  l2GPRegularizationOperator(boost::shared_ptr<linearOperator<ARRAY_TYPE> > _op){
		  op = _op;

	  }
	  l2GPRegularizationOperator(boost::shared_ptr<linearOperator<ARRAY_TYPE> > _op,
			  boost::shared_ptr<ARRAY_TYPE> _prior): gpRegularizationOperator(_prior){
		  op = _op;

	  }
	  virtual REAL add_gradient(ARRAY_TYPE* x, ARRAY_TYPE* g){
		  ARRAY_TYPE tmp(op->get_codomain_dimensions());
		  ARRAY_TYPE q(op->get_domain_dimensions());
		  ARRAY_TYPE* x2 = x;
		  if (this->prior.get()){
			  x2 = new ARRAY_TYPE;
			  *x2 = *x;
			  axpy(REAL(-1),this->prior.get(),x2);

		  }
		  op->mult_M(x2,&tmp,false);

		  REAL cost = op->get_weight()*dot(&tmp,&tmp);
		  op->mult_MH(&tmp,&q,false);
		  axpy(op->get_weight(),&q,g);
		  if (this->prior.get()){ delete x2; }
		  return cost;
	  }
	  boost::shared_ptr<linearOperator<ARRAY_TYPE> > op;
	  virtual void set_domain_dimensions(std::vector<unsigned int> *dims){
		  op->set_domain_dimensions(dims);
		  if (op->get_codomain_dimensions()->size() == 0){
			  std::cout << "WARNING: Codomain dimension not set. Setting to domain_dimension" << std::endl;
			  op->set_codomain_dimensions(dims);
		  }
	  }
	};

	class l1GroupGPRegularizationOperator : public gpRegularizationOperator {
	 public:
	  l1GroupGPRegularizationOperator(std::vector<boost::shared_ptr<linearOperator<ARRAY_TYPE> > >_ops){
		  group = _ops;
		  threshold = REAL(1e-8);

	  }
	l1GroupGPRegularizationOperator(std::vector<boost::shared_ptr<linearOperator<ARRAY_TYPE> > >_ops,
			boost::shared_ptr<ARRAY_TYPE> _prior): gpRegularizationOperator(_prior){
		  group = _ops;
		  threshold = REAL(1e-8);

	  }
	  virtual REAL add_gradient(ARRAY_TYPE* x, ARRAY_TYPE* g){
		  std::vector<boost::shared_ptr<ARRAY_TYPE> > data;
		  ARRAY_TYPE gData(group.front()->get_codomain_dimensions());
		  gData.clear();


		  ARRAY_TYPE* x2 = x;
		  if (this->prior.get()){
			  x2 = new ARRAY_TYPE;
			  *x2 = *x;
			  axpy(REAL(-1),this->prior.get(),x2);

		  }


		  for (int i = 0; i < group.size(); i++ ){
			  boost::shared_ptr<linearOperator<ARRAY_TYPE> > op = group[i];
			  boost::shared_ptr<ARRAY_TYPE> tmp(new ARRAY_TYPE(op->get_codomain_dimensions().get()));
			  op->mult_M(x2,tmp.get());

			  data.push_back(tmp);

			  ARRAY_TYPE tmp2 = *tmp;
			  tmp2 *= *tmp;
			  gData += tmp2;

		  }
		  if (this->prior.get()){
		  			  delete x2;
		  }
		  gData.sqrt();
		  REAL cost = group.front()->get_weight()*asum(&gData);
		  clamp_min(&gData,threshold);
		  gData.reciprocal();

		  ARRAY_TYPE q(group.front()->get_domain_dimensions());



		  for (int i = 0; i < group.size(); i++ ){
			  boost::shared_ptr<linearOperator<ARRAY_TYPE> > op = group[i];
			  boost::shared_ptr<ARRAY_TYPE> tmp = data[i];
			  *tmp *= gData;


			  op->mult_MH(tmp.get(),&q,false);
			  axpy(op->get_weight(),&q,g);
		  }

		  return cost;
	  }

	  void set_threshold(REAL _threshold){
		  threshold = _threshold;
	  }

	 protected:
	  std::vector<boost::shared_ptr<linearOperator<ARRAY_TYPE> > > group;
	  REAL threshold;

	  virtual void set_domain_dimensions(std::vector<unsigned int> *dims){

		  for (int i = 0; i < group.size(); i++ ){
			  boost::shared_ptr<linearOperator<ARRAY_TYPE> > op = group[i];
			  op->set_domain_dimensions(dims);
			  if (op->get_codomain_dimensions()->size() == 0){
				  std::cout << "WARNING: Codomain dimension not set. Setting to domain_dimension" << std::endl;
				  op->set_codomain_dimensions(dims);
			  }
		  }
	  }


	 };

	std::vector< boost::shared_ptr<gpRegularizationOperator> > operators;
	std::vector< boost::shared_ptr<linearOperator<ARRAY_TYPE> > >  current_group;


};
}
