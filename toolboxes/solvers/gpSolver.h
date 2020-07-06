#pragma once

#include "linearOperatorSolver.h"
#include "real_utilities.h"
#include "complext.h"

#include <vector>
#include <iostream>

namespace Gadgetron{

  /* Using adaptive step size from Zhou et al, 2006, Computational Optimization and Applications,
   * DOI: 10.1007/s10589-006-6446-0
   */

  template <class ARRAY_TYPE> class gpSolver : public linearOperatorSolver<ARRAY_TYPE>
  {
  protected:
    typedef typename ARRAY_TYPE::element_type ELEMENT_TYPE;
    typedef typename realType<ELEMENT_TYPE>::Type REAL;

  public:

    virtual void set_domain_dimensions(std::vector<size_t> *dims ){
      for (int i = 0;  i < operators.size(); i++) operators[i]->set_domain_dimensions(dims);
    }
    virtual ~gpSolver(){}

    virtual void add_nonlinear_operator(boost::shared_ptr< generalOperator<ARRAY_TYPE> > op ){
      operators.push_back(op);
    }

    virtual void add_regularization_operator(boost::shared_ptr< linearOperator<ARRAY_TYPE> > op ){
      add_regularization_operator(op,2);
    }
    virtual void add_regularization_operator(boost::shared_ptr< linearOperator<ARRAY_TYPE> > op, int L_norm ){
      if (L_norm==1){
        operators.push_back(boost::shared_ptr<gpRegularizationOperator>(new l1GPRegularizationOperator(op)));
      }else{
        operators.push_back(op);
      }
    }

    virtual void add_regularization_operator(boost::shared_ptr< linearOperator<ARRAY_TYPE> > op, boost::shared_ptr<ARRAY_TYPE> prior, int L_norm=2 ){
      if (L_norm==1){
        operators.push_back(boost::shared_ptr<gpRegularizationOperator>(new l1GPRegularizationOperator(op,prior)));
      }else{
        operators.push_back(boost::shared_ptr<gpRegularizationOperator>(new l2GPRegularizationOperator(op,prior)));
      }
    }

    virtual void add_regularization_group_operator ( boost::shared_ptr< linearOperator<ARRAY_TYPE> > op )
    {
      current_group.push_back(op);
    }

    virtual void add_group(int L_norm=1)
    {
      if(current_group.size()==0){
        throw std::runtime_error( "Error: gpBBSolver::add_group : no regularization group operators added" );
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
        throw std::runtime_error( "Error: gpBBSolver::add_group : no regularization group operators added" );

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


    virtual void add_gradient(ARRAY_TYPE* x, ARRAY_TYPE* g){
      for (int i = 0; i < operators.size(); i++){
        boost::shared_ptr<generalOperator<ARRAY_TYPE> > op = operators[i];
        op->gradient(x,g,true);
      }

    }



    class gpRegularizationOperator : public generalOperator<ARRAY_TYPE> {
    public:
      gpRegularizationOperator() : generalOperator<ARRAY_TYPE>(){
      }

      gpRegularizationOperator(std::vector<size_t> *dims): generalOperator<ARRAY_TYPE>(){this->set_domain_dimensions(dims);};
      gpRegularizationOperator(
                               boost::shared_ptr<ARRAY_TYPE> _prior): generalOperator<ARRAY_TYPE>(){
        prior = _prior;
      }

      gpRegularizationOperator(boost::shared_ptr<ARRAY_TYPE> _prior,std::vector<size_t> *dims): generalOperator<ARRAY_TYPE>(){
        prior = _prior;
        set_domain_dimensions(dims);
      }



    protected:

      boost::shared_ptr<ARRAY_TYPE> prior;

    };


    class l2GPRegularizationOperator : public gpRegularizationOperator {
    public:
      l2GPRegularizationOperator(boost::shared_ptr<linearOperator<ARRAY_TYPE> > _op){
        op = _op;
      }
      l2GPRegularizationOperator(
                                 boost::shared_ptr<linearOperator<ARRAY_TYPE> > _op,
                                 boost::shared_ptr<ARRAY_TYPE> _prior): gpRegularizationOperator(_prior){op = _op;}
      virtual void gradient(ARRAY_TYPE* x, ARRAY_TYPE* g,bool accumulate=false){
        ARRAY_TYPE* x2 = x;
        if (this->prior.get()){
          x2 = new ARRAY_TYPE(*x);
          axpy(ELEMENT_TYPE(-1),this->prior.get(),x2);
        }
        op->gradient(x,g,accumulate);
      }

      virtual REAL magnitude(ARRAY_TYPE* x){
        ARRAY_TYPE tmp(op->get_codomain_dimensions());
        ARRAY_TYPE* x2 = x;
        if (this->prior.get()){
          x2 = new ARRAY_TYPE(*x);
          *x2 -= *this->prior;
        }
        else clear(&tmp);
        op->mult_M(x2,&tmp);
        if (this->prior.get()) delete x2;
        return std::sqrt(op->get_weight())*real(dot(&tmp,&tmp));
      }
    protected:
      boost::shared_ptr<linearOperator<ARRAY_TYPE> > op;
    };

    class l1GPRegularizationOperator : public gpRegularizationOperator {
    public:
      l1GPRegularizationOperator(boost::shared_ptr<linearOperator<ARRAY_TYPE> > _op){

        op = _op;
      }
      l1GPRegularizationOperator(
                                 boost::shared_ptr<linearOperator<ARRAY_TYPE> > _op,
                                 boost::shared_ptr<ARRAY_TYPE> _prior): gpRegularizationOperator(_prior){op = _op;}


      virtual void gradient(ARRAY_TYPE* x, ARRAY_TYPE* g, bool accumulate=false){
        ARRAY_TYPE tmp(op->get_codomain_dimensions());
        ARRAY_TYPE q(op->get_domain_dimensions());
        ARRAY_TYPE* x2 = x;

        if (!accumulate) clear(g);
        if (this->prior.get()){
          x2 = new ARRAY_TYPE;
          *x2 = *x;
          axpy(ELEMENT_TYPE(-1),this->prior.get(),x2);

        }
        op->mult_M(x2,&tmp);
        sgn_inplace(&tmp);
        op->mult_MH(&tmp,&q,false);
        axpy(ELEMENT_TYPE(op->get_weight()),&q,g);
        if (this->prior.get()) delete x2;
      }

      virtual REAL magnitude(ARRAY_TYPE* x){
        ARRAY_TYPE tmp(op->get_codomain_dimensions());
        ARRAY_TYPE* x2 = x;
        if (this->prior.get()){
          x2 = new ARRAY_TYPE(*x);
          *x2 -= *this->prior;
        }
        op->mult_M(x2,&tmp);
        if (this->prior.get()) delete x2;
        return op->get_weight()*asum(&tmp);
      }


      virtual void set_domain_dimensions(std::vector<size_t> *dims){
        generalOperator<ARRAY_TYPE>::set_domain_dimensions(dims);
        op->set_domain_dimensions(dims);
        if (op->get_codomain_dimensions()->size() == 0){
          GDEBUG_STREAM("WARNING: Codomain dimension not set. Setting to domain_dimension" << std::endl);
          op->set_codomain_dimensions(dims);
        }
      }
      boost::shared_ptr<linearOperator<ARRAY_TYPE> > op;
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
      virtual void gradient(ARRAY_TYPE* x, ARRAY_TYPE* g,bool accumulate=false){
        std::vector<boost::shared_ptr<ARRAY_TYPE> > data;
        ARRAY_TYPE gData(group.front()->get_codomain_dimensions());
        clear(&gData);

        if (!accumulate) clear(g);
        ARRAY_TYPE* x2 = x;
        if (this->prior.get()){
          x2 = new ARRAY_TYPE(*x);
          *x2 -= *this->prior;
        }

        for (int i = 0; i < group.size(); i++ ){
          boost::shared_ptr<linearOperator<ARRAY_TYPE> > op = group[i];
          boost::shared_ptr<ARRAY_TYPE> tmp(new ARRAY_TYPE(op->get_codomain_dimensions().get()));
          op->mult_M(x2,tmp.get());
          data.push_back(tmp);
          ARRAY_TYPE tmp2 = *tmp;
          tmp2 *= *tmp; //Square data
          gData += tmp2;
        }
        if (this->prior.get()){
          delete x2;
        }
        sqrt_inplace(&gData);
        //REAL cost = group.front()->get_weight()*asum(&gData);
        clamp_min(&gData,threshold);
        reciprocal_inplace(&gData);

        ARRAY_TYPE q(group.front()->get_domain_dimensions());

        for (int i = 0; i < group.size(); i++ ){
          boost::shared_ptr<linearOperator<ARRAY_TYPE> > op = group[i];
          boost::shared_ptr<ARRAY_TYPE> tmp = data[i];
          *tmp *= gData;
          op->mult_MH(tmp.get(),&q,false);
          axpy(ELEMENT_TYPE(op->get_weight()),&q,g);
        }
      }


      virtual REAL magnitude(ARRAY_TYPE* x){
        ARRAY_TYPE gData(group.front()->get_codomain_dimensions());
        clear(&gData);
        ARRAY_TYPE* x2 =x;
        if (this->prior.get()){
          x2 = new ARRAY_TYPE(*x);
          *x2 -= *this->prior;
        }
        for (int i = 0; i < group.size(); i++ ){
          boost::shared_ptr<linearOperator<ARRAY_TYPE> > op = group[i];
          ARRAY_TYPE tmp(op->get_codomain_dimensions().get());
          op->mult_M(x2,&tmp);
          tmp *= tmp;
          gData += tmp;
        }
        if (this->prior.get()) delete x2;
        sqrt_inplace(&gData);
        return group.front()->get_weight()*asum(&gData);
      }

      void set_threshold(REAL _threshold){
        threshold = _threshold;
      }

    protected:

      std::vector<boost::shared_ptr<linearOperator<ARRAY_TYPE> > > group;
      REAL threshold;

      virtual void set_domain_dimensions(std::vector<size_t> *dims){
        generalOperator<ARRAY_TYPE>::set_domain_dimensions(dims);
        for (int i = 0; i < group.size(); i++ ){
          boost::shared_ptr<linearOperator<ARRAY_TYPE> > op = group[i];
          op->set_domain_dimensions(dims);
          if (op->get_codomain_dimensions()->size() == 0){
            GDEBUG_STREAM("WARNING: Codomain dimension not set. Setting to domain_dimension" << std::endl);
            op->set_codomain_dimensions(dims);
          }
        }
      }
    };

    std::vector< boost::shared_ptr< generalOperator<ARRAY_TYPE> > > operators;
    std::vector< boost::shared_ptr<linearOperator<ARRAY_TYPE> > >  current_group;
  };
}
