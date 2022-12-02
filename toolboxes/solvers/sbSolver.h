/*
  An implementation of the "Generalized Split Bregman Algorithm" - sec. 3.2. of the paper
  "The Split Bregman Method for L1-Regularized Problems" by Tom Goldstein and Stanley Osher. 
  Siam J. Imaging Sciences. Vol. 2, No. 2, pp. 323-343.
 */

#pragma once

#include "linearOperatorSolver.h"
#include "vector_td_utilities.h"
#include "encodingOperatorContainer.h"
#include "identityOperator.h"

#include <vector>
#include <iostream>
#include <set>

namespace Gadgetron{

template< class ARRAY_TYPE_REAL,
class ARRAY_TYPE_ELEMENT,
class INNER_SOLVER >

class sbSolver : public linearOperatorSolver<ARRAY_TYPE_ELEMENT>
{

protected:

	typedef typename ARRAY_TYPE_REAL::element_type REAL;
	typedef typename ARRAY_TYPE_ELEMENT::element_type ELEMENT_TYPE;

	class sbRegularizationOperator{

	public:

		sbRegularizationOperator() {}
		sbRegularizationOperator(boost::shared_ptr< linearOperator<ARRAY_TYPE_ELEMENT> > op) { reg_op=op; }
		virtual ~sbRegularizationOperator(){}

		virtual void initialize(boost::shared_ptr< std::vector<size_t> > image_dims, REAL normalization_factor = REAL(1))
		{
			d_k = boost::shared_ptr<ARRAY_TYPE_ELEMENT>(new ARRAY_TYPE_ELEMENT(*reg_op->get_codomain_dimensions()));
			b_k = boost::shared_ptr<ARRAY_TYPE_ELEMENT>(new ARRAY_TYPE_ELEMENT(*reg_op->get_codomain_dimensions()));
			clear(d_k.get());
			clear(b_k.get());
			if(prior.get()){
				p_M = boost::shared_ptr<ARRAY_TYPE_ELEMENT>(new ARRAY_TYPE_ELEMENT(*reg_op->get_codomain_dimensions()));
				reg_op->mult_M(prior.get(),p_M.get());
				*p_M *= normalization_factor;
			}
		}

		virtual void update_encoding_space(ARRAY_TYPE_ELEMENT* encoding_space)
		{
			*encoding_space = *d_k;
			*encoding_space -= *b_k;
			if(prior.get())
				*encoding_space += *p_M;
		}

		virtual void deinitialize()
		{
			d_k.reset();
			b_k.reset();
			p_M.reset();
		}

		REAL get_weight(){ return reg_op->get_weight(); }
		void set_weight(REAL weight){ reg_op->set_weight(weight); }

		virtual void update_dk(ARRAY_TYPE_ELEMENT*) = 0;
		virtual void update_dk_bk(ARRAY_TYPE_ELEMENT*) = 0;

		virtual boost::shared_ptr< std::vector<size_t> > get_codomain_dimensions(){
			return reg_op->get_codomain_dimensions();
		}

		virtual void set_prior(boost::shared_ptr<ARRAY_TYPE_ELEMENT> image){ prior=image; }

		boost::shared_ptr< linearOperator< ARRAY_TYPE_ELEMENT> > reg_op;
		boost::shared_ptr<ARRAY_TYPE_ELEMENT> d_k;
		boost::shared_ptr<ARRAY_TYPE_ELEMENT> b_k;
		boost::shared_ptr<ARRAY_TYPE_ELEMENT> p_M;
		boost::shared_ptr<ARRAY_TYPE_ELEMENT> prior;
	};

	class sbL1RegularizationOperator : public sbRegularizationOperator
	{
	public:
		sbL1RegularizationOperator(boost::shared_ptr< linearOperator< ARRAY_TYPE_ELEMENT> > op) : sbRegularizationOperator(op) {}

		virtual void update_dk(ARRAY_TYPE_ELEMENT* u_k)
		{
			ARRAY_TYPE_ELEMENT tmp(*this->b_k);
			this->reg_op->mult_M(u_k,&tmp,true);
			if (this->prior.get())
				tmp -= *(this->p_M);
			shrink1(&tmp,REAL(1)/this->reg_op->get_weight(),this->d_k.get());
		}

		virtual void update_dk_bk(ARRAY_TYPE_ELEMENT* u_k)
		{
			this->reg_op->mult_M(u_k,this->b_k.get(),true);
			if (this->prior.get())
				*(this->b_k) -= *(this->p_M);
			shrink1(this->b_k.get(),REAL(1)/this->reg_op->get_weight(),this->d_k.get());
			*this->b_k -= *this->d_k;
		}
	};


	class sbL0RegularizationOperator : public sbRegularizationOperator
	{
	public:
		sbL0RegularizationOperator(boost::shared_ptr< linearOperator< ARRAY_TYPE_ELEMENT> > op,REAL _p = REAL(0.5)) : sbRegularizationOperator(op), p(_p) {}

		virtual void update_dk(ARRAY_TYPE_ELEMENT* u_k)
		{
			ARRAY_TYPE_ELEMENT tmp(*this->b_k);
			this->reg_op->mult_M(u_k,&tmp,true);
			if (this->prior.get())
				tmp -= *(this->p_M);
			pshrink(&tmp,REAL(1)/this->reg_op->get_weight(),p,this->d_k.get());
		}

		virtual void update_dk_bk(ARRAY_TYPE_ELEMENT* u_k)
		{
			this->reg_op->mult_M(u_k,this->b_k.get(),true);
			if (this->prior.get())
				*(this->b_k) -= *(this->p_M);
			pshrink(this->b_k.get(),REAL(1)/this->reg_op->get_weight(),p,this->d_k.get());
			*this->b_k -= *this->d_k;
		}
	protected:
		REAL p;
	};


	class sbL1GroupRegularizationOperator : public sbRegularizationOperator
	{
	public:
		sbL1GroupRegularizationOperator(std::vector<boost::shared_ptr< linearOperator<ARRAY_TYPE_ELEMENT> > > group)
	: sbRegularizationOperator()
	{
			op_cont = boost::shared_ptr<encodingOperatorContainer<ARRAY_TYPE_ELEMENT> >
			(new encodingOperatorContainer<ARRAY_TYPE_ELEMENT>);
			for (int i = 0; i < group.size(); i++)
				op_cont->add_operator(group[i]);
			reg_ops = group;
			this->reg_op = op_cont;
	}

		virtual void update_encoding_space(ARRAY_TYPE_ELEMENT* encoding_space)
		{
			for (int i=0; i < reg_ops.size(); i++){
				ARRAY_TYPE_ELEMENT tmp(*codom_dims,encoding_space->get_data_ptr()+op_cont->get_offset(i));
				tmp = *d_ks[i];
				tmp -= *b_ks[i];
				if (this->prior.get())
					tmp += *p_Ms[i];
			}
		}

		virtual void initialize(boost::shared_ptr< std::vector<size_t> > image_dims, REAL normalization_factor = REAL(1))
		{
			codom_dims = reg_ops.front()->get_codomain_dimensions();
			d_ks = std::vector< boost::shared_ptr<ARRAY_TYPE_ELEMENT> >(reg_ops.size());
			b_ks = std::vector< boost::shared_ptr<ARRAY_TYPE_ELEMENT> >(reg_ops.size());
			if (this->prior.get())
				p_Ms = std::vector< boost::shared_ptr<ARRAY_TYPE_ELEMENT> >(reg_ops.size());
			for (int i=0; i<reg_ops.size(); i++){
				d_ks[i] = boost::shared_ptr<ARRAY_TYPE_ELEMENT>(new ARRAY_TYPE_ELEMENT(*codom_dims));
				clear(d_ks[i].get());
				b_ks[i] = boost::shared_ptr<ARRAY_TYPE_ELEMENT>(new ARRAY_TYPE_ELEMENT(*codom_dims));
				clear(b_ks[i].get());
				if (this->prior.get()){
					p_Ms[i] = boost::shared_ptr<ARRAY_TYPE_ELEMENT>(new ARRAY_TYPE_ELEMENT(*codom_dims));
					reg_ops[i]->mult_M(this->prior.get(),p_Ms[i].get());
					*p_Ms[i] *= normalization_factor;
				}
			}
		}

		virtual void deinitialize()
		{
			d_ks.clear();
			b_ks.clear();
			p_Ms.clear();
		}

		virtual void update_dk(ARRAY_TYPE_ELEMENT* u_k)
		{
			ARRAY_TYPE_REAL s_k(*codom_dims);
			ARRAY_TYPE_ELEMENT *tmp = new ARRAY_TYPE_ELEMENT[reg_ops.size()];
			for (int i=0; i<reg_ops.size(); i++) {
				tmp[i] = *b_ks[i];
				this->reg_ops[i]->mult_M(u_k,&tmp[i],true);
				if (this->prior.get())
					tmp[i] -= *p_Ms[i];
				(i==0) ? s_k = *abs_square<ELEMENT_TYPE>(&tmp[i]) : s_k += *abs_square<ELEMENT_TYPE>(&tmp[i]);
			}
			sqrt_inplace(&s_k);
			for (int i=0; i<reg_ops.size(); i++) {
				shrinkd(&tmp[i],&s_k,REAL(1)/reg_ops[i]->get_weight(),d_ks[i].get());
			}
			delete[] tmp;
		}

		virtual void update_dk_bk(ARRAY_TYPE_ELEMENT* u_k)
		{
			ARRAY_TYPE_REAL s_k(*codom_dims);
			for (int i=0; i<reg_ops.size(); i++) {
				this->reg_ops[i]->mult_M(u_k,b_ks[i].get(),true);
				if (this->prior.get())
					*b_ks[i] -= *p_Ms[i];
				(i==0) ? s_k = *abs_square(b_ks[i].get()) : s_k += *abs_square(b_ks[i].get());
			}
			sqrt_inplace(&s_k);
			for (int i=0; i<reg_ops.size(); i++) {
				shrinkd(b_ks[i].get(),&s_k,REAL(1)/reg_ops[i]->get_weight(),d_ks[i].get());
				*b_ks[i] -= *d_ks[i];
			}
		}

		virtual boost::shared_ptr< std::vector<size_t> > get_codomain_dimensions(){
			return reg_ops.front()->get_codomain_dimensions();
		}

	protected:
		std::vector<boost::shared_ptr< linearOperator<ARRAY_TYPE_ELEMENT> > > reg_ops;
		std::vector< boost::shared_ptr<ARRAY_TYPE_ELEMENT> > d_ks;
		std::vector< boost::shared_ptr<ARRAY_TYPE_ELEMENT> > b_ks;
		std::vector< boost::shared_ptr<ARRAY_TYPE_ELEMENT> > p_Ms;
		boost::shared_ptr<encodingOperatorContainer<ARRAY_TYPE_ELEMENT> > op_cont;
		boost::shared_ptr< std::vector<size_t> > codom_dims;
	};


	class sbL0GroupRegularizationOperator : public sbRegularizationOperator
	{
	public:
		sbL0GroupRegularizationOperator(std::vector<boost::shared_ptr< linearOperator<ARRAY_TYPE_ELEMENT> > > group, REAL _p = REAL(0.5))
	: sbRegularizationOperator(), p(_p)
	{
			op_cont = boost::shared_ptr<encodingOperatorContainer<ARRAY_TYPE_ELEMENT> >
			(new encodingOperatorContainer<ARRAY_TYPE_ELEMENT>);
			for (int i = 0; i < group.size(); i++)
				op_cont->add_operator(group[i]);
			reg_ops = group;
			this->reg_op = op_cont;
	}

		virtual void update_encoding_space(ARRAY_TYPE_ELEMENT* encoding_space)
		{
			for (int i=0; i < reg_ops.size(); i++){
				ARRAY_TYPE_ELEMENT tmp(*codom_dims,encoding_space->get_data_ptr()+op_cont->get_offset(i));
				tmp = *d_ks[i];
				tmp -= *b_ks[i];
				if (this->prior.get())
					tmp += *p_Ms[i];
			}
		}

		virtual void initialize(boost::shared_ptr< std::vector<size_t> > image_dims, REAL normalization_factor = REAL(1))
		{
			codom_dims = reg_ops.front()->get_codomain_dimensions();
			d_ks = std::vector< boost::shared_ptr<ARRAY_TYPE_ELEMENT> >(reg_ops.size());
			b_ks = std::vector< boost::shared_ptr<ARRAY_TYPE_ELEMENT> >(reg_ops.size());
			if (this->prior.get())
				p_Ms = std::vector< boost::shared_ptr<ARRAY_TYPE_ELEMENT> >(reg_ops.size());
			for (int i=0; i<reg_ops.size(); i++){
				d_ks[i] = boost::shared_ptr<ARRAY_TYPE_ELEMENT>(new ARRAY_TYPE_ELEMENT(*codom_dims));
				clear(d_ks[i].get());
				b_ks[i] = boost::shared_ptr<ARRAY_TYPE_ELEMENT>(new ARRAY_TYPE_ELEMENT(*codom_dims));
				clear(b_ks[i].get());
				if (this->prior.get()){
					p_Ms[i] = boost::shared_ptr<ARRAY_TYPE_ELEMENT>(new ARRAY_TYPE_ELEMENT(*codom_dims));
					reg_ops[i]->mult_M(this->prior.get(),p_Ms[i].get());
					*p_Ms[i] *= normalization_factor;
				}
			}
		}

		virtual void deinitialize()
		{
			d_ks.clear();
			b_ks.clear();
			p_Ms.clear();
		}

		virtual void update_dk(ARRAY_TYPE_ELEMENT* u_k)
		{
			ARRAY_TYPE_REAL s_k(*codom_dims);
			ARRAY_TYPE_ELEMENT *tmp = new ARRAY_TYPE_ELEMENT[reg_ops.size()];
			for (int i=0; i<reg_ops.size(); i++) {
				tmp[i] = *b_ks[i];
				this->reg_ops[i]->mult_M(u_k,&tmp[i],true);
				if (this->prior.get())
					tmp[i] -= *p_Ms[i];
				(i==0) ? s_k = *abs_square<ELEMENT_TYPE>(&tmp[i]) : s_k += *abs_square<ELEMENT_TYPE>(&tmp[i]);
			}
			sqrt_inplace(&s_k);
			for (int i=0; i<reg_ops.size(); i++) {
				pshrinkd(&tmp[i],&s_k,REAL(1)/reg_ops[i]->get_weight(),p,d_ks[i].get());
			}
			delete[] tmp;
		}

		virtual void update_dk_bk(ARRAY_TYPE_ELEMENT* u_k)
		{
			ARRAY_TYPE_REAL s_k(*codom_dims);
			for (int i=0; i<reg_ops.size(); i++) {
				this->reg_ops[i]->mult_M(u_k,b_ks[i].get(),true);
				if (this->prior.get())
					*b_ks[i] -= *p_Ms[i];
				(i==0) ? s_k = *abs_square(b_ks[i].get()) : s_k += *abs_square(b_ks[i].get());
			}
			sqrt_inplace(&s_k);
			for (int i=0; i<reg_ops.size(); i++) {
				pshrinkd(b_ks[i].get(),&s_k,REAL(1)/reg_ops[i]->get_weight(),p,d_ks[i].get());
				*b_ks[i] -= *d_ks[i];
			}
		}

		virtual boost::shared_ptr< std::vector<size_t> > get_codomain_dimensions(){
			return reg_ops.front()->get_codomain_dimensions();
		}

	protected:
		std::vector<boost::shared_ptr< linearOperator<ARRAY_TYPE_ELEMENT> > > reg_ops;
		std::vector< boost::shared_ptr<ARRAY_TYPE_ELEMENT> > d_ks;
		std::vector< boost::shared_ptr<ARRAY_TYPE_ELEMENT> > b_ks;
		std::vector< boost::shared_ptr<ARRAY_TYPE_ELEMENT> > p_Ms;
		boost::shared_ptr<encodingOperatorContainer<ARRAY_TYPE_ELEMENT> > op_cont;
		boost::shared_ptr< std::vector<size_t> > codom_dims;
		REAL p;
	};

	class sbL2RegularizationOperator : public sbRegularizationOperator
	{
	public:

		sbL2RegularizationOperator(boost::shared_ptr< linearOperator<ARRAY_TYPE_ELEMENT> > op) : sbRegularizationOperator(op) {}

		virtual void update_dk(ARRAY_TYPE_ELEMENT* u_k)
		{
			*this->d_k = *this->b_k;
			this->reg_op->mult_M(u_k,this->d_k.get(),true);
			if (this->prior.get()){
				*this->d_k -= *this->p_M;
			}
			*(this->d_k) *= REAL(1)/(1+this->reg_op->get_weight());
		}

		virtual void update_dk_bk(ARRAY_TYPE_ELEMENT* u_k){
			update_dk(u_k);
			*(this->b_k) = *(this->d_k);
			*(this->b_k) *= this->reg_op->get_weight();
		}
	};

	class sbNonNegativityOperator : public sbRegularizationOperator
	{
	public:
		sbNonNegativityOperator(): sbRegularizationOperator(){
			this->reg_op = boost::shared_ptr<identityOperator<ARRAY_TYPE_ELEMENT> >(new identityOperator<ARRAY_TYPE_ELEMENT>);
		}

		virtual void initialize(boost::shared_ptr< std::vector<size_t> > image_dims,
				REAL normalization_factor = REAL(1))
		{
			this->reg_op->set_domain_dimensions(image_dims.get());
			this->reg_op->set_codomain_dimensions(image_dims.get());
			sbRegularizationOperator::initialize( image_dims, normalization_factor);
		}

		virtual void update_encoding_space(ARRAY_TYPE_ELEMENT* encoding_space){
			*encoding_space = *(this->d_k);
			clamp_min(encoding_space,REAL(0));
			*encoding_space += *(this->b_k);
		}

		virtual void update_dk(ARRAY_TYPE_ELEMENT* u_k){
			*(this->d_k) = *u_k;
			*(this->d_k) -= (*(this->b_k));
			clamp_min(this->d_k.get(),REAL(0));
		}

		virtual void update_dk_bk(ARRAY_TYPE_ELEMENT* u_k){
			update_dk(u_k);
			*(this->b_k) += *(this->d_k);
			*(this->b_k) -= *u_k;
		}
	};

	public:

	// Constructor
	//

	sbSolver() : linearOperatorSolver<ARRAY_TYPE_ELEMENT>()
	{
		normalization_mode_ = SB_NORMALIZE_TO_IMAGE_SPACE_IDENTITY;
		tolerance_ = REAL(0);
		outer_iterations_ = 10;
		inner_iterations_ = 1;
		num_reg_operators_ = 0;
		inner_solver_ = boost::shared_ptr<INNER_SOLVER>( new INNER_SOLVER() );
		non_negativity_filter_weight_ = REAL(0);
		use_x0_ = false;
	}

	// Destructor
	//

	virtual ~sbSolver() {}

	// Add regularization operator to group (for isotropic regularization)
	//

	virtual void add_regularization_group_operator( boost::shared_ptr< linearOperator<ARRAY_TYPE_ELEMENT> > op )
	{
		if( !op.get() ){
			throw std::runtime_error( "Error: sbSolver::add_regularization_group_operator : NULL operator provided" );
		}
		current_group_.push_back(op);
	}

	// Add isotroic regularization group (multiple groups allowed)
	//

	virtual void add_group(int L_norm=1)
	{
		if(current_group_.size()==0){
			throw std::runtime_error( "Error: sbSolver::add_group : no regularization group operators added" );
		}
		if (L_norm==2){
			for (int i=0; i<current_group_.size(); i++){
				regularization_operators_.push_back(boost::shared_ptr<sbL2RegularizationOperator>(new sbL2RegularizationOperator(current_group_[i])));
			}
		} else if (L_norm==0){
			boost::shared_ptr<sbL0GroupRegularizationOperator> group(new sbL0GroupRegularizationOperator(current_group_));
			regularization_operators_.push_back(group);
		}else if (L_norm ==1){
			boost::shared_ptr<sbL1GroupRegularizationOperator> group(new sbL1GroupRegularizationOperator(current_group_));
			regularization_operators_.push_back(group);
		} else throw std::runtime_error("Illega L-norm used in add_group");
		current_group_.clear();
	}

	virtual void add_group( boost::shared_ptr<ARRAY_TYPE_ELEMENT> prior, int L_norm=1 )
	{
		if(current_group_.size()==0){
			throw std::runtime_error( "Error: sbSolver::add_group : no regularization group operators added" );
		}
		if (L_norm==2){
			for (int i=0; i<current_group_.size(); i++){
				regularization_operators_.push_back(boost::shared_ptr<sbL2RegularizationOperator>(new sbL2RegularizationOperator(current_group_[i])));
				regularization_operators_.back()->set_prior(prior);
			}
		} else if (L_norm==0){
			boost::shared_ptr<sbL0GroupRegularizationOperator> group(new sbL0GroupRegularizationOperator(current_group_));
			group->set_prior(prior);
			regularization_operators_.push_back(group);
		} else if (L_norm==1){
			boost::shared_ptr<sbL1GroupRegularizationOperator> group(new sbL1GroupRegularizationOperator(current_group_));
			group->set_prior(prior);
			regularization_operators_.push_back(group);
		} else throw std::runtime_error("Illega L-norm used in add_group");
		current_group_.clear();
	}

	virtual void add_regularization_operator(boost::shared_ptr< linearOperator<ARRAY_TYPE_ELEMENT> > op, int L_norm=1 ){
		if (L_norm==1){
			regularization_operators_.push_back(boost::shared_ptr<sbL1RegularizationOperator>(new sbL1RegularizationOperator(op)));
		}else if (L_norm == 0){
			regularization_operators_.push_back(boost::shared_ptr<sbL0RegularizationOperator>(new sbL0RegularizationOperator(op)));
		}else{
			regularization_operators_.push_back(boost::shared_ptr<sbL2RegularizationOperator>(new sbL2RegularizationOperator(op)));
		}
	}

	virtual void add_regularization_operator(boost::shared_ptr< linearOperator<ARRAY_TYPE_ELEMENT> > op, boost::shared_ptr<ARRAY_TYPE_ELEMENT> prior, int L_norm=1 ){
		if (L_norm==1){
			regularization_operators_.push_back(boost::shared_ptr<sbL1RegularizationOperator>(new sbL1RegularizationOperator(op)));
			regularization_operators_.back()->set_prior(prior);
		}else if (L_norm == 0){
			regularization_operators_.push_back(boost::shared_ptr<sbL0RegularizationOperator>(new sbL0RegularizationOperator(op)));
			regularization_operators_.back()->set_prior(prior);
		}else{
			regularization_operators_.push_back(boost::shared_ptr<sbL2RegularizationOperator>(new sbL2RegularizationOperator(op)));
			regularization_operators_.back()->set_prior(prior);
		}
	}

	// Specify normalization mode.
	// The default mode is to use image space normalization.

	enum SB_normalization_mode{
		SB_NO_NORMALIZATION,
		SB_NORMALIZE_TO_IMAGE_SPACE_IDENTITY
	};

	virtual void set_normalization_mode( SB_normalization_mode mode ){
		normalization_mode_ = mode;
	}

	// Set termination criterium tolerance
	//

	virtual void set_tc_tolerance( REAL tolerance )
	{
		if( tolerance < REAL(0) )
			this->solver_warning( "Warning: sbSolver::set_tc_tolerence : tolerance cannot be negative. Ignored." );
		else tolerance_ = tolerance;
	}

	virtual void set_non_negativity_filter(REAL nnf){
		non_negativity_filter_weight_ = nnf;
	}

	// Set/get maximum number of outer Split-Bregman iterations
	//

	virtual void set_max_outer_iterations( unsigned int iterations ) { outer_iterations_ = iterations; }
	virtual unsigned int get_max_outer_iterations() { return outer_iterations_; }

	// Set/get maximum number of inner Split-Bregman iterations
	//

	virtual void set_max_inner_iterations( unsigned int iterations ) { inner_iterations_ = iterations; }
	virtual unsigned int get_max_inner_iterations() { return inner_iterations_; }

	virtual void set_use_inner_x0(bool use){ use_x0_=use; }

	// Get the inner solver
	//

	virtual boost::shared_ptr<INNER_SOLVER> get_inner_solver() { return inner_solver_; }

	// Provide the user an option to access u_k right after its update.
	//

	virtual bool post_linear_solver_callback( ARRAY_TYPE_ELEMENT* ) { return true; }

	//
	// Main solver interface
	//

	virtual boost::shared_ptr<ARRAY_TYPE_ELEMENT> solve( ARRAY_TYPE_ELEMENT *_f )
    		{
		// Check that operators etc. have been provided and consistent in dimensionality
		//
		validate_solver();

		// Define u_k
		//
		boost::shared_ptr<ARRAY_TYPE_ELEMENT> u_k( new ARRAY_TYPE_ELEMENT(*this->encoding_operator_->get_domain_dimensions()) );

		// Use x0 (if provided) as starting solution estimate
		//
		if( this->get_x0().get() )
			*u_k = *(this->get_x0());
		else
			clear(u_k.get());

		// Normalize and _then_ initialize (the order matters)
		boost::shared_ptr<ARRAY_TYPE_ELEMENT> f(new ARRAY_TYPE_ELEMENT(*_f));
		REAL normalization_factor = normalize_data( f.get() );
		initialize( normalization_factor );

		// Invoke the core solver
		//
		core( tolerance_, outer_iterations_, inner_iterations_, f, u_k);

		// Clean up memory occupied by the operator container and inner solver
		deinitialize();

		// Undo normalization
		*u_k /= normalization_factor;

		// ... and return the result
		//
		return u_k;
    		}

	protected:

	//
	// Everything beyond this point is internal to the implementation
	// and not intended to be exposed as a public interface
	//

	// Validate operator
	//

	virtual void validate_encoding_operator()
	{
		boost::shared_ptr< linearOperator<ARRAY_TYPE_ELEMENT> > op = this->get_encoding_operator();

		if( !op.get() ){
			throw std::runtime_error( "Error: sbSolver::validate_encoding_operator : operator not set" );
		}

		boost::shared_ptr< std::vector<size_t> > op_dims = op->get_domain_dimensions();
		if( op_dims->size() == 0 ){
			throw std::runtime_error( "Error: sbSolver::validate_encoding_operator : encoding operator must have specified domain dimensions" );
		}

		op_dims = op->get_codomain_dimensions();
		if( op_dims->size() == 0 ){
			throw std::runtime_error( "Error: sbSolver::validate_encoding_operator : encoding operator must have specified codomain dimensions" );
		}
	}

	// Validate regularization operator
	//

	virtual void validate_regularization_operators( std::vector<size_t> *image_dims )
	{
		if( image_dims->size() == 0 ){
			throw std::runtime_error( "Error: sbSolver::validate_regularization_operators : empty dimensions vector provided" );
		}

		for( unsigned int i=0; i<this->regularization_operators_.size(); i++ ){

			boost::shared_ptr< linearOperator<ARRAY_TYPE_ELEMENT> > op = regularization_operators_[i]->reg_op;
			boost::shared_ptr< std::vector<size_t> > op_dims = op->get_domain_dimensions();
			boost::shared_ptr< std::vector<size_t> > op_codims = op->get_codomain_dimensions();
			if (!op_codims.get()){
				throw std::runtime_error("Error: sbSolver::validate_regularization_operators : operator codomain dimension not set");
			}

			if( !op.get() ){
				throw std::runtime_error( "Error: sbSolver::validate_regularization_operators : invalid operator provided" );
			}

			if( *op_dims != *image_dims ){
				throw std::runtime_error( "Error: sbSolver::validate_regularization_operators : operator domain dimensions mismatch between encoding and regularization operators" );
			}
		}
	}

	// Check that the solver is set up properly
	virtual void validate_solver()
	{
		// Some tests to check if we are ready to go...
		//

		validate_encoding_operator();
		boost::shared_ptr< std::vector<size_t> > op_dims = this->encoding_operator_->get_domain_dimensions();
		validate_regularization_operators(op_dims.get());
	}

	// Initialize solver
	virtual void initialize( REAL normalization_factor = REAL(1) )
	{
		// Get image dimensions
		boost::shared_ptr< std::vector<size_t> > image_dims =
				this->encoding_operator_->get_domain_dimensions();

		if (non_negativity_filter_weight_ > REAL(0)){
			regularization_operators_.push_back(boost::shared_ptr<sbNonNegativityOperator>(new sbNonNegativityOperator));
			regularization_operators_.back()->set_weight(non_negativity_filter_weight_);
		}

		// Set up inner solver
		//

		enc_op_container_ = boost::shared_ptr<encodingOperatorContainer<ARRAY_TYPE_ELEMENT> >( new encodingOperatorContainer<ARRAY_TYPE_ELEMENT>() );
		inner_solver_->set_encoding_operator( enc_op_container_ );
		enc_op_container_->add_operator( this->encoding_operator_ );

		// Invoke initialization on all regularization operators
		//

		for (int i=0; i < regularization_operators_.size(); i++){
			regularization_operators_[i]->initialize(image_dims, normalization_factor);
			enc_op_container_->add_operator( regularization_operators_[i]->reg_op );
		}
	}

	// Clean up operator memory in the inner solver
	// Also restore the weights we temporarily changed

	virtual void deinitialize()
	{
		enc_op_container_ = boost::shared_ptr<encodingOperatorContainer<ARRAY_TYPE_ELEMENT> >( new encodingOperatorContainer<ARRAY_TYPE_ELEMENT>);
		inner_solver_->set_encoding_operator( enc_op_container_ );
		for (int i=0; i < regularization_operators_.size(); i++){
			regularization_operators_[i]->deinitialize();
		}
		if (non_negativity_filter_weight_ > REAL(0)){
			regularization_operators_.pop_back();
		}
	}

	// The core of the Split Bregman solver.
	//

	virtual void core( REAL tolerance, unsigned int outer_iterations, unsigned int inner_iterations,
			boost::shared_ptr<ARRAY_TYPE_ELEMENT> f,
			boost::shared_ptr<ARRAY_TYPE_ELEMENT> u_k )
	{
		// Image space dimensions
		boost::shared_ptr< std::vector<size_t> > image_dims =
				this->encoding_operator_->get_domain_dimensions();

		// Keep a copy of the "previous" u_k to compute the outer loop change of u_k
		//

		ARRAY_TYPE_ELEMENT u_k_prev;
		if( tolerance > REAL(0) || this->output_mode_ >= solver<ARRAY_TYPE_ELEMENT, ARRAY_TYPE_ELEMENT>::OUTPUT_VERBOSE ){
			u_k_prev = *u_k;
		}

		//
		// Outer loop
		//

		for( unsigned int outer_iteration=0; outer_iteration<outer_iterations; outer_iteration++ ) {

			if( this->output_mode_ >= solver<ARRAY_TYPE_ELEMENT, ARRAY_TYPE_ELEMENT>::OUTPUT_VERBOSE )
				GDEBUG_STREAM(std::endl << "SB outer loop iteration " << outer_iteration << std::endl << std::endl);

			//
			// Inner loop
			//

			for( unsigned int inner_iteration=0; inner_iteration<inner_iterations; inner_iteration++ ) {

				if( this->output_mode_ >= solver<ARRAY_TYPE_ELEMENT, ARRAY_TYPE_ELEMENT>::OUTPUT_VERBOSE )
					GDEBUG_STREAM(std::endl << "SB inner loop iteration " << inner_iteration << std::endl << std::endl);

				{ // Brackets used to free 'data' below as soon as it goes out of scope

					// Setup input vector to the encoding operator container (argument to the inner solver's solve)
					//

					ARRAY_TYPE_ELEMENT data(*enc_op_container_->get_codomain_dimensions());
					ARRAY_TYPE_ELEMENT tmp(f->dimensions(), data.get_data_ptr());

					tmp = *f;

					// Next add the regularization operators' data, d_k - b_k
					//

					for( unsigned int i=0; i< regularization_operators_.size(); i++ ){
						boost::shared_ptr<sbRegularizationOperator > op = regularization_operators_[i];
						tmp.create( *op->get_codomain_dimensions(), data.get_data_ptr()+enc_op_container_->get_offset(i+1) );
						op->update_encoding_space(&tmp);
					}

					// Solve for u_k
					//

					{
						if (use_x0_){
							get_inner_solver()->set_x0(u_k);
						}

						boost::shared_ptr<ARRAY_TYPE_ELEMENT> tmp_u_k =
								get_inner_solver()->solve( &data );

						// Invoke the post inner solver callback
						post_linear_solver_callback( tmp_u_k.get() );

						// Compute change in u_k
						if( this->output_mode_ >= solver<ARRAY_TYPE_ELEMENT, ARRAY_TYPE_ELEMENT>::OUTPUT_VERBOSE ){
							*u_k -= *tmp_u_k;
							GDEBUG_STREAM("u_k delta l2-norm (inner loop): " << nrm2(u_k.get()) << std::endl);
						}

						// Update u_k
						*u_k = *tmp_u_k;
					}
				}

				// Update d_k (and b_k in final inner iteration)
				//

				for( unsigned int i=0; i< regularization_operators_.size(); i++ ){
					boost::shared_ptr<sbRegularizationOperator > op = regularization_operators_[i];
					if( inner_iteration < inner_iterations-1 )
						op->update_dk(u_k.get());
					else
						op->update_dk_bk(u_k.get());
				}
			} // end of inner loop

			// Output change in u_k
			if( tolerance > REAL(0) || this->output_mode_ >= solver<ARRAY_TYPE_ELEMENT, ARRAY_TYPE_ELEMENT>::OUTPUT_VERBOSE ){
				u_k_prev *= ELEMENT_TYPE(-1);
				u_k_prev += *u_k;
				REAL delta = nrm2(&u_k_prev);

				if( this->output_mode_ >= solver<ARRAY_TYPE_ELEMENT, ARRAY_TYPE_ELEMENT>::OUTPUT_VERBOSE )
					GDEBUG_STREAM("u_k delta l2-norm (outer loop): " << delta << std::endl << std::endl);

				if( delta < tolerance )
					break;

				u_k_prev = *u_k;
			}
		} // end of outer loop
	}

	virtual REAL normalize_data( ARRAY_TYPE_ELEMENT *f )
	{
		REAL image_scale = REAL(1);

		if( normalization_mode_ == SB_NORMALIZE_TO_IMAGE_SPACE_IDENTITY ){

			//
			// Normalize to an average energy of "one intensity unit per image element"
			//

			boost::shared_ptr< linearOperator<ARRAY_TYPE_ELEMENT> > op = this->encoding_operator_;
			ARRAY_TYPE_ELEMENT tmp( *op->get_domain_dimensions() );
			op->mult_MH( f, &tmp );
			REAL sum = asum( &tmp );
			image_scale = REAL(tmp.get_number_of_elements())/sum;
			*f *= image_scale;
		}

		return image_scale;
	}

	protected:
	SB_normalization_mode normalization_mode_;
	REAL tolerance_;
	unsigned int outer_iterations_, inner_iterations_;
	unsigned int num_reg_operators_;
	std::vector< boost::shared_ptr<sbRegularizationOperator> > regularization_operators_;
	std::vector< boost::shared_ptr<linearOperator<ARRAY_TYPE_ELEMENT> > > current_group_;
	boost::shared_ptr<INNER_SOLVER> inner_solver_;
	boost::shared_ptr<encodingOperatorContainer<ARRAY_TYPE_ELEMENT> > enc_op_container_;
	std::vector<unsigned int> weights_backup_;
	REAL non_negativity_filter_weight_;
	bool use_x0_;
};
}
