#pragma once

#include "gpSolver.h"



namespace Gadgetron {



template<class ARRAY> class FISTASolver : public gpSolver<ARRAY> {

	typedef typename ARRAY::element_type ELEMENT_TYPE;
	typedef typename realType<ELEMENT_TYPE>::Type REAL;


public:
	virtual boost::shared_ptr<ARRAY> solve(ARRAY* in){

			if( this->encoding_operator_.get() == 0 ){
				BOOST_THROW_EXCEPTION(runtime_error("Error: gpBBSolver::compute_rhs : no encoding operator is set" ));
			}

			// Get image space dimensions from the encoding operator
			//

			boost::shared_ptr< std::vector<unsigned int> > image_dims = this->encoding_operator_->get_domain_dimensions();
			if( image_dims->size() == 0 ){
				BOOST_THROW_EXCEPTION(runtime_error("Error: gpBBSolver::compute_rhs : encoding operator has not set domain dimension" ));
			}

			ARRAY * y = new ARRAY;
			y->create(image_dims.get());

			ARRAY y_old(image_dims.get());



			ARRAY * g = new ARRAY;
			g->create(image_dims.get());
			ARRAY *  g_old = new ARRAY;
			g_old->create(image_dims.get());

			if (this->x0_.get()){
				*y = *(this->x0_.get());
			} else  {
				clear(y);
			}

			//ARRAY x_old(*y);
			//ARRAY x(*y);


			ARRAY encoding_space;
			REAL reg_res,data_res;
			encoding_space.create(in->get_dimensions().get());
			if( this->output_mode_ >= solver<ARRAY,ARRAY>::OUTPUT_VERBOSE ){
				std::cout << "Iterating..." << std::endl;
			}

			REAL t = 1;

			for (int i = 0; i < iterations_; i++){


				if ((i==0) && (!this->x0_.get())){
					clear(g);

					this->encoding_operator_->mult_MH(in,g);
					if (precond_.get()) {
						precond_->apply(g,g);
						precond_->apply(g,g);
					}

					*g *=  -this->encoding_operator_->get_weight();
					data_res = real(dot(in,in));
					reg_res=REAL(0);
				} else {
					this->encoding_operator_->mult_M(y,&encoding_space);
					axpy(REAL(-1),in,&encoding_space);
					data_res = real(dot(&encoding_space,&encoding_space));
					this->encoding_operator_->mult_MH(&encoding_space,g);
					if (precond_.get()) {
						precond_->apply(g,g);
						precond_->apply(g,g);
					}
					*g *=  this->encoding_operator_->get_weight();
				}

				this->add_gradient(y,g); // Adds the gradient from all the regularization operators

				if( this->output_mode_ >= solver<ARRAY,ARRAY>::OUTPUT_VERBOSE ){
					std::cout << "Data residual: " << data_res << std::endl;
				}

				//if (non_negativity_constraint_) solver_non_negativity_filter(y,g);
				ELEMENT_TYPE nabla;
				if (i==0){
					ARRAY tmp_encoding = *in;
					this->encoding_operator_->mult_M(g,&tmp_encoding);
					if (this->x0_.get()){
						nabla = dot(&encoding_space,&tmp_encoding)/dot(&tmp_encoding,&tmp_encoding);
					} else {
						nabla = -dot(in,&tmp_encoding)/dot(&tmp_encoding,&tmp_encoding);
					}
				} else {
					y_old -= *y;
					*g_old -= *g;
					ELEMENT_TYPE xx = dot(&y_old,&y_old);
					ELEMENT_TYPE gx = dot(g_old,&y_old);

					ELEMENT_TYPE nabla1 = xx/gx;

					/* This is the code that enables the adaptive step size.
				 REAL gg = dot(g_old,&y_old);
				 REAL nabla2 = gx/gg;
				 if ((nabla2/nabla1) < 0.5) nabla = nabla2;
				 else nabla = nabla1;*/
					nabla = nabla1;
					//nabla = REAL(1e-4);
				}

				y_old = *y;
				axpy(-nabla,g,y);

				prox(y);
/*
				REAL tnew = (REAL(1)+std::sqrt(REAL(1)+REAL(4)*t*t))/REAL(2);
				x = *y;

				*y *= (t-REAL(1))/tnew;

				axpy(-(t-REAL(1)/tnew),&x_old,y);
				t = tnew;

				x_old = x;
				*/
				//iteration_callback(y,i,data_res,reg_res);



				REAL grad_norm = nrm2(g);

				if( this->output_mode_ >= solver<ARRAY,ARRAY>::OUTPUT_VERBOSE ){
					std::cout << "Iteration " <<i << ". Gradient norm: " <<  grad_norm << std::endl;
				}

				//if (non_negativity_constraint_) clamp_min(x,REAL(0));
				if (grad_norm < tc_tolerance_)  break;

				*g_old = *g;
			  /**g_old = y_old;
				*g_old -= *y;
				*g_old *= REAL(1)/nabla;*/

			}
			delete g,g_old;

			return boost::shared_ptr<ARRAY>(y);
		}

	void set_frame(boost::shared_ptr<linearOperator<ARRAY> > frame){
		tightFrame = frame;
	}


	virtual void set_max_iterations( unsigned int iterations ) { iterations_ = iterations; }
	virtual unsigned int get_max_iterations() { return iterations_; }

	// Set/get tolerance threshold for termination criterium
	//
	virtual void set_tc_tolerance( REAL tolerance ) { tc_tolerance_ = tolerance; }
	virtual REAL get_tc_tolerance() { return tc_tolerance_; }


	virtual void set_dump_residual(bool dump_res){
		dump_residual = dump_res;
	}
	// Set preconditioner
	//

	virtual void set_preconditioner( boost::shared_ptr< cgPreconditioner<ARRAY> > precond ) {
		precond_ = precond;
	}
protected:
	void prox(ARRAY* x){
		if (tightFrame.get()){
			ARRAY tmp(tightFrame->get_codomain_dimensions());
			clear(&tmp);
			tightFrame->mult_M(x,&tmp,false);
			shrink1(&tmp,tightFrame->get_weight());
			tightFrame->mult_MH(&tmp,x,false);
		}

	}

	unsigned int iterations_;

	REAL tc_tolerance_;
	REAL threshold;
	bool dump_residual;
	// Preconditioner
	boost::shared_ptr< cgPreconditioner<ARRAY> > precond_;

	boost::shared_ptr<linearOperator<ARRAY> > tightFrame; //This has to be a a tightframe (or basis), meaning M^H M x = x




};
}
