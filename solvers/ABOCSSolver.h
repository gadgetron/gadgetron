#pragma once

namespace Gadgetron{
/**
 *
* Accelerated barrier optimization compressed sensing (ABOCS) reconstruction for cone-beam CT: Phantom studies
* Tianye Niu and Lei Zhu, Med. Phys. 39, 4588 (2012), DOI:10.1118/1.4729837
* Note we use a slightly different solver
*/
//Curiously recurring template pattern....
template<class GPBB> class ABOCSSolver : public GPBB{

	protected:
		typedef typename GPBB::ARRAY_CLASS ARRAY;
		typedef typename ARRAY::element_type ELEMENT_TYPE;
		typedef typename realType<ELEMENT_TYPE>::Type REAL;
		REAL beta; // Basically the maximal height of the Log Barrier
		REAL eps; // Data residual / position of the log barrier



	public:

		ABOCSSolver() : GPBB(), beta(REAL(1e-6)), eps(REAL(1e10)){

		}

		virtual void set_beta(REAL _beta){ beta = _beta;}

		virtual void set_eps(REAL _eps){ eps = _eps;}

		virtual ~ABOCSSolver(){};
		virtual boost::shared_ptr<ARRAY> solve(ARRAY* in)
			{
				if( this->encoding_operator_.get() == 0 ){
					BOOST_THROW_EXCEPTION(runtime_error("Error: gpBBSolver::compute_rhs : no encoding operator is set" ));
				}

				// Get image space dimensions from the encoding operator
				//

				boost::shared_ptr< std::vector<unsigned int> > image_dims = this->encoding_operator_->get_domain_dimensions();
				if( image_dims->size() == 0 ){
					BOOST_THROW_EXCEPTION(runtime_error("Error: gpBBSolver::compute_rhs : encoding operator has not set domain dimension" ));
				}

				ARRAY * x = new ARRAY;
				x->create(image_dims.get());

				ARRAY x_old(image_dims.get());

				ARRAY * g = new ARRAY;
				g->create(image_dims.get());
				ARRAY *  g_old = new ARRAY;
				g_old->create(image_dims.get());

				if (this->x0_.get()){
					*x = *(this->x0_.get());
				} else  {
					clear(x);
				}

				ARRAY encoding_space;
				REAL reg_res,data_res;
				encoding_space.create(in->get_dimensions().get());
				if( this->output_mode_ >= solver<ARRAY,ARRAY>::OUTPUT_VERBOSE ){
					std::cout << "Iterating..." << std::endl;
				}
				for (int i = 0; i < this->iterations_; i++){
					if ((i==0) && (!this->x0_.get())){
						clear(g);

						this->encoding_operator_->mult_MH(in,g);
						if (this->precond_.get()) {
							this->precond_->apply(g,g);
							this->precond_->apply(g,g);
						}
						data_res = real(dot(in,in));
						REAL abocs = std::max(beta,eps-data_res);
						*g *=  -this->encoding_operator_->get_weight()/abocs;

						reg_res=REAL(0);
					} else {
						this->encoding_operator_->mult_M(x,&encoding_space);
						axpy(REAL(-1),in,&encoding_space);
						data_res = real(dot(&encoding_space,&encoding_space));
						this->encoding_operator_->mult_MH(&encoding_space,g);
						if (this->precond_.get()) {
							this->precond_->apply(g,g);
							this->precond_->apply(g,g);
						}
						REAL abocs = std::max(beta,eps-data_res);
						*g *=  this->encoding_operator_->get_weight()/abocs;
					}

					this->add_gradient(x,g); // Adds the gradient from all the regularization operators

					if( this->output_mode_ >= solver<ARRAY,ARRAY>::OUTPUT_VERBOSE ){
						std::cout << "Data residual: " << data_res << std::endl;
					}

					if (this->non_negativity_constraint_) this->solver_non_negativity_filter(x,g);
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
						x_old -= *x;
						*g_old -= *g;
						ELEMENT_TYPE xx = dot(&x_old,&x_old);
						ELEMENT_TYPE gx = dot(g_old,&x_old);

						ELEMENT_TYPE nabla1 = xx/gx;

						/* This is the code that enables the adaptive step size.
							 REAL gg = dot(g_old,&x_old);
							 REAL nabla2 = gx/gg;
							 if ((nabla2/nabla1) < 0.5) nabla = nabla2;
							 else nabla = nabla1;*/
						nabla = nabla1;
					}

					ARRAY * tmp;
					tmp=g_old;
					g_old=g;
					g=tmp;

					x_old = *x;
					REAL grad_norm = nrm2(g_old);

					if( this->output_mode_ >= solver<ARRAY,ARRAY>::OUTPUT_VERBOSE ){
						std::cout << "Iteration " <<i << ". Gradient norm: " <<  grad_norm << std::endl;
					}
					this->iteration_callback(x,i,data_res,reg_res);
					axpy(-nabla,g_old,x);
					if (this->non_negativity_constraint_) clamp_min(x,REAL(0));
					if (grad_norm < this->tc_tolerance_)  break;
				}
				delete g,g_old;

				return boost::shared_ptr<ARRAY>(x);
			}




};


}
