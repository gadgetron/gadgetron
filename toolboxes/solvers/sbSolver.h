/*
	An implementation of the "Generalized Split Bregman Algorithm" - sec. 3.2. of the paper
	"The Split Bregman Method for L1-Regularized Problems" by Tom Goldstein and Stanley Osher. 
	Siam J. Imaging Sciences. Vol. 2, No. 2, pp. 323–343.
*/

#pragma once

#include "solver.h"
#include "matrixOperator.h"
#include "vector_td_utilities.h"

#include <vector>
#include <iostream>

template <class REAL, class ELEMENT_TYPE, class ARRAY_TYPE_REAL, class ARRAY_TYPE_ELEMENT> 
class sbSolver : public solver<ARRAY_TYPE_ELEMENT>
{
public:

	sbSolver( int output_mode = solver<ARRAY_TYPE_ELEMENT>::OUTPUT_SILENT ) : solver<ARRAY_TYPE_ELEMENT>( output_mode ) { 
		tolerance_ = get_zero<REAL>();
		outer_iterations_ = 10;
		inner_iterations_ = 3;
	}

	virtual ~sbSolver() {}

	virtual int set_inner_solver( boost::shared_ptr< solver<ARRAY_TYPE_ELEMENT> > solver ) {
		inner_solver_ = solver;
		return 0;
	}

	virtual int set_encoding_operator( boost::shared_ptr< matrixOperator<REAL, ARRAY_TYPE_ELEMENT> > op ) {
		encoding_operator_ = op;
		return 0;
	}

	virtual int add_regularization_operator( boost::shared_ptr< matrixOperator<REAL, ARRAY_TYPE_ELEMENT> > op ) {
		regularization_operators_.push_back(op);
		return 0;
	}

	virtual int add_regularization_group_operator( boost::shared_ptr< matrixOperator<REAL, ARRAY_TYPE_ELEMENT> > op ) {
		__regularization_group_operators_.push_back(op);
		return 0;
	}

	virtual int add_group() {
		regularization_group_operators_.push_back(__regularization_group_operators_);
		__regularization_group_operators_.clear();
		return 0;
	}

	virtual void set_tolerance( REAL tolerance ) {
		if( tolerance < get_zero<REAL>() ) this->solver_error( "sbSolver::set_tolerence : tolerance cannot be negative" );
		else tolerance_ = tolerance;
	}

	virtual void set_outer_iterations( unsigned int iterations ) {
		outer_iterations_ = iterations;
	}

	virtual void set_inner_iterations( unsigned int iterations ) {
		inner_iterations_ = iterations;
	}

	virtual void set_image_dimensions( boost::shared_ptr< std::vector<unsigned int> > dims ){
		image_dims_ = dims;
	}

	virtual bool solver_clear_real( ARRAY_TYPE_REAL* ) = 0;
	virtual bool solver_clear_element( ARRAY_TYPE_ELEMENT* ) = 0;
	virtual bool solver_sqrt( ARRAY_TYPE_REAL* ) = 0;
	virtual bool solver_scal( ELEMENT_TYPE, ARRAY_TYPE_ELEMENT* ) = 0;
	virtual bool solver_axpy_real( REAL, ARRAY_TYPE_REAL*, ARRAY_TYPE_REAL* ) = 0;
	virtual bool solver_axpy_element( ELEMENT_TYPE, ARRAY_TYPE_ELEMENT*, ARRAY_TYPE_ELEMENT* ) = 0;
	virtual REAL solver_asum( ARRAY_TYPE_ELEMENT* ) = 0;
	virtual boost::shared_ptr<ARRAY_TYPE_REAL> solver_norm( ARRAY_TYPE_ELEMENT* ) = 0;
	virtual boost::shared_ptr<ARRAY_TYPE_REAL> solver_norm_squared( ARRAY_TYPE_ELEMENT* ) = 0;
	virtual bool solver_shrink1( REAL, ARRAY_TYPE_ELEMENT*, ARRAY_TYPE_ELEMENT* ) = 0;
	virtual bool solver_shrinkd( REAL, ARRAY_TYPE_REAL*, ARRAY_TYPE_ELEMENT*, ARRAY_TYPE_ELEMENT* ) = 0;

	virtual boost::shared_ptr<ARRAY_TYPE_ELEMENT> solve( ARRAY_TYPE_ELEMENT *_f )
	{		
		// Some tests to see if we are ready to go...
		//
		if( !inner_solver_.get() ){
			this->solver_error( "sbSolver::solve : inner solver has not been set" );
			return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
		}

		if( !encoding_operator_.get() ){
			this->solver_error( "sbSolver::solve : encoding operator has not been set" );
			return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
		}

		if( regularization_operators_.size() == 0 && regularization_group_operators_.size() == 0 ){
			this->solver_error( "sbSolver::solve : at least one matrix regularizer must be added" );
			return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
		}

		if( !image_dims_.get() ){
			this->solver_error( "sbSolver::solve : image dimensions have not been set" );
			return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
		}

		// Make a copy of _f - we are going to "normalize" the input shortly
		//
		ARRAY_TYPE_ELEMENT f(*_f);
		if( !f.get_data_ptr() ){
			this->solver_error( "sbSolver::solve : memory allocation of f failed" );
			return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
		}

		// Initialize u_k to E^H f 
		//
		boost::shared_ptr<ARRAY_TYPE_ELEMENT> u_k = boost::shared_ptr<ARRAY_TYPE_ELEMENT>(new ARRAY_TYPE_ELEMENT());
		if( u_k.get() ) u_k->create( image_dims_.get() );

		if( !u_k.get() || !u_k->get_data_ptr() ){
			this->solver_error( "sbSolver::solve : memory allocation of u_k failed" );
			return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
		}    

		if( encoding_operator_->mult_MH( &f, u_k.get() ) < 0 ){
			this->solver_error( "sbSolver::solve : adjoint encoding operation failed on f" );
			return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
		}

		// Normalize u_k and f
		//
		ELEMENT_TYPE image_scale=get_one<ELEMENT_TYPE>();
		{
			// Normalize to an average energy of "one intensity unit per image element"
			REAL sum = solver_asum( u_k.get() );
			image_scale = mul<REAL>(( (REAL) (u_k->get_number_of_elements())/sum), get_one<ELEMENT_TYPE>() );
			solver_scal( image_scale, u_k.get() );
			solver_scal( image_scale, &f );
		}

		// Keep a copy of E^H f for subsequent rhs computations in the inner solver
		// 
		ARRAY_TYPE_ELEMENT EHf(*u_k.get());
		if( !EHf.get_data_ptr() ){
			this->solver_error( "sbSolver::solve : memory allocation of EHf failed" );
			return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
		}

		// Scale EHf with mu
		//
		if( !solver_scal( mul<REAL>(encoding_operator_->get_weight(), get_one<ELEMENT_TYPE>()), &EHf ) ){
			this->solver_error( "sbSolver::solve : error scaling EHf with mu" );
			return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
		}

		// Keep a copy of the "previous" u_k to compute the outer loop change of u_k
		// 
		ARRAY_TYPE_ELEMENT u_k_prev;
		if( tolerance_ > get_zero<REAL>() || this->output_mode_ >= solver<ARRAY_TYPE_ELEMENT>::OUTPUT_VERBOSE ){
			u_k_prev = *u_k.get();
			if( !u_k_prev.get_data_ptr() ){
				this->solver_error( "sbSolver::solve : memory allocation of u_k_prev failed" );
				return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
			}
		}

		// Initialize d_k and b_k arrays to zero images
		//
		//

		unsigned int num_reg_operators = regularization_operators_.size();
		for( unsigned int i=0; i<regularization_group_operators_.size(); i++ )
			num_reg_operators += regularization_group_operators_[i].size();

		boost::shared_array< boost::shared_ptr<ARRAY_TYPE_ELEMENT> > d_k = 
			boost::shared_array< boost::shared_ptr<ARRAY_TYPE_ELEMENT> >( new boost::shared_ptr<ARRAY_TYPE_ELEMENT>[num_reg_operators] );

		boost::shared_array< boost::shared_ptr<ARRAY_TYPE_ELEMENT> > b_k = 
			boost::shared_array< boost::shared_ptr<ARRAY_TYPE_ELEMENT> >( new boost::shared_ptr<ARRAY_TYPE_ELEMENT>[num_reg_operators] );

		if( !d_k.get() || !b_k.get() ){
			this->solver_error( "sbSolver::solve : memory allocation of d_k or b_k failed" );
			return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
		}

		for( unsigned int i=0; i<num_reg_operators; i++ ){

			d_k[i] = boost::shared_ptr<ARRAY_TYPE_ELEMENT>(new ARRAY_TYPE_ELEMENT());

			if( d_k[i] ) d_k[i]->create( image_dims_.get() );

			if( !d_k[i]->get_data_ptr() ){
				this->solver_error( "sbSolver::solve : memory allocation of d_k failed" );
				return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
			}

			if( !solver_clear_element( d_k[i].get() )){
				this->solver_error( "sbSolver::solve : failed to clear internal memory buffer d_k" );
				return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
			}

			b_k[i] = boost::shared_ptr<ARRAY_TYPE_ELEMENT>(new ARRAY_TYPE_ELEMENT());

			if( b_k[i] ) b_k[i]->create( image_dims_.get() );

			if( !b_k[i]->get_data_ptr() ){
				this->solver_error( "sbSolver::solve : memory allocation of b_k failed" );
				return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
			}

			if( !solver_clear_element( b_k[i].get() )){
				this->solver_error( "sbSolver::solve : failed to clear internal memory buffer b_k" );
				return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
			} 
		}

		// Outer loop
		//
		for( unsigned int outer_iteration=0; outer_iteration<outer_iterations_; outer_iteration++ ) {

			if( this->output_mode_ >= solver<ARRAY_TYPE_ELEMENT>::OUTPUT_MAX )
				std::cout << std::endl << "SB outer loop iteration " << outer_iteration << std::endl << std::endl;

			// Inner loop
			//
			for( unsigned int inner_iteration=0; inner_iteration<inner_iterations_; inner_iteration++ ) {

				if( this->output_mode_ >= solver<ARRAY_TYPE_ELEMENT>::OUTPUT_MAX )
					std::cout << std::endl << "SB inner loop iteration " << inner_iteration << std::endl << std::endl;

				// Form rhs for inner loop solver (initializes to adjoint encoding operator)
				// 
				ARRAY_TYPE_ELEMENT rhs(EHf);
				if( !rhs.get_data_ptr() ){
					this->solver_error( "sbSolver::solve : memory allocation of rhs failed" );
					return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
				}

				unsigned int operator_idx = 0;

				// Add regularization operators to rhs
				//
				for( unsigned int i=0; i<regularization_operators_.size(); i++ ){

					ARRAY_TYPE_ELEMENT tmp_diff, reg_out;
					if( tmp_diff.create( image_dims_.get() ) < 0 || reg_out.create( image_dims_.get() ) < 0 ){
						this->solver_error( "sbSolver::solve : memory allocation for regularization operator failed in rhs computation" );
						return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
					}    

					tmp_diff = *d_k[operator_idx];

					if( !solver_axpy_element( get_zero<ELEMENT_TYPE>()-get_one<ELEMENT_TYPE>(), b_k[operator_idx].get(), &tmp_diff )){
						this->solver_error( "sbSolver::solve : computation of regularization argument failed in rhs computation" );
						return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
					}    

					if( regularization_operators_.at(i)->mult_MH( &tmp_diff, &reg_out ) < 0 ){
						this->solver_error( "sbSolver::solve : application of regularization operator failed in rhs computation" );
						return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
					}    

					if( !solver_axpy_element( mul<REAL>(regularization_operators_.at(i)->get_weight(), get_one<ELEMENT_TYPE>()), &reg_out, &rhs )){
						this->solver_error( "sbSolver::solve : accumulation in rhs computation failed (1)" );
						return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
					}
					operator_idx++;
				}	

				// Add regularization group operators to rhs
				//
				for( unsigned int i=0; i<regularization_group_operators_.size(); i++ ){
					for( unsigned int j=0; j<regularization_group_operators_[i].size(); j++ ){

						ARRAY_TYPE_ELEMENT tmp_diff, reg_out;
						if( tmp_diff.create( image_dims_.get() ) < 0 || reg_out.create( image_dims_.get() ) < 0 ){
							this->solver_error( "sbSolver::solve : memory allocation for group regularization operator failed in rhs computation" );
							return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
						}    

						tmp_diff = *d_k[operator_idx];

						if( !solver_axpy_element( get_zero<ELEMENT_TYPE>()-get_one<ELEMENT_TYPE>(), b_k[operator_idx].get(), &tmp_diff )){
							this->solver_error( "sbSolver::solve : computation of group regularization argument failed in rhs computation" );
							return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
						}    

						if( regularization_group_operators_[i].at(j)->mult_MH( &tmp_diff, &reg_out ) < 0 ){
							this->solver_error( "sbSolver::solve : application of group regularization operator failed in rhs computation" );
							return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
						}    

						if( !solver_axpy_element( mul<REAL>(regularization_group_operators_[i].at(j)->get_weight(), get_one<ELEMENT_TYPE>()), &reg_out, &rhs )){
							this->solver_error( "sbSolver::solve : accumulation in rhs computation failed (2)" );
							return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
						}
						operator_idx++;
					}
				}	

				// Solve for u_k
				//
				{
					boost::shared_ptr<ARRAY_TYPE_ELEMENT> tmp = inner_solver_->solve(&rhs);

					// Compute change in u_k
					//
					if( this->output_mode_ >= solver<ARRAY_TYPE_ELEMENT>::OUTPUT_VERBOSE ){
						if( !solver_axpy_element( get_zero<ELEMENT_TYPE>()-get_one<ELEMENT_TYPE>(), tmp.get(), u_k.get() )){
							this->solver_error( "sbSolver::solve : error computing inner loop u_k delta" );
							return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
						}
						std::cout << std::endl << "u_k delta (inner loop): " << solver_asum(u_k.get()) << std::endl;
					}

					// Update u_k
					u_k = tmp;
				}

				operator_idx = 0;

				// Update d_k (and in the final inner iteration b_k also)
				//
				for( unsigned int i=0; i<regularization_operators_.size(); i++ ){

					ARRAY_TYPE_ELEMENT tmp_sum, reg_out;
					if( tmp_sum.create( image_dims_.get() ) < 0 || reg_out.create( image_dims_.get() ) < 0 ){
						this->solver_error( "sbSolver::solve : memory allocation for regularization operator failed in {d_k,b_k} update" );
						return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
					}

					tmp_sum = *b_k[operator_idx];

					if( regularization_operators_.at(i)->mult_M( u_k.get(), &reg_out ) < 0 ){
						this->solver_error( "sbSolver::solve : application of regularization operator failed in {d_k,b_k} update" );
						return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
					}

					if( !solver_axpy_element( get_one<ELEMENT_TYPE>(), &reg_out, &tmp_sum )){
						this->solver_error( "sbSolver::solve : computation of shrinkage_1 argument for d_k failed" );
						return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
					}

					// Update of d_k
					if( !solver_shrink1( reciprocal<REAL>(regularization_operators_.at(i)->get_weight()), &tmp_sum, d_k[operator_idx].get() )){
						this->solver_error( "sbSolver::solve : shrinkage_1 of d_k failed" );
						return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
					}

					// Update of b_k (only in the last inner iteration)
					if( inner_iteration == inner_iterations_-1 ){

						if( !solver_axpy_element( get_zero<ELEMENT_TYPE>()-get_one<ELEMENT_TYPE>(), d_k[operator_idx].get(), &reg_out )){
							this->solver_error( "sbSolver::solve : computation of update argument to b_k failed" );
							return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
						}

						// Update of b_k
						if( !solver_axpy_element( get_one<ELEMENT_TYPE>(), &reg_out, b_k[operator_idx].get() )){
							this->solver_error( "sbSolver::solve : update of b_k failed" );
							return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
						}
					}
					operator_idx++;
				}

				for( unsigned int i=0; i<regularization_group_operators_.size(); i++ ){

					unsigned int k = regularization_group_operators_[i].size();
					ARRAY_TYPE_ELEMENT *sums = new ARRAY_TYPE_ELEMENT[k], *reg_out = new ARRAY_TYPE_ELEMENT[k];

					if( !sums || !reg_out ){
						this->solver_error( "sbSolver::solve : host memory allocation for temporary arrays failed" );
						return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
					}

					ARRAY_TYPE_REAL s_k;
					if( s_k.create(image_dims_.get()) < 0 ){
						this->solver_error( "sbSolver::solve : memory allocation for s_k failed" );
						return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
					}

					if( !solver_clear_real(&s_k) ){
						this->solver_error( "sbSolver::solve : failed to clear s_k" );
						return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
					} 	  

					for( unsigned int j=0; j<k; j++ ){

						if( sums[j].create( image_dims_.get() ) < 0 || reg_out[j].create( image_dims_.get() ) < 0 ){
							this->solver_error( "sbSolver::solve : memory allocation for regularization operator failed in {d_k,b_k} update" );
							return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
						}

						sums[j] = *b_k[operator_idx+j];

						if( regularization_group_operators_[i].at(j)->mult_M( u_k.get(), &reg_out[j] ) < 0 ){
							this->solver_error( "sbSolver::solve : application of regularization operator failed in {d_k,b_k} update" );
							return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
						}

						if( !solver_axpy_element( get_one<ELEMENT_TYPE>(), &reg_out[j], &sums[j] )){
							this->solver_error( "sbSolver::solve : computation of shrinkage_d argument for d_k failed" );
							return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
						}

						boost::shared_ptr<ARRAY_TYPE_REAL> tmp_s_k = solver_norm_squared(&sums[j]);
						if( !solver_axpy_real( get_one<REAL>(), tmp_s_k.get(), &s_k )){
							this->solver_error( "sbSolver::solve : accumulation of s_k failed" );
							return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
						}
					}

					if( !solver_sqrt(&s_k) ){
						this->solver_error( "sbSolver::solve : sqrt of s_k failed" );
						return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
					}

					for( unsigned int j=0; j<k; j++ ){

						// Update of d_k
						if( !solver_shrinkd( reciprocal<REAL>(regularization_group_operators_[i].at(j)->get_weight()), &s_k, &sums[j], d_k[operator_idx+j].get() )){
							this->solver_error( "sbSolver::solve : shrinkage_d of d_k failed" );
							return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
						}

						// Update of b_k (only in the last inner iteration)
						if( inner_iteration == inner_iterations_-1 ){

							if( !solver_axpy_element( get_zero<ELEMENT_TYPE>()-get_one<ELEMENT_TYPE>(), d_k[operator_idx+j].get(), &reg_out[j] )){
								this->solver_error( "sbSolver::solve : computation of update argument to b_k failed" );
								return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
							}

							if( !solver_axpy_element( get_one<ELEMENT_TYPE>(), &reg_out[j], b_k[operator_idx+j].get() )){
								this->solver_error( "sbSolver::solve : update of b_k failed" );
								return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
							}
						}
					}
					operator_idx += k;
					delete[] sums; delete[] reg_out;
				}
			} // end of inner loop

			// Output change in u_k
			if( tolerance_ > get_zero<REAL>() || this->output_mode_ >= solver<ARRAY_TYPE_ELEMENT>::OUTPUT_VERBOSE ){

				if( !solver_scal( get_zero<ELEMENT_TYPE>()-get_one<ELEMENT_TYPE>(), &u_k_prev ) ){
					this->solver_error( "sbSolver::solve : error computing inner loop u_k delta (scale)" );
					return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
				}

				if( !solver_axpy_element( get_one<ELEMENT_TYPE>(), u_k.get(), &u_k_prev )){
					this->solver_error( "sbSolver::solve : error computing inner loop u_k delta (axpy)" );
					return boost::shared_ptr<ARRAY_TYPE_ELEMENT>();
				}
				
				REAL delta = solver_asum(&u_k_prev);
				
				if( this->output_mode_ >= solver<ARRAY_TYPE_ELEMENT>::OUTPUT_VERBOSE )
					std::cout << std::endl << "u_k delta (outer loop): " << delta << std::endl << std::endl;

				if( delta < tolerance_ )
					break;

				u_k_prev = *u_k.get();
			}

		} // end of outer loop

		// Undo intermediate scaling of u_k ...
		solver_scal( reciprocal<ELEMENT_TYPE>(image_scale), u_k.get() );

		// ... and return the result
		//

		return u_k;
	}

protected:
	REAL tolerance_;
	unsigned int outer_iterations_, inner_iterations_;
	boost::shared_ptr< std::vector<unsigned int> > image_dims_;
	boost::shared_ptr< solver<ARRAY_TYPE_ELEMENT> > inner_solver_;
	boost::shared_ptr< matrixOperator<REAL, ARRAY_TYPE_ELEMENT> > encoding_operator_;
	std::vector< boost::shared_ptr< matrixOperator<REAL, ARRAY_TYPE_ELEMENT> > > regularization_operators_;
	std::vector< boost::shared_ptr< matrixOperator<REAL, ARRAY_TYPE_ELEMENT> > > __regularization_group_operators_;
	std::vector< std::vector< boost::shared_ptr< matrixOperator<REAL, ARRAY_TYPE_ELEMENT> > > > regularization_group_operators_;
};
