#pragma once

#include "cgPreconditioner.h"
#include "complext.h"
#include "gpSolver.h"
#include "linearOperatorSolver.h"
#include "real_utilities.h"

#include <iostream>
#include <numeric>
#include <vector>

namespace Gadgetron {
/** Nonlinear conjugate gradient solver.
 * Adapted from Y.H. Dai & Y. Yuan 2001 "An Efficient Hybrid Conjugate Gradient Method for Unconstrained Optimization"
 * Annals of Operations Research, March 2001, Volume 103, Issue 1-4, pp 33-47
 *
 */

template <class ARRAY_TYPE> class nlcgSolver : public gpSolver<ARRAY_TYPE> {

  protected:
    typedef typename ARRAY_TYPE::element_type ELEMENT_TYPE;
    typedef typename realType<ELEMENT_TYPE>::Type REAL;
    typedef ARRAY_TYPE ARRAY_CLASS;
    typedef gpSolver<ARRAY_TYPE> GP;
    typedef typename gpSolver<ARRAY_TYPE>::l1GPRegularizationOperator l1GPRegularizationOperator;

  public:
    nlcgSolver() : gpSolver<ARRAY_TYPE>() {
        iterations_ = 10;
        tc_tolerance_ = (REAL)1e-7;
        non_negativity_constraint_ = false;
        dump_residual = false;
        threshold = REAL(1e-7);

        rho = 0.5f;
    }

    virtual ~nlcgSolver() {}

    virtual void set_rho(REAL _rho) { rho = _rho; }

    virtual boost::shared_ptr<ARRAY_TYPE> solve(ARRAY_TYPE* in) {
        if (this->encoding_operator_.get() == 0) {
            throw std::runtime_error("Error: nlcgSolver::compute_rhs : no encoding operator is set");
        }

        // Get image space dimensions from the encoding operator
        //

        std::vector<size_t> image_dims = this->encoding_operator_->get_domain_dimensions();
        if (image_dims.empty()) {
            throw std::runtime_error("Error: nlcgSolver::compute_rhs : encoding operator has not set domain dimension");
        }

        ARRAY_TYPE* x = new ARRAY_TYPE(image_dims); // The image. Will be returned inside a shared_ptr

        ARRAY_TYPE g(image_dims);     // Contains the gradient of the current step
        ARRAY_TYPE g_old(image_dims); // Contains the gradient of the previous step

        ARRAY_TYPE g_linear(image_dims); // Contains the linear part of the gradient;

        // If a prior image was given, use it for the initial guess.
        if (this->x0_.get()) {
            *x = *(this->x0_.get());
        } else {
            clear(x);
        }

        // Contains the encoding space of the linear regularization operators
        std::vector<ARRAY_TYPE> regEnc;

        // Initialize encoding space
        for (int i = 0; i < this->regularization_operators_.size(); i++) {
            regEnc.push_back(ARRAY_TYPE(this->regularization_operators_[i]->get_codomain_dimensions()));
            if (reg_priors[i].get()) {
                regEnc.back() = *reg_priors[i];
                regEnc.back() *= -std::sqrt(this->regularization_operators_[i]->get_weight());
            }
        }
        std::vector<ARRAY_TYPE> regEnc2 = regEnc;

        ARRAY_TYPE d(image_dims); // Search direction.
        clear(&d);

        ARRAY_TYPE encoding_space(
            in->get_dimensions()); // Contains the encoding space, or, equivalently, the residual vector

        ARRAY_TYPE g_step(image_dims); // Linear part of the gradient of the step d will be stored here

        ARRAY_TYPE encoding_space2(in->get_dimensions());
        REAL reg_res, data_res;

        if (this->output_mode_ >= solver<ARRAY_TYPE, ARRAY_TYPE>::OUTPUT_VERBOSE) {
            GDEBUG_STREAM("Iterating..." << std::endl);
        }
        REAL grad_norm0;

        for (int i = 0; i < iterations_; i++) {
            if (i == 0) {
                if (this->x0_.get()) {
                    this->encoding_operator_->mult_M(x, &encoding_space);

                } else
                    clear(&encoding_space);
                encoding_space -= *in;
                this->encoding_operator_->mult_MH(&encoding_space, &g_linear);

                g_linear *= this->encoding_operator_->get_weight();
                data_res =
                    std::sqrt(this->encoding_operator_->get_weight()) * real(dot(&encoding_space, &encoding_space));

                calc_regMultM(x, regEnc);
                for (int n = 0; n < regEnc.size(); n++)
                    if (reg_priors[n].get())
                        axpy(ELEMENT_TYPE(-std::sqrt(this->regularization_operators_[n]->get_weight())),
                             reg_priors[n].get(), &regEnc[n]);
                add_linear_gradient(regEnc, &g_linear);
                g = g_linear;
                this->add_gradient(x, &g);

                reg_res = REAL(0);

            } else {
                data_res = real(dot(&encoding_space, &encoding_space));
            }

            if (non_negativity_constraint_)
                solver_non_negativity_filter(x, &g);
            if (i == 0)
                grad_norm0 = nrm2(&g);
            REAL grad_norm = nrm2(&g);
            if (this->output_mode_ >= solver<ARRAY_TYPE, ARRAY_TYPE>::OUTPUT_VERBOSE) {

                GDEBUG_STREAM("Iteration " << i << ". Relative gradient norm: " << grad_norm / grad_norm0 << std::endl);
            }

            if (i == 0) {
                d -= g;
                if (this->precond_.get()) {
                    this->precond_->apply(&d, &d);
                    this->precond_->apply(&d, &d);
                }

            } else {

                g_step = g; // Not using g_step for anything right now, so let's use it for our beta calculation
                if (this->precond_.get()) {
                    this->precond_->apply(&g_step, &g_step); // Perform first half of the preconditioning
                    this->precond_->apply(&g_old, &g_old);
                }

                ELEMENT_TYPE g_old_norm = dot(&g_old, &g_old);
                ELEMENT_TYPE ggold = dot(&g_step, &g_old);
                g_old -= g_step;
                REAL gg = real(dot(&g_step, &g_step));
                ELEMENT_TYPE gy = -dot(&d, &g_old);
                // ELEMENT_TYPE beta = -dot(g,g_old)/g_old_norm; //PRP ste[
                // ELEMENT_TYPE theta = gy/g_old_norm;

                REAL betaDy = -gg / real(dot(&d, &g_old));
                REAL betaHS = real(dot(&g_step, &g_old)) / real(dot(&d, &g_old));
                REAL beta = std::max(REAL(0), std::min(betaDy, betaHS)); // Hybrid step size from Dai and Yuan 2001

                d *= beta;

                if (this->precond_.get())
                    this->precond_->apply(&g_step, &g_step); // Perform the rest of the preconditioning

                d -= g_step;
                GDEBUG_STREAM("Beta " << beta << std::endl);
            }

            this->encoding_operator_->mult_M(&d, &encoding_space2);

            calc_regMultM(&d, regEnc2);

            this->encoding_operator_->mult_MH(&encoding_space2, &g_step);
            g_step *= this->encoding_operator_->get_weight();

            add_linear_gradient(regEnc2, &g_step);

            REAL gd = real(dot(&g, &d));

            REAL alpha0 = REAL(1);

            // In the linear or semi-linear case, we can calculate the ideal step size.
            if (this->operators.size() == 0)
                alpha0 = -real(dot(&encoding_space, &encoding_space2) + calc_dot(regEnc, regEnc2)) /
                         real(dot(&encoding_space2, &encoding_space2) + calc_dot(regEnc2, regEnc2));

            REAL alpha;
            REAL old_norm = functionValue(&encoding_space, regEnc, x);

            g_old = g;

            {
                FunctionEstimator f(&encoding_space, &encoding_space2, &regEnc, &regEnc2, x, &d, &g_linear, &g_step,
                                    this);
                alpha = backtracking(f, alpha0, gd, rho, old_norm);
                // alpha=cg_linesearch(f,alpha0,gd,old_norm);
                if (alpha == 0) {
                    std::cerr << "Linesearch failed, returning current iteration" << std::endl;
                    return boost::shared_ptr<ARRAY_TYPE>(x);
                }
            }

            GDEBUG_STREAM("Alpha: " << alpha << std::endl);

            if (non_negativity_constraint_) {
                // Restore encoding space and gradient. Why not keep a copy? Memory!
                axpy(ELEMENT_TYPE(-alpha), &encoding_space2, &encoding_space);
                reg_axpy(-alpha, regEnc2, regEnc);
                axpy(ELEMENT_TYPE(-alpha), &g_step, &g_linear);

                ARRAY_TYPE x2 = *x;
                axpy(ELEMENT_TYPE(alpha), &d, &x2);

                clamp_min(&x2, REAL(0));

                d = x2;
                d -= *x;
                gd = real(dot(&g, &d));
                x2 = *x;
                alpha0 = 1;
                this->encoding_operator_->mult_M(&d, &encoding_space2);
                calc_regMultM(&d, regEnc2);

                this->encoding_operator_->mult_MH(&encoding_space2, &g_step);
                g_step *= this->encoding_operator_->get_weight();
                add_linear_gradient(regEnc2, &g_step);

                FunctionEstimator f(&encoding_space, &encoding_space2, &regEnc, &regEnc2, x, &d, &g_linear, &g_step,
                                    this);
                // alpha=gold(f,0,alpha0*1.5);
                // alpha = wolfesearch(f,alpha0,gd,rho,old_norm);
                alpha = backtracking(f, alpha0, gd, rho, old_norm);

                // alpha = cg_linesearch(f,alpha0,gd,old_norm);
                axpy(ELEMENT_TYPE(alpha), &d, x);
                if (alpha == 0) {
                    std::cerr << "Linesearch failed, returning current iteration" << std::endl;
                    return boost::shared_ptr<ARRAY_TYPE>(x);
                }
            } else {
                axpy(ELEMENT_TYPE(alpha), &d, x);
            }

            GDEBUG_STREAM("Function value: " << functionValue(&encoding_space, regEnc, x) << std::endl);

            g = g_linear;

            this->add_gradient(x, &g);

            iteration_callback(x, i, data_res, reg_res);

            if (grad_norm / grad_norm0 < tc_tolerance_)
                break;
        }

        return boost::shared_ptr<ARRAY_TYPE>(x);
    }

    // Set preconditioner
    //
    /*virtual void set_preconditioner( boost::shared_ptr< cgPreconditioner<ARRAY_TYPE> > precond ) {
  precond_ = precond;
  }*/

    // Set/get maximally allowed number of iterations
    //
    virtual void set_max_iterations(unsigned int iterations) { iterations_ = iterations; }
    virtual unsigned int get_max_iterations() { return iterations_; }

    // Set/get tolerance threshold for termination criterium
    //
    virtual void set_tc_tolerance(REAL tolerance) { tc_tolerance_ = tolerance; }
    virtual REAL get_tc_tolerance() { return tc_tolerance_; }

    virtual void set_non_negativity_constraint(bool non_negativity_constraint) {
        non_negativity_constraint_ = non_negativity_constraint;
    }

    virtual void set_dump_residual(bool dump_res) { dump_residual = dump_res; }
    // Set preconditioner
    //

    virtual void set_preconditioner(boost::shared_ptr<cgPreconditioner<ARRAY_TYPE>> precond) { precond_ = precond; }

    virtual void add_regularization_operator(boost::shared_ptr<linearOperator<ARRAY_TYPE>> op) {
        if (!op.get()) {
            throw std::runtime_error(
                "Error: linearOperatorSolver::add_regularization_operator : NULL operator provided");
        }
        this->regularization_operators_.push_back(op);
        reg_priors.push_back(boost::shared_ptr<ARRAY_TYPE>((ARRAY_TYPE*)0));
    }

    virtual void add_regularization_operator(boost::shared_ptr<linearOperator<ARRAY_TYPE>> op,
                                             boost::shared_ptr<ARRAY_TYPE> prior) {
        if (!op.get()) {
            throw std::runtime_error(
                "Error: linearOperatorSolver::add_regularization_operator : NULL operator provided");
        }

        this->regularization_operators_.push_back(op);
        reg_priors.push_back(prior);
    }

    virtual void add_regularization_operator(boost::shared_ptr<linearOperator<ARRAY_TYPE>> op, int L_norm) {
        if (L_norm == 1) {

            this->operators.push_back(
                boost::shared_ptr<l1GPRegularizationOperator>(new l1GPRegularizationOperator(op)));
        } else {
            add_regularization_operator(op);
        }
    }

    virtual void add_regularization_operator(boost::shared_ptr<linearOperator<ARRAY_TYPE>> op,
                                             boost::shared_ptr<ARRAY_TYPE> prior, int L_norm) {
        if (L_norm == 1) {
            this->operators.push_back(
                boost::shared_ptr<l1GPRegularizationOperator>(new l1GPRegularizationOperator(op, prior)));
        } else {
            add_regularization_operator(op, prior);
        }
    }

  protected:
    typedef typename std::vector<boost::shared_ptr<linearOperator<ARRAY_TYPE>>>::iterator csIterator;
    typedef typename std::vector<std::vector<boost::shared_ptr<linearOperator<ARRAY_TYPE>>>>::iterator csGroupIterator;

    virtual void iteration_callback(ARRAY_TYPE*, int i, REAL, REAL) {};

    ELEMENT_TYPE calc_dot(std::vector<ARRAY_TYPE>& x, std::vector<ARRAY_TYPE>& y) {
        ELEMENT_TYPE res(0);
        for (int i = 0; i < x.size(); i++)
            res += dot(&x[i], &y[i]);
        return res;
    }

    void add_linear_gradient(std::vector<ARRAY_TYPE>& elems, ARRAY_TYPE* g) {
        ARRAY_TYPE tmp(g->get_dimensions());
        for (int i = 0; i < elems.size(); i++) {
            this->regularization_operators_[i]->mult_MH(&elems[i], &tmp);
            axpy(ELEMENT_TYPE(std::sqrt(this->regularization_operators_[i]->get_weight())), &tmp, g);
        }
    }

    void calc_regMultM(ARRAY_TYPE* x, std::vector<ARRAY_TYPE>& elems) {
        for (int i = 0; i < elems.size(); i++) {
            this->regularization_operators_[i]->mult_M(x, &elems[i]);
            elems[i] *= std::sqrt(this->regularization_operators_[i]->get_weight());
        }
    }

    void reg_axpy(REAL alpha, std::vector<ARRAY_TYPE>& x, std::vector<ARRAY_TYPE>& y) {
        for (int i = 0; i < x.size(); i++) {
            axpy(ELEMENT_TYPE(alpha), &x[i], &y[i]);
        }
    }

    class FunctionEstimator {
      public:
        FunctionEstimator(ARRAY_TYPE* _encoding_space, ARRAY_TYPE* _encoding_step, std::vector<ARRAY_TYPE>* _regEnc,
                          std::vector<ARRAY_TYPE>* _regEnc_step, ARRAY_TYPE* _x, ARRAY_TYPE* _d, ARRAY_TYPE* _g,
                          ARRAY_TYPE* _g_step, nlcgSolver<ARRAY_TYPE>* _parent) {
            encoding_step = _encoding_step;
            encoding_space = _encoding_space;
            regEnc = _regEnc;
            regEnc_step = _regEnc_step;
            x = _x;
            xtmp = *x;
            d = _d;
            parent = _parent;
            alpha_old = 0;
            g = _g;
            g_step = _g_step;
        }

        REAL operator()(REAL alpha) {
            axpy(ELEMENT_TYPE(alpha - alpha_old), encoding_step, encoding_space);

            axpy(ELEMENT_TYPE(alpha - alpha_old), g_step, g);
            parent->reg_axpy(alpha - alpha_old, *regEnc_step, *regEnc);
            axpy(ELEMENT_TYPE(alpha - alpha_old), d, &xtmp);

            alpha_old = alpha;
            REAL res = parent->functionValue(encoding_space, *regEnc, &xtmp);
            return res;
        }

        ELEMENT_TYPE dir_deriv() {
            ARRAY_TYPE g_tmp = *g;
            parent->add_gradient(&xtmp, &g_tmp);
            return dot(d, &g_tmp);
        }

      private:
        REAL alpha_old;
        ARRAY_TYPE* encoding_step;
        ARRAY_TYPE* encoding_space;
        std::vector<ARRAY_TYPE>* regEnc;
        std::vector<ARRAY_TYPE>* regEnc_step;
        ARRAY_TYPE *x, *d;
        ARRAY_TYPE *g, *g_step;

        nlcgSolver<ARRAY_TYPE>* parent;
        ARRAY_TYPE xtmp;
    };
    friend class FunctionEstimator;

    /***
     * @brief Gold section search algorithm. Only works with unimodal functions, which we assume we're dealing with, at
     * least locally
     * @param f Functor to calculate the function to minimize
     * @param a Start of the bracketing
     * @param d End of bracketing
     * @return Value minimizing the function f.
     */
    REAL gold(FunctionEstimator& f, REAL a, REAL d) {
        const REAL gold = 1.0 / (1.0 + std::sqrt(5.0)) / 2;

        REAL b = d - (d - a) * gold;
        REAL c = (d - a) * gold - a;

        REAL fa = f(a);
        REAL fb = f(b);
        REAL fc = f(c);
        REAL fd = f(d);
        REAL tol = 1e-6;

        while (abs(a - d) > tol * (abs(b) + abs(c))) {
            if (fb > fc) {
                a = b;
                fa = fb;
                b = c;
                fb = fc;
                c = b * gold + (1.0 - gold) * d;
                fc = f(c);
            } else {
                d = c;
                fd = fc;
                c = b;
                fc = fb;
                b = c * gold + (1 - gold) * a;
                fb = f(b);
            }
        }
        if (fb < fc) {
            f(b);
            return b;
        } else {
            f(c);
            return c;
        }
    }

    /***
     * Armijo type linesearch
     * @param f
     * @param alpha0
     * @param gd
     * @param rho
     * @param old_norm
     * @return
     */
    REAL backtracking(FunctionEstimator& f, const REAL alpha0, const REAL gd, const REAL rho, const REAL old_norm) {
        REAL alpha;
        REAL delta = 0.1;
        REAL sigma = 0.9;
        // REAL precision = 0.0003; //Estimated precision of function evaluation
        REAL precision = 1e-4f; // Estimated precision of function evaluation
        bool wolfe = false;
        int k = 0;

        while (!wolfe) {
            alpha = alpha0 * std::pow(rho, k);
            // if (f(alpha) <= old_norm+alpha*delta*gd) wolfe = true;//Strong Wolfe condition..
            REAL fa = f(alpha);
            ELEMENT_TYPE dir_deriv = f.dir_deriv();
            if (((2 * delta - 1.0) * real(gd) >= real(dir_deriv)) && (fa < (old_norm + precision)))
                wolfe = true; // Approx Wolfe condition from Hager, W. and Zhang, H.SIAM Journal on Optimization 2005
                              // 16:1, 170-192
            if (abs(dir_deriv) > sigma * abs(gd))
                wolfe = false; // Strong Wolfe condition..
            k++;
            if (alpha == 0) {
                // GDEBUG_STREAM("Backtracking search failed, switching to slow wolfe-search" << std::endl);
                // return wolfesearch(f,alpha0,gd,rho,old_norm);
                return 0;
            }
        }

        return alpha;
    }

    /***
     * Line search taken from Numerical Optimization (Wright and Nocedal 1999).
     * Adapted from the scipy optimize algorithm.
     * Like the gold-section method it works quite poorly in practice.
     * @param f
     * @param alpha0
     * @param gd
     * @param rho
     * @param old_norm
     * @return
     */
    REAL wolfesearch(FunctionEstimator& f, const REAL alpha_init, const REAL gd, const REAL rho, const REAL old_norm) {
        using std::abs;
        using std::sqrt;
        REAL delta = 0.01;
        unsigned int k = 0;
        REAL alpha0 = alpha_init;
        REAL f0 = f(alpha0);

        if (f0 <= old_norm + alpha0 * delta * gd) { // Strong Wolfe condition..
            return alpha0;
        }

        REAL alpha1 = -gd * alpha0 * alpha0 / 2.0 / (f0 - old_norm - gd * alpha0);
        // GDEBUG_STREAM("F0 " <<f0 << " old " << old_norm << " gd " << gd <<std::endl);
        GDEBUG_STREAM("Alpha0: " << alpha0 << std::endl);
        // GDEBUG_STREAM("Alpha1: "  << alpha1 << std::endl);
        REAL f1 = f(alpha1);

        if (f1 <= old_norm + alpha1 * delta * gd) { // Strong Wolfe condition..
            return alpha1;
        }

        while (alpha1 > 0) {
            double factor = alpha0 * alpha0 * alpha1 * alpha1 * (alpha1 - alpha0);
            double a =
                alpha0 * alpha0 * (f1 - old_norm - gd * alpha1) - alpha1 * alpha1 * (f0 - old_norm - gd * alpha0);
            a /= factor;

            double b = -alpha0 * alpha0 * alpha0 * (f1 - old_norm - gd * alpha1) +
                       alpha1 * alpha1 * alpha1 * (f0 - old_norm - gd * alpha0);
            b /= factor;

            double alpha2 = (-b + std::sqrt(std::abs(b * b - 3 * a * gd))) / (3 * a);
            REAL f2 = f(alpha2);
            // GDEBUG_STREAM("a " << a << "b " << b << std::endl);
            GDEBUG_STREAM("Alpha1: " << alpha1 << std::endl);
            GDEBUG_STREAM("Alpha2: " << alpha2 << std::endl);
            if (f2 < old_norm + alpha2 * delta * gd) { // Strong Wolfe condition..
                return alpha2;
            }

            if (((alpha1 - alpha2) > (alpha1 / 2.0)) || ((1.0 - alpha2 / alpha1) < 0.96)) {
                alpha2 = alpha1 / 2.0;
            }

            alpha0 = alpha1;
            alpha1 = alpha2;
            f0 = f1;
            f1 = f2;
            k++;
        }

        throw std::runtime_error("Wolfe line search failed");
    }

    /***
     * CG linesearch adapted from  Hager, W. and Zhang, H.SIAM Journal on Optimization 2005 16:1, 170-192
     * @param f
     * @param alpha0
     * @param gd
     * @param rho
     * @param old_norm
     * @return
     */
    REAL cg_linesearch(FunctionEstimator& f, const REAL alpha0, const REAL gd, const REAL old_norm) {
        REAL delta = 0.1;
        REAL sigma = 0.9;
        REAL nabla = 0.66;
        // REAL precision = 0.0003; //Estimated precision of function evaluation
        REAL precision = 1e-4f; // Estimated precision of function evaluation

        REAL a = 0;
        REAL b = alpha0;

        REAL ak = a;
        REAL bk = b;
        REAL fa = old_norm;
        ELEMENT_TYPE a_deriv = gd;
        REAL fb = f(alpha0);
        ELEMENT_TYPE b_deriv = f.dir_deriv();

        while (abs(a - b) > 0) {
            if ((((2 * delta - 1.0) * real(gd) >= real(b_deriv)) &&
                 (fb < old_norm + precision)) && // Check Approximate Wolfe conditions
                (abs(b_deriv) <= sigma * abs(gd))) {
                f(b);
                return b;
            }

            if ((((2 * delta - 1.0) * real(gd) >= real(a_deriv)) &&
                 (fa < old_norm + precision)) && // Check Approximate Wolfe conditions
                (abs(a_deriv) <= sigma * abs(gd))) {
                f(a);
                return a;
            }

            secant2(a, b, f, old_norm + precision);
            if ((b - a) > nabla * (bk - ak)) {
                REAL c = (a + b) / 2;
                interval_update(a, b, c, f, old_norm);
            }
            if (a != ak) {
                fa = f(a);
                a_deriv = f.dir_deriv();
            }

            if (b != bk) {
                fb = f(b);
                b_deriv = f.dir_deriv();
            }

            ak = a;
            bk = b;

            GDEBUG_STREAM("a: " << a << " b: " << b << std::endl);
        }
        return 0;
        // throw std::runtime_error("CG_linesearch failed");
    }

    void secant2(REAL& a, REAL& b, FunctionEstimator& f, REAL old_norm) {
        REAL fa = f(a);
        ELEMENT_TYPE dfa = f.dir_deriv();
        REAL fb = f(b);
        ELEMENT_TYPE dfb = f.dir_deriv();

        REAL c = real((a * dfb - b * dfa) / (dfb - dfa));

        REAL fc = f(c);
        ELEMENT_TYPE dfc = f.dir_deriv();

        REAL A = a;
        REAL B = b;

        interval_update(A, B, c, f, old_norm);

        if (c == B) {
            c = real((b * dfc - c * dfb) / (dfc - dfb));
            interval_update(A, B, c, f, old_norm);
        }
        if (c == A) {
            c = real((a * dfc - c * dfa) / (dfc - dfa));
            interval_update(A, B, c, f, old_norm);
        }

        a = A;
        b = B;
    }

    void interval_update(REAL& a, REAL& b, REAL c, FunctionEstimator& f, REAL old_norm) {
        REAL theta = 0.5;
        if (c < a || c > b)
            return; // C not in interval
        REAL fc = f(c);
        ELEMENT_TYPE dfc = f.dir_deriv();

        if (real(dfc) >= 0) {
            b = c;
            return;
        }
        if (fc < old_norm) {
            a = c;
            return;
        }
        b = c;
        while (true) {
            REAL d = (1 - theta) * a + theta * b;
            REAL fd = f(d);
            ELEMENT_TYPE dfd = f.dir_deriv();

            if (real(dfd) >= 0) {
                b = d;
                return;
            }
            if (fd < old_norm) {
                a = d;
            } else
                b = d;

            GDEBUG_STREAM("Interval a: " << a << " b: " << b << std::endl);
        }
    }

    REAL functionValue(ARRAY_TYPE* encoding_space, std::vector<ARRAY_TYPE>& regEnc, ARRAY_TYPE* x) {
        REAL res = this->encoding_operator_->get_weight() * abs(dot(encoding_space, encoding_space));

        for (int i = 0; i < this->operators.size(); i++) {
            res += this->operators[i]->magnitude(x);
        }

        res += abs(calc_dot(regEnc, regEnc));
        return res;
    }

  protected:
    // Preconditioner
    // boost::shared_ptr< cgPreconditioner<ARRAY_TYPE> > precond_;
    // Maximum number of iterations
    unsigned int iterations_;
    bool non_negativity_constraint_;
    REAL tc_tolerance_;
    REAL threshold;
    bool dump_residual;
    REAL rho;

    // Preconditioner

    std::vector<boost::shared_ptr<ARRAY_TYPE>> reg_priors;
    boost::shared_ptr<cgPreconditioner<ARRAY_TYPE>> precond_;
};
} // namespace Gadgetron
