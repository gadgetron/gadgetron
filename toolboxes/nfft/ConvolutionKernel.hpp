
#include "ConvolutionKernel.h"

#include <cmath>

namespace Gadgetron
{
    __inline__ __device__ __host__
    double bessi0(double x)
    {
        double denominator;
        double numerator;
        double z;

        if (x == 0.0) {
            return 1.0;
        } else {
            z = x * x;
            numerator = (z * (z * (z * (z * (z * (z * (z * (z * (z * (z * (z * (z * (z *
                                                                                     (z * 0.210580722890567e-22 +
                                                                                      0.380715242345326e-19) +
                                                                                     0.479440257548300e-16) +
                                                                                0.435125971262668e-13) +
                                                                           0.300931127112960e-10) +
                                                                      0.160224679395361e-7) +
                                                                 0.654858370096785e-5) + 0.202591084143397e-2) +
                                                       0.463076284721000e0) + 0.754337328948189e2) +
                                             0.830792541809429e4) + 0.571661130563785e6) +
                                   0.216415572361227e8) + 0.356644482244025e9) +
                         0.144048298227235e10);

            denominator = (z * (z * (z - 0.307646912682801e4) +
                                0.347626332405882e7) - 0.144048298227235e10);
        }

        return -numerator / denominator;
    }

    __inline__ __device__ __host__
    float bessi0(float x)
    {
        float denominator;
        float numerator;
        float z;

        if (x == 0.0f) {
            return 1.0f;
        } else {
            z = x * x;
            numerator = (z * (z * (z * (z * (z * (z * (z * (z * (z * (z * (z * (z * (z *
                                                                                     (z * 0.210580722890567e-22f +
                                                                                      0.380715242345326e-19f) +
                                                                                     0.479440257548300e-16f) +
                                                                                0.435125971262668e-13f) +
                                                                           0.300931127112960e-10f) +
                                                                      0.160224679395361e-7f) +
                                                                 0.654858370096785e-5f) + 0.202591084143397e-2f) +
                                                       0.463076284721000e0f) + 0.754337328948189e2f) +
                                             0.830792541809429e4f) + 0.571661130563785e6f) +
                                   0.216415572361227e8f) + 0.356644482244025e9f) +
                         0.144048298227235e10f);

            denominator = (z * (z * (z - 0.307646912682801e4f) +
                                0.347626332405882e7f) - 0.144048298227235e10f);
        }

        return -numerator / denominator;
    }


    // Kaiser Bessel according to Beatty et. al. IEEE TMI 2005;24(6):799-808.
    // There is a slight difference wrt Jackson's formulation, IEEE TMI 1991;10(3):473-478.
    __inline__ __device__ __host__
    double KaiserBessel(double u, double matrix_size_os,
                        double one_over_W, double beta)
    {
        double _tmp = 2.0 * u * one_over_W;
        double tmp = _tmp * _tmp;
        double arg = beta * std::sqrt(1.0 - tmp);
        double bessi = bessi0(arg);
        double ret = matrix_size_os * bessi * one_over_W;
        return ret;
    }

    __inline__ __device__ __host__
    float KaiserBessel(float u, float matrix_size_os,
                       float one_over_W, float beta) {
        float _tmp = 2.0f * u * one_over_W;
        float tmp = _tmp * _tmp;
        float arg = beta * std::sqrt(1.0f - tmp);
        float bessi = bessi0(arg);
        float ret = matrix_size_os * bessi * one_over_W;
        return ret;
    }


    template<class REAL> __inline__ __device__ __host__
    REAL KaiserBessel(
        const Gadgetron::vector_td<REAL, 1> &u,
        const Gadgetron::vector_td<REAL, 1> &matrix_size_os,
        REAL one_over_W,
        const vector_td<REAL, 1> &beta)
    {
        REAL phi_x = KaiserBessel(u.vec[0], matrix_size_os.vec[0], one_over_W, beta[0]);
        return phi_x;
    }

    template<class REAL> __inline__ __device__ __host__
    REAL KaiserBessel(
        const Gadgetron::vector_td<REAL, 2> &u,
        const Gadgetron::vector_td<REAL, 2> &matrix_size_os,
        REAL one_over_W,
        const vector_td<REAL, 2> &beta)
    {
        REAL phi_x = KaiserBessel(u.vec[0], matrix_size_os.vec[0], one_over_W, beta[0]);
        REAL phi_y = KaiserBessel(u.vec[1], matrix_size_os.vec[1], one_over_W, beta[1]);
        return phi_x * phi_y;
    }

    template<class REAL> __inline__ __device__ __host__
    REAL KaiserBessel(
        const Gadgetron::vector_td<REAL, 3> &u,
        const Gadgetron::vector_td<REAL, 3> &matrix_size_os,
        REAL one_over_W,
        const vector_td<REAL, 3> &beta)
    {
        REAL phi_x = KaiserBessel(u.vec[0], matrix_size_os.vec[0], one_over_W, beta[0]);
        REAL phi_y = KaiserBessel(u.vec[1], matrix_size_os.vec[1], one_over_W, beta[1]);
        REAL phi_z = KaiserBessel(u.vec[2], matrix_size_os.vec[2], one_over_W, beta[2]);
        return phi_x * phi_y * phi_z;
    }

    template<class REAL> __inline__ __device__ __host__
    REAL KaiserBessel(
        const Gadgetron::vector_td<REAL, 4> &u,
        const Gadgetron::vector_td<REAL, 4> &matrix_size_os,
        REAL one_over_W,
        const vector_td<REAL, 4> &beta)
    {
        REAL phi_x = KaiserBessel(u.vec[0], matrix_size_os.vec[0], one_over_W, beta[0]);
        REAL phi_y = KaiserBessel(u.vec[1], matrix_size_os.vec[1], one_over_W, beta[1]);
        REAL phi_z = KaiserBessel(u.vec[2], matrix_size_os.vec[2], one_over_W, beta[2]);
        REAL phi_w = KaiserBessel(u.vec[3], matrix_size_os.vec[3], one_over_W, beta[3]);
        return phi_x * phi_y * phi_z * phi_w;
    }


    template<class REAL, unsigned int D, template<class, unsigned int> class K>
    ConvolutionKernel<REAL, D, K>::ConvolutionKernel(REAL width)
      : width_(width)
      , radius_(width / REAL(2))
    {

    }


    template<class REAL, unsigned int D, template<class, unsigned int> class K>
    ConvolutionKernel<REAL, D, K>::~ConvolutionKernel()
    {

    }  


    template<class REAL, unsigned int D, template<class, unsigned int> class K>
    inline REAL ConvolutionKernel<REAL, D, K>::get(const vector_td<REAL, D>& u) const
    {
        if (weak_greater(u, vector_td<REAL, D>(this->radius_)))
            return REAL(0);
        return this->compute(u);
    }


    template<class REAL, unsigned int D, template<class, unsigned int> class K>
    inline REAL ConvolutionKernel<REAL, D, K>::get(REAL r, size_t ax) const
    {
        r = abs(r);
        if (r > this->radius_)
            return REAL(0);
        return this->compute(r, ax);
    }


    template<class REAL, unsigned int D, template<class, unsigned int> class K>
    inline REAL ConvolutionKernel<REAL, D, K>::compute(const vector_td<REAL, D>& u) const
    {
        return static_cast<const K<REAL, D>*>(this)->compute(u);
    }


    template<class REAL, unsigned int D, template<class, unsigned int> class K>
    inline REAL ConvolutionKernel<REAL, D, K>::compute(REAL r, size_t ax) const
    {
        return static_cast<const K<REAL, D>*>(this)->compute(r, ax);
    }


    template<class REAL, unsigned int D, template<class, unsigned int> class K>
    inline REAL ConvolutionKernel<REAL, D, K>::lookup(const vector_td<REAL, D>& u) const
    {
        return static_cast<const K<REAL, D>*>(this)->lookup(u);
    }


    template<class REAL, unsigned int D, template<class, unsigned int> class K>
    inline REAL ConvolutionKernel<REAL, D, K>::lookup(REAL r, size_t ax) const
    {
        return static_cast<const K<REAL, D>*>(this)->lookup(r, ax);
    }


    template<class REAL, unsigned int D, template<class, unsigned int> class K>
    inline REAL ConvolutionKernel<REAL, D, K>::get_width() const
    {
        return width_;
    }


    template<class REAL, unsigned int D, template<class, unsigned int> class K>
    inline REAL ConvolutionKernel<REAL, D, K>::get_radius() const
    {
        return radius_;
    }


    template<class REAL, unsigned int D>
    KaiserKernel<REAL, D>::KaiserKernel(
        const vector_td<unsigned int, D>& matrix_size,
        const vector_td<unsigned int, D>& matrix_size_os,
        REAL width)
      : ConvolutionKernel<REAL, D, KaiserKernel>(width)
    {
        this->matrix_size_ = matrix_size;
        this->matrix_size_os_ = matrix_size_os;
        this->os_factor_ = vector_td<REAL, D>(matrix_size_os) /
                           vector_td<REAL, D>(matrix_size);
        this->width_inv_ = REAL(1) / this->width_;
        this->beta_ = this->compute_beta();
    }


    template<class REAL, unsigned int D>
    KaiserKernel<REAL, D>::KaiserKernel(
        const vector_td<unsigned int, D>& matrix_size,
        REAL os_factor,
        REAL width)
      : ConvolutionKernel<REAL, D, KaiserKernel>(width)
    {
        this->matrix_size_ = matrix_size;
        this->matrix_size_os_ = vector_td<unsigned int, D>(
            vector_td<REAL, D>(matrix_size) * os_factor);
        this->os_factor_ = vector_td<REAL, D>(os_factor);
        this->width_inv_ = REAL(1) / this->width_;
        this->beta_ = this->compute_beta();
    }


    template<class REAL, unsigned int D>
    inline REAL KaiserKernel<REAL, D>::compute(const vector_td<REAL, D>& u) const
    {
        return KaiserBessel<REAL>(
            u,
            vector_td<REAL, D>(this->matrix_size_os_),
            this->width_inv_,
            this->beta_);
    }


    template<class REAL, unsigned int D>
    inline REAL KaiserKernel<REAL, D>::compute(REAL r, size_t ax) const
    {
        return KaiserBessel(
            r,
            static_cast<REAL>(this->matrix_size_os_[ax]),
            this->width_inv_,
            this->beta_[ax]);
    }


    template<class REAL, unsigned int D>
    inline vector_td<REAL, D> KaiserKernel<REAL, D>::get_beta() const
    {
        return this->beta_;
    }


    template<class REAL, unsigned int D>
    vector_td<REAL, D> KaiserKernel<REAL, D>::compute_beta() const
    {
        // Square utility.
        auto sqr = [](auto x) { return x * x; };

        // Compute Kaiser-Bessel beta parameter according to 
        // Beatty et al. IEEE TMI 2005;24(6):799-808.
        vector_td<REAL, D> beta;
        for (unsigned int d = 0; d < D; d++)
        {
            beta[d] = REAL(M_PI) * std::sqrt(
                sqr(this->width_) / sqr(this->os_factor_[d]) *
                sqr(this->os_factor_[d] - REAL(0.5)) - REAL(0.8));
        }

        return beta;
    }


    template<class REAL, unsigned int D>
    JincKernel<REAL, D>::JincKernel(
        float kernelWidth)
      : ConvolutionKernel<REAL, D, JincKernel>(0.0)
    {
        // this->matrix_size_ = matrix_size;
        // this->matrix_size_os_ = matrix_size_os;
        // this->os_factor_ = vector_td<REAL, D>(matrix_size_os) /
        //                    vector_td<REAL, D>(matrix_size);
        // The width is not chosen by the user in this kernel; instead it
        // depends on the oversampling ratio. To keep the kernel isotropic, we
        // take the maximum oversampling factor. This should be very similar in
        // all axes anyway.
        // std::vector<REAL> temp;
        
        // for (auto d=0; d<D; d++)
        //     temp.push_back(kernelWidth);
            //this -> os_factor_.push_back(kernelWidth);
        //this -> os_factor_ = from_std_vector<size_t, 2>(temp);
        this->radius_ = this->norm_radius * kernelWidth;
        this->width_ = this->radius_ * REAL(2);
    }


    // template<class REAL, unsigned int D>
    // JincKernel<REAL, D>::JincKernel(
    //     float kernelWidth)
    //   : ConvolutionKernel<REAL, D, JincKernel>(0.0)
    // {
    //     // this->matrix_size_ = matrix_size;
    //     // this->matrix_size_os_ = vector_td<unsigned int, D>(
    //     //     vector_td<REAL, D>(matrix_size) * os_factor);
    //     // this->os_factor_ = vector_td<REAL, D>(os_factor);
    //     // The width is not chosen by the user in this kernel; instead it
    //     // depends on the oversampling ratio. To keep the kernel isotropic, we
    //     // take the maximum oversampling factor. This should be very similar in
    //     // all axes anyway.
    //     this->radius_ = this->norm_radius * kernelWidth;
    //     this->width_ = this->radius_ * REAL(2);
    // }


    template<class REAL, unsigned int D>
    inline REAL JincKernel<REAL, D>::compute(const vector_td<REAL, D>& u) const
    {
        // This kernel is implemented with circular symmetry, so get the radius
        // for specified coordinates, normalize it and compute based on that.
        REAL r = REAL(0);
        for (unsigned int d = 0; d < D; d++)
        {
            r += u[d] * u[d];
        }
        r = sqrt(r);
        return this->compute(r);
    }


    template<class REAL, unsigned int D>
    inline REAL JincKernel<REAL, D>::compute(REAL r, size_t ax) const
    {
        // Normalize radius.
        r /= this->radius_;

        // Beyond kernel radius, the value is zero.
        if (r >= this->norm_radius)
            return REAL(0);

        // Compute kernel value using polynomial fit.
        // In a perfect world where CUDA supports const-qualified arrays, we
        // would do this.
        // REAL value = this->poly_coefs[0];
        // for (unsigned int i = 1; i <= this->poly_order; i++)
        // {
        //     value += pow(r, static_cast<REAL>(i)) *
        //         this->poly_coefs[i];
        // }

        // Compute kernel value using polynomial fit.
        // Unfortunately, we need to do it this way.
        REAL value = this->poly_coef_0 +
                     this->poly_coef_1 * pow(r, REAL(1)) +
                     this->poly_coef_2 * pow(r, REAL(2)) + 
                     this->poly_coef_3 * pow(r, REAL(3)) +
                     this->poly_coef_4 * pow(r, REAL(4)) +
                     this->poly_coef_5 * pow(r, REAL(5));
            
        // Clip negative values. There shouldn't be any, but just in case.
        if (value < REAL(0))
            value = REAL(0);

        return value;
    }

} // namespace Gadgetron
