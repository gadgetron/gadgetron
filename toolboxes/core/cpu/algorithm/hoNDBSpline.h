/** \file       hoNDBSpline.h
    \brief      N-dimensional inteprolation BSpline implemenation

                The source code is partially from http://bigwww.epfl.ch/thevenaz/interpolation/
                by Philippe Th�venaz

                References:

                [1] P. Th�venaz, T. Blu, M. Unser, "Interpolation Revisited," IEEE Trans on Medical Imaging, Vol 19, 7, 739-758, July 2000.
                [2] M. Unser, A. Aldroubi and M. Eden, "B-Spline Signal Processing: Part I--Theory," IEEE Trans on Signal Processing, Vol 41, 2, 821-832, Feb 1993.
                [3] M. Unser, A. Aldroubi and M. Eden, "B-Spline Signal Processing: Part II--Efficient Design and Applications," IEEE Trans on Signal Processing, Vol 41, 2, 834-848, Feb 1993.

    \author     Hui Xue
*/

#pragma once

#include "hoNDArray.h"
#include "hoNDImage.h"

namespace Gadgetron
{
    template <typename T, unsigned int D, typename coord_type = double>
    class hoNDBSpline
    {
    public:

        typedef hoNDBSpline<T, D> Self;

        typedef T element_type;
        typedef T value_type;

        /// type for bspline computation, can be 'float' or 'double'
        typedef typename realType<T>::Type bspline_float_type;

        typedef hoNDArray<T> ArrayType;
        typedef hoNDImage<T, D> ImageType;

        hoNDBSpline() {}
        ~hoNDBSpline() {}

        /// compute BSpline coefficient
        bool computeBSplineCoefficients(const hoNDArray<T>& data, unsigned int SplineDegree, hoNDArray<T>& coeff);
        bool computeBSplineCoefficients(const hoNDImage<T, D>& data, unsigned int SplineDegree, hoNDArray<T>& coeff);

        bool computeBSplineCoefficients(const T* data, const std::vector<size_t>& dimension, unsigned int SplineDegree, T* coeff);
        bool computeBSplineCoefficients(const T* data, size_t len, unsigned int SplineDegree, T* coeff);
        bool computeBSplineCoefficients(const T* data, size_t sx, size_t sy, unsigned int SplineDegree, T* coeff);
        bool computeBSplineCoefficients(const T* data, size_t sx, size_t sy, size_t sz, unsigned int SplineDegree, T* coeff);
        bool computeBSplineCoefficients(const T* data, size_t sx, size_t sy, size_t sz, size_t st, unsigned int SplineDegree, T* coeff);
        bool computeBSplineCoefficients(const T* data, size_t sx, size_t sy, size_t sz, size_t st, size_t sp, unsigned int SplineDegree, T* coeff);
        bool computeBSplineCoefficients(const T* data, size_t sx, size_t sy, size_t sz, size_t st, size_t sp, size_t sq, unsigned int SplineDegree, T* coeff);
        bool computeBSplineCoefficients(const T* data, size_t sx, size_t sy, size_t sz, size_t st, size_t sp, size_t sq, size_t sr, unsigned int SplineDegree, T* coeff);
        bool computeBSplineCoefficients(const T* data, size_t sx, size_t sy, size_t sz, size_t st, size_t sp, size_t sq, size_t sr, size_t ss, unsigned int SplineDegree, T* coeff);
        bool computeBSplineCoefficients(const T* data, size_t sx, size_t sy, size_t sz, size_t st, size_t sp, size_t sq, size_t sr, size_t ss, size_t su, unsigned int SplineDegree, T* coeff);

        /// evaluate BSpline
        /// derivative: can be 0/1/2, for 0-order, first-order and second-order derivative

        T evaluateBSpline(const T* coeff, const std::vector<size_t>& dimension, unsigned int SplineDegree, 
                        const std::vector<unsigned int>& derivative, 
                        const coord_type* pos);

        T evaluateBSpline(const T* coeff, const std::vector<size_t>& dimension, unsigned int SplineDegree, 
                        const std::vector<unsigned int>& derivative, 
                        const std::vector<coord_type>& pos);

        T evaluateBSpline(const T* coeff, size_t len, unsigned int SplineDegree, 
                        unsigned int dx, 
                        coord_type x);

        T evaluateBSpline(const T* coeff, size_t sx, size_t sy, unsigned int SplineDegree, 
                        unsigned int dx, unsigned int dy, 
                        coord_type x, coord_type y);

        T evaluateBSpline(const T* coeff, size_t sx, size_t sy, size_t sz, unsigned int SplineDegree, 
                        unsigned int dx, unsigned int dy, unsigned int dz, 
                        coord_type x, coord_type y, coord_type z);

        T evaluateBSpline(const T* coeff, size_t sx, size_t sy, size_t sz, size_t st, unsigned int SplineDegree, 
                        unsigned int dx, unsigned int dy, unsigned int dz, unsigned int dt, 
                        coord_type x, coord_type y, coord_type z, coord_type t);

        T evaluateBSpline(const T* coeff, size_t sx, size_t sy, size_t sz, size_t st, size_t sp, unsigned int SplineDegree, 
                        unsigned int dx, unsigned int dy, unsigned int dz, unsigned int dt, unsigned int dp, 
                        coord_type x, coord_type y, coord_type z, coord_type t, coord_type p);

        T evaluateBSpline(const T* coeff, size_t sx, size_t sy, size_t sz, size_t st, size_t sp, size_t sq, unsigned int SplineDegree, 
                        unsigned int dx, unsigned int dy, unsigned int dz, unsigned int dt, unsigned int dp, unsigned int dq, 
                        coord_type x, coord_type y, coord_type z, coord_type t, coord_type p, coord_type q);

        T evaluateBSpline(const T* coeff, size_t sx, size_t sy, size_t sz, size_t st, size_t sp, size_t sq, size_t sr, unsigned int SplineDegree, 
                        unsigned int dx, unsigned int dy, unsigned int dz, unsigned int dt, unsigned int dp, unsigned int dq, unsigned int dr, 
                        coord_type x, coord_type y, coord_type z, coord_type t, coord_type p, coord_type q, coord_type r);

        T evaluateBSpline(const T* coeff, size_t sx, size_t sy, size_t sz, size_t st, size_t sp, size_t sq, size_t sr, size_t ss, unsigned int SplineDegree, 
                        unsigned int dx, unsigned int dy, unsigned int dz, unsigned int dt, unsigned int dp, unsigned int dq, unsigned int dr, unsigned int ds, 
                        coord_type x, coord_type y, coord_type z, coord_type t, coord_type p, coord_type q, coord_type r, coord_type s);

        T evaluateBSpline(const T* coeff, size_t sx, size_t sy, size_t sz, size_t st, size_t sp, size_t sq, size_t sr, size_t ss, size_t su, unsigned int SplineDegree, 
                        unsigned int dx, unsigned int dy, unsigned int dz, unsigned int dt, unsigned int dp, unsigned int dq, unsigned int dr, unsigned int ds, unsigned int du, 
                        coord_type x, coord_type y, coord_type z, coord_type t, coord_type p, coord_type q, coord_type r, coord_type s, coord_type u);


        /// evaluate BSpline with pre-computed weights
        T evaluateBSpline(const T* coeff, size_t sx, size_t sy, unsigned int SplineDegree, 
                        bspline_float_type* xWeight, bspline_float_type* yWeight, 
                        coord_type x, coord_type y);

        T evaluateBSpline(const T* coeff, size_t sx, size_t sy, size_t sz, unsigned int SplineDegree, 
                        bspline_float_type* xWeight, bspline_float_type* yWeight, bspline_float_type* zWeight, 
                        coord_type x, coord_type y, coord_type z);

        T evaluateBSpline(const T* coeff, size_t sx, size_t sy, size_t sz, size_t st, unsigned int SplineDegree, 
                        bspline_float_type* xWeight, bspline_float_type* yWeight, bspline_float_type* zWeight, bspline_float_type* tWeight, 
                        coord_type x, coord_type y, coord_type z, coord_type t);

        T evaluateBSpline(const T* coeff, const std::vector<size_t>& dimension, unsigned int SplineDegree,
                        bspline_float_type** weight, const coord_type* pos);

        T evaluateBSpline(const T* coeff, const std::vector<size_t>& dimension, unsigned int SplineDegree,
                        bspline_float_type** weight, const std::vector<coord_type>& pos);

        /// compute the BSpline based derivative for an ND array
        /// derivative indicates the order of derivatives for every dimension
        bool computeBSplineDerivative(const hoNDArray<T>& data, const hoNDArray<T>& coeff, unsigned int SplineDegree, const std::vector<unsigned int>& derivative, hoNDArray<T>& deriv);
        bool computeBSplineDerivative(const hoNDImage<T,D>& data, const hoNDArray<T>& coeff, unsigned int SplineDegree, const 
          std::vector<unsigned int>& derivative, hoNDImage<T,D>& deriv);
        /// evaluate BSpline at an array of points (not at array grid)
        /// pts: N x D array, every line is a point to evaluate BSpline
        /// deriv: N x 1 array, evaluation results
        bool computeBSplineDerivativePoints(const hoNDArray<T>& pts, const hoNDArray<T>& coeff, unsigned int SplineDegree, const std::vector<unsigned int>& derivative, hoNDArray<T>& deriv);

        /// print out the image information
        void print(std::ostream& os) const;

    protected:

        /// these BSpline coefficients paramerters are modified from http://bigwww.epfl.ch/thevenaz/interpolation/
        static void ConvertToInterpolationCoefficients(
                                                        T               c[],            /* input samples --> output coefficients */
                                                        size_t          DataLength,     /* number of samples or coefficients */
                                                        bspline_float_type          z[],            /* poles */
                                                        long            NbPoles,        /* number of poles */
                                                        bspline_float_type          Tolerance       /* admissible relative error */ 
                                                      );

        static T InitialCausalCoefficient(
                                            T           c[],                /* coefficients */
                                            size_t      DataLength,         /* number of coefficients */
                                            bspline_float_type      z,                  /* actual pole */
                                            bspline_float_type      Tolerance           /* admissible relative error */
                                         );

        static T InitialAntiCausalCoefficient(
                                                T           c[],                /* coefficients */
                                                size_t      DataLength,         /* number of samples or coefficients */
                                                bspline_float_type      z                   /* actual pole */
                                             );

        static void Pole(bspline_float_type* pole, unsigned int SplineDegree, unsigned int& NbPoles);

        /// BSpline function
        /// this function implements the symmetrical BSpline function of SplineDegree n
        /// Equation 2.6 of reference [2]
        static bspline_float_type BSpline(bspline_float_type x, unsigned int SplineDegree);

        /// compute the discrete BSpline value
        /// this function is modified from the source code at http://bigwww.epfl.ch/thevenaz/interpolation/
        static void BSplineDiscrete(bspline_float_type x, unsigned int SplineDegree, bspline_float_type* weight, long long* xIndex);

        /// compute the discrete BSpline value with the first order derivative
        static void BSplineDiscreteFirstOrderDerivative(bspline_float_type x, unsigned int SplineDegree, bspline_float_type* weight, long long* xIndex);
        /// compute the discrete BSpline value with the second order derivative
        static void BSplineDiscreteSecondOrderDerivative(bspline_float_type x, unsigned int SplineDegree, bspline_float_type* weight, long long* xIndex);

        /// compute BSpline interpolation locations
        /// xIndex has at least SplineDegree elements
        static void BSplineInterpolationLocation(bspline_float_type x, unsigned int SplineDegree, long long* xIndex);

        /// apply mirror boundary condition for interpolation locations
        static void BSplineInterpolationMirrorBoundaryCondition(unsigned int SplineDegree, long long* xIndex, size_t Width);

        /// compute the derivative of BSpline
        /// first order derivative dBSpline(x, SplineDegree)/dx = BSpline(x+0.5, SplineDegree-1) - BSpline(x-0.5, SplineDegree-1)
        static bspline_float_type BSplineFirstOrderDerivative(bspline_float_type x, unsigned int SplineDegree);
        /// second order derivative d2BSpline(x, SplineDegree)/dx2 = BSpline(x+1, SplineDegree-2) + BSpline(x-1, SplineDegree-2) - 2*BSpline(x, SplineDegree-2)
        static bspline_float_type BSplineSecondOrderDerivative(bspline_float_type x, unsigned int SplineDegree);

        /// compute BSpline interpolation locations and weights
        static void computeBSplineInterpolationLocationsAndWeights(size_t len, unsigned int SplineDegree, unsigned int dx, coord_type x, bspline_float_type* weight, long long* xIndex);
    };
}

#include "hoNDBSpline.hxx"
