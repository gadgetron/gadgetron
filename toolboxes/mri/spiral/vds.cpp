#include "vds.h"

#include "hoNDArray_math.h"
#include "vector_td_utilities.h"
#include <boost/math/constants/constants.hpp>
#include <boost/range/algorithm.hpp>
#include <boost/range/combine.hpp>
#include <math.h>
#include <stdio.h>

constexpr double GAMMA = 4258.0; /* Hz/G */
constexpr double PI = 3.1416;    // boost::math::constants::pi<double>();

/* #define TESTCODE 	For testing as regular C code... */

/*
  %
  %	VARIABLE DENSITY SPIRAL GENERATION:
  %	----------------------------------
  %
  %	This is a general description of how the following C code
  %	works.  This text is taken from a matlab script, vds.m, from
  %	which the C code was derived.  However, note that the C code
  %	runs considerably faster.
  %
  %
  %	Function generates variable density spiral which traces
  %	out the trajectory
  %
  %			k(t) = r(t) exp(i*q(t)), 		[1]
  %
  %	Where q IS THE SAME AS theta, and r IS THE SAME AS kr.
  %
  %		r and q are chosen to satisfy:
  %
  %		1) Maximum gradient amplitudes and slew rates.
  %		2) Maximum gradient due to FOV, where FOV can
  %		   vary with k-space radius r, as
  %
  %			FOV(r) = F0 + F1*r + F2*r*r 		[2]
  %
  %
  %	INPUTS:
  %	-------
  %	smax = maximum slew rate G/cm/s
  %	gmax = maximum gradient G/cm (limited by Gmax or FOV)
  %	T = sampling period (s) for gradient AND acquisition.
  %	N = number of interleaves.
  %	F0,F1,F2 = FOV coefficients with respect to r - see above.
  %	rmax= value of k-space radius at which to stop (cm^-1).
  %		rmax = 1/(2*resolution);
  %
  %
  %	OUTPUTS:
  %	--------
  %	k = k-space trajectory (kx+iky) in cm-1.
  %	g = gradient waveform (Gx+iGy) in G/cm.
  %	s = derivative of g (Sx+iSy) in G/cm/s.
  %	time = time points corresponding to above (s).
  %	r = k-space radius vs time (used to design spiral)
  %	theta = atan2(ky,kx) = k-space angle vs time.
  %
  %
  %	METHODS:
  %	--------
  %	Let r1 and r2 be the first derivatives of r in [1].
  %	Let q1 and q2 be the first derivatives of theta in [1].
  %	Also, r0 = r, and q0 = theta - sometimes both are used.
  %	F = F(r) defined by F0,F1,F2.
  %
  %	Differentiating [1], we can get G = a(r0,r1,q0,q1,F)
  %	and differentiating again, we get S = b(r0,r1,r2,q0,q1,q2,F)
  %
  %	(functions a() and b() are reasonably easy to obtain.)
  %
  %	FOV limits put a constraint between r and q:
  %
  %		dr/dq = N/(2*pi*F)				[3]
  %
  %	We can use [3] and the chain rule to give
  %
  %		q1 = 2*pi*F/N * r1				[4]
  %
  %	and
  %
  %		q2 = 2*pi/N*dF/dr*r1^2 + 2*pi*F/N*r2		[5]
  %
  %
  %
  %	Now using [4] and [5], we can substitute for q1 and q2
  %	in functions a() and b(), giving
  %
  %		G = c(r0,r1,F)
  %	and 	S = d(r0,r1,r2,F,dF/dr)
  %
  %
  %	Using the fact that the spiral should be either limited
  %	by amplitude (Gradient or FOV limit) or slew rate, we can
  %	solve
  %		|c(r0,r1,F)| = |Gmax|  				[6]
  %
  %	analytically for r1, or
  %
  %	  	|d(r0,r1,r2,F,dF/dr)| = |Smax|	 		[7]
  %
  %	analytically for r2.
  %
  %	[7] is a quadratic equation in r2.  The smaller of the
  %	roots is taken, and the real part of the root is used to
  %	avoid possible numeric errors - the roots should be real
  %	always.
  %
  %	The choice of whether or not to use [6] or [7], and the
  %	solving for r2 or r1 is done by calcthetadotdot().
  %
  %	Once the second derivative of theta(q) or r is obtained,
  %	it can be integrated to give q1 and r1, and then integrated
  %	again to give q and r.  The gradient waveforms follow from
  %	q and r.
  %
  %	Brian Hargreaves -- Sept 2000.
  %
  %
*/

namespace Gadgetron {

/* ----------------------------------------------------------------------- */
void calcthetadotdot(double slewmax, double gradmax, double kr, double krdot, double Tgsample, double Tdsample,
                     int Ninterleaves, double* fov, int numfov, double* thetadotdot, double* krdotdot, double krmax)
/*
 * Function calculates the 2nd derivative of kr and theta at each
 * sample point within calc_vds().  ie, this is the iterative loop
 * for calc_vds.  See the text at the top of this file for more details
 * */

// double slewmax;		/*	Maximum slew rate, G/cm/s		*/
// double gradmax;		/* 	maximum gradient amplitude, G/cm	*/
// double kr;		/* 	Current kr. */
// double krdot;		/*	Current krdot. */
// double Tgsample;	/*	Gradient Sample period (s) 	*/
// double Tdsample;	/*	Data Sample period (s) 		*/
// int Ninterleaves;	/*	Number of interleaves			*/
// double *fov;		/*	FOV coefficients		*/
// int numfov;		/*	Number of FOV coefficients		*/
// double *thetadotdot;	/*	[output] 2nd derivative of theta.	*/
// double *krdotdot;	/*	[output] 2nd derivative of kr		*/

/* ----------------------------------------------------------------------- */
{
    double fovval = 0;    /* FOV for this value of kr	*/
    double dfovdrval = 0; /* dFOV/dkr for this value of kr	*/
    double gmaxfov;       /* FOV-limited Gmax.	*/
    double maxkrdot;
    int count;

    double tpf;   /* Used to simplify expressions. */
    double tpfsq; /* 	" 		"        */

    double qdfA, qdfB, qdfC; /* Quadratic formula coefficients */
    double rootparta, rootpartb;

    /* Calculate the actual FOV and dFOV/dkr for this R,
     * based on the fact that the FOV is expressed
     * as a polynomial in kr.*/

    for (count = 0; count < numfov; count++) {
        fovval = fovval + fov[count] * pow(kr / krmax, count);
        if (count > 0)
            dfovdrval = dfovdrval + count * fov[count] * pow(kr / krmax, count - 1) / krmax;
    }

    /* Calculate FOV limit on gmax.  This is the rate of motion along
     * a trajectory, and really should not be a limitation.  Thus,
     * it is reasonable to comment out the following lines. */

    gmaxfov = 1 / GAMMA / fovval / Tdsample;
    if (gradmax > gmaxfov)
        gradmax = gmaxfov;

    /* Maximum dkr/dt, based on gradient amplitude.  */

    maxkrdot = sqrt(pow(GAMMA * gradmax, 2) / (1 + pow(2 * PI * fovval * kr / Ninterleaves, 2)));

    /* These two are just to simplify expressions below */
    tpf = 2 * PI * fovval / Ninterleaves;
    tpfsq = pow(tpf, 2);

    if (krdot > maxkrdot) /* Then choose krdotdot so that krdot is in range */
    {
        *krdotdot = (maxkrdot - krdot) / Tgsample;
    } else /* Choose krdotdot based on max slew rate limit. */
    {

        /* Set up for quadratic formula solution. */

        qdfA = 1 + tpfsq * kr * kr;
        qdfB = 2 * tpfsq * kr * krdot * krdot + 2 * tpfsq / fovval * dfovdrval * kr * kr * krdot * krdot;
        qdfC = pow(tpfsq * kr * krdot * krdot, 2) + 4 * tpfsq * pow(krdot, 4) +
               pow(tpf * dfovdrval / fovval * kr * krdot * krdot, 2) +
               4 * tpfsq * dfovdrval / fovval * kr * pow(krdot, 4) - pow(GAMMA * slewmax, 2);

        rootparta = -qdfB / (2 * qdfA);
        rootpartb = qdfB * qdfB / (4 * qdfA * qdfA) - qdfC / qdfA;

        if (rootpartb < 0) /* Safety check - if complex, take real part.*/

            *krdotdot = rootparta;

        else
            *krdotdot = rootparta + sqrt(rootpartb);

        /* Could check resulting slew rate here, as in q2r21.m. */
    }

    /* Calculate thetadotdot */

    *thetadotdot = tpf * dfovdrval / fovval * krdot * krdot + tpf * (*krdotdot);
}

// /* ----------------------------------------------------------------------- */
// void
// calc_vds(double slewmax, double gradmax, double Tgsample, double Tdsample, int Ninterleaves,
//          double *fov, int numfov, double krmax,
//          int ngmax, double **xgrad, double **ygrad, int *numgrad)

// /*	Function designs a variable-density spiral gradient waveform
//  *	that is defined by a number of interleaves, resolution (or max number
//  *	of samples), and field-of-view.
//  *	The field-of-view is a polynomial function of the
//  *	k-space radius, so fov is an array of coefficients so that
//  *
//  *	FOV = fov[0]+fov[1]*kr+fov[2]*kr^2+ ... +fov[numfov-1]*kr^(numfov-1)
//  *
//  * 	Gradient design is subject to a constant-slew-rate-limit model,
//  * 	with maximum slew rate slewmax, and maximum gradient amplitude
//  * 	of gradmax.
//  *
//  * 	Tgsample is the gradient sampling rate, and Tdsample is the data
//  * 	sampling rate.  It is highly recommended to OVERSAMPLE the gradient
//  * 	in the design to make the integration more stable.
//  *
//  * */

// //double slewmax;		/*	Maximum slew rate, G/cm/s		*/
// //double gradmax;		/* 	maximum gradient amplitude, G/cm	*/
// //double Tgsample;	/*	Gradient Sample period (s)		*/
// //double Tdsample;	/*	Data Sample period (s)			*/
// //int Ninterleaves;	/*	Number of interleaves			*/
// //double *fov;		/*	FOV coefficients		*/
// //int numfov;		/*	Number of FOV coefficients		*/
// //double krmax;		/*	Maximum k-space extent (/cm)		*/
// //int ngmax;		/*	Maximum number of gradient samples	*/
// //double **xgrad;		/* 	[output] X-component of gradient (G/cm) */
// //double **ygrad;		/*	[output] Y-component of gradient (G/cm)	*/
// //int *numgrad;		/* 	[output] Number of gradient samples */

// /* ----------------------------------------------------------------------- */
// {
//     int gradcount = 0;

//     double kr = 0;            /* Current value of kr	*/
//     double krdot = 0;        /* Current value of 1st derivative of kr */
//     double krdotdot = 0;        /* Current value of 2nd derivative of kr */

//     double theta = 0;            /* Current value of theta */
//     double thetadot = 0;        /* Current value of 1st derivative of theta */
//     double thetadotdot = 0;        /* Current value of 2nd derivative of theta */

//     double lastkx = 0;        /* x-component of last k-location. */
//     double lastky = 0;        /* y-component of last k-location */
//     double kx, ky;            /* x and y components of current k-location */

//     double *gxptr, *gyptr;        /* Pointers to gradient variables. */

//     /* First just find the gradient length. */

//     while ((kr < krmax) && (gradcount < ngmax)) {
//         calcthetadotdot(slewmax, gradmax, kr, krdot, Tgsample, Tdsample,
//                         Ninterleaves, fov, numfov, &thetadotdot, &krdotdot);

//         /* Integrate to obtain new values of kr, krdot, theta and thetadot:*/

//         thetadot = thetadot + thetadotdot * Tgsample;
//         theta = theta + thetadot * Tgsample;

//         krdot = krdot + krdotdot * Tgsample;
//         kr = kr + krdot * Tgsample;

//         gradcount++;

//     }

//     /* Allocate memory for gradients. */

//     *numgrad = gradcount;

//     //*xgrad = (double *)malloc(*numgrad*sizeof(double));
//     //*ygrad = (double *)malloc(*numgrad*sizeof(double));

//     *xgrad = new double[*numgrad];
//     *ygrad = new double[*numgrad];

//     /* Reset parameters */

//     kr = 0;
//     krdot = 0;
//     theta = 0;
//     thetadot = 0;
//     gradcount = 0;
//     gxptr = *xgrad;
//     gyptr = *ygrad;

//     /* Now re-calculate gradient to find length. */

//     while ((kr < krmax) && (gradcount < ngmax)) {
//         calcthetadotdot(slewmax, gradmax, kr, krdot, Tgsample, Tdsample,
//                         Ninterleaves, fov, numfov, &thetadotdot, &krdotdot);

//         /* Integrate to obtain new values of kr, krdot, theta and thetadot:*/

//         thetadot = thetadot + thetadotdot * Tgsample;
//         theta = theta + thetadot * Tgsample;

//         krdot = krdot + krdotdot * Tgsample;
//         kr = kr + krdot * Tgsample;

//         /* Define current gradient values from kr and theta. */

//         kx = kr * cos(theta);
//         ky = kr * sin(theta);
//         *gxptr++ = (1 / GAMMA / Tgsample) * (kx - lastkx);
//         *gyptr++ = (1 / GAMMA / Tgsample) * (ky - lastky);
//         lastkx = kx;
//         lastky = ky;

//         gradcount++;
//     }

// }
void calc_vds(double slewmax, double gradmax, double Tgsample, double Tdsample, int Ninterleaves, double* fov,
              int numfov, double krmax, int ngmax, double** xgrad, double** ygrad, int* numgrad, int numgradmax)

/*	Function designs a variable-density spiral gradient waveform
 *	that is defined by a number of interleaves, resolution (or max number
 *	of samples), and field-of-view.
 *	The field-of-view is a polynomial function of the
 *	k-space radius, so fov is an array of coefficients so that
 *
 *	FOV = fov[0]+fov[1]*kr+fov[2]*kr^2+ ... +fov[numfov-1]*kr^(numfov-1)
 *
 * 	Gradient design is subject to a constant-slew-rate-limit model,
 * 	with maximum slew rate slewmax, and maximum gradient amplitude
 * 	of gradmax.
 *
 * 	Tgsample is the gradient sampling rate, and Tdsample is the data
 * 	sampling rate.  It is highly recommended to OVERSAMPLE the gradient
 * 	in the design to make the integration more stable.
 *
 *  Returns false if failed
 */

// double slewmax;		/*	Maximum slew rate, G/cm/s		*/
// double gradmax;		/* 	maximum gradient amplitude, G/cm	*/
// double Tgsample;		/*	Gradient Sample period (s)		*/
// double Tdsample;		/*	Data Sample period (s)			*/
// int Ninterleaves;		/*	Number of interleaves			*/
// double *fov;			/*	FOV coefficients		*/
// int numfov;			/*	Number of FOV coefficients		*/
// double krmax;			/*	Maximum k-space extent (/cm)		*/
// int ngmax;			/*	Maximum number of gradient samples	*/
// double **xgrad;		/* 	[output] X-component of gradient (G/cm) */
// double **ygrad;		/*	[output] Y-component of gradient (G/cm)	*/
// int *numgrad;			/* 	[output] Number of gradient samples */
// int numgradmax;		/*  rbhpg Maximum number of gradient samples, as arrays are allocated by calling
// function
// */

/* ----------------------------------------------------------------------- */
{
    int gradcount = 0;

    double kr = 0;       /* Current value of kr	*/
    double krdot = 0;    /* Current value of 1st derivative of kr */
    double krdotdot = 0; /* Current value of 2nd derivative of kr */

    double theta = 0;       /* Current value of theta */
    double thetadot = 0;    /* Current value of 1st derivative of theta */
    double thetadotdot = 0; /* Current value of 2nd derivative of theta */

    double lastkx = 0; /* x-component of last k-location. */
    double lastky = 0; /* y-component of last k-location */
    double kx, ky;     /* x and y components of current k-location */

    double *gxptr, *gyptr; /* Pointers to gradient variables. */

    /* First just find the gradient length. */

    while ((kr < krmax) && (gradcount < ngmax)) {
        calcthetadotdot(slewmax, gradmax, kr, krdot, Tgsample, Tdsample, Ninterleaves, fov, numfov, &thetadotdot,
                        &krdotdot, krmax);

        /* Integrate to obtain new values of kr, krdot, theta and thetadot:*/

        thetadot = thetadot + thetadotdot * Tgsample;
        theta = theta + thetadot * Tgsample;

        krdot = krdot + krdotdot * Tgsample;
        kr = kr + krdot * Tgsample;

        gradcount++;
    }

    /* Allocate memory for gradients. */

    *numgrad = gradcount;

    // KSN start
    // don't allocate memory for gradients, that responsibility is shifted to the parent function
    // *xgrad = new double[*numgrad*sizeof(double)];
    // *ygrad = new double[*numgrad*sizeof(double)];

    // KSN end

    /* Reset parameters */

    kr = 0;
    krdot = 0;
    theta = 0;
    thetadot = 0;
    gradcount = 0;
    gxptr = *xgrad;
    gyptr = *ygrad;
    int gradcount_real = 0;
    /* Now re-calculate gradient to find length. */

    while ((kr < krmax) && (gradcount < ngmax)) {
        calcthetadotdot(slewmax, gradmax, kr, krdot, Tgsample, Tdsample, Ninterleaves, fov, numfov, &thetadotdot,
                        &krdotdot, krmax);

        /* Integrate to obtain new values of kr, krdot, theta and thetadot:*/

        thetadot = thetadot + thetadotdot * Tgsample;
        theta = theta + thetadot * Tgsample;

        krdot = krdot + krdotdot * Tgsample;
        kr = kr + krdot * Tgsample;
        // std::cout<<"krdotdot = "<<krdotdot<<"thetadotdot = %d"<<thetadotdot<<std::endl;
        /* Define current gradient values from kr and theta. */
        if ((gradcount - 2) % 8 == 0) {
            kx = kr * cos(theta);
            ky = kr * sin(theta);
            *gxptr++ = (1 / GAMMA / Tdsample) * (kx - lastkx);
            *gyptr++ = (1 / GAMMA / Tdsample) * (ky - lastky);
            lastkx = kx;
            lastky = ky;
            gradcount_real++;
        }

        gradcount++;
    }
    *numgrad = gradcount_real;

}

/* ----------------------------------------------------------------------- */
void calc_traj(double* xgrad, double* ygrad, int ngrad, int Nints, double Tgsamp, double krmax, double** x_trajectory,
               double** y_trajectory,
               double** weights) //, double** y_weights)
/*
 *inputs:
 *      xgrad   X gradient waveform
 *      ygrad   Y gradient waveform
 *      ngrad   number of gradient samples
 *      Nints   number of interleaves
 *      Tgsamp  sampling time for gradients
 *
 *outputs:
 *      x_trajectory    X position in k-space
 *      y_trajectory    Y position in k-space
 *      x_weights       X weighting
 *      y_weights       Y weighting
 *
 **/
{

    *x_trajectory = new double[(ngrad * Nints)];
    *y_trajectory = new double[(ngrad * Nints)];
    *weights = new double[(ngrad * Nints)];

    double* txptr = *x_trajectory;
    double* typtr = *y_trajectory;
    double* wptr = *weights;

    for (int inter = 0; inter < Nints; inter++) {
        double rotation = (inter * 2 * PI) / Nints;
        double x_tr = 0;
        double y_tr = 0;
        float x_temp, y_temp;
        for (int gradcount = 0; gradcount < ngrad; gradcount++) {
            if (gradcount > 0) {
                x_tr += (GAMMA)*xgrad[gradcount - 1] * Tgsamp;
                y_tr += (GAMMA)*ygrad[gradcount - 1] * Tgsamp;
            }

            x_temp = (x_tr * cos(rotation)) + (y_tr * sin(rotation));
            y_temp = -(x_tr * sin(rotation)) + (y_tr * cos(rotation));
            *(txptr++) = x_temp / krmax;
            *(typtr++) = y_temp / krmax;

            // abs(g(:)
            double abs_w = sqrt((pow(xgrad[gradcount], 2)) + (pow(ygrad[gradcount], 2)));
            double ang_g = xgrad[gradcount] == 0.0 ? PI / 2 : atan2(ygrad[gradcount], xgrad[gradcount]);

            double ang_t = x_tr == 0.0 ? PI / 2 : atan2(y_tr, x_tr);

            double tp_w = sin(ang_g - ang_t);
            tp_w = abs_w * abs(tp_w);
            *wptr++ = tp_w;
        }
    }
}
hoNDArray<float> calculate_weights(const hoNDArray<floatd2>& gradients, const hoNDArray<floatd2>& trajectories) {

    hoNDArray<float> weights(gradients.dimensions());

    boost::transform(gradients, trajectories, weights.begin(), [](const floatd2& gradient, const floatd2& trajectory) {
        auto abs_w = norm(gradient);
        auto ang_g = atan2(gradient[1], gradient[0]);
        auto ang_t = atan2(trajectory[1], trajectory[0]);
        return abs(sin(ang_g - ang_t)) * abs_w * 2;
    });

    return weights;
}

hoNDArray<float> calculate_weights_Hoge(const hoNDArray<floatd2>& gradients, const hoNDArray<floatd2>& trajectories) {

    hoNDArray<float> weights(gradients.dimensions());

    boost::transform(gradients, trajectories, weights.begin(), [](const floatd2& gradient, const floatd2& trajectory) {
        auto abs_g = norm(gradient);
        auto abs_t = norm(trajectory);
        auto ang_g = atan2(gradient[1], gradient[0]);
        auto ang_t = atan2(trajectory[1], trajectory[0]);
        return abs(cos(ang_g - ang_t)) * abs_g * abs_t;
    });

    return weights;
}

hoNDArray<floatd2> calculate_trajectories(const hoNDArray<floatd2>& gradients, float sample_time, float krmax) {
    const int nints = gradients.get_size(1);
    const int ngrad = gradients.get_size(0);

    auto trajectory = hoNDArray<floatd2>(gradients.dimensions());

    for (int interleave = 0; interleave < nints; interleave++) {
        trajectory(0, interleave)[0] = float(GAMMA) * gradients(0, interleave)[0] * sample_time;
        trajectory(0, interleave)[1] = float(GAMMA) * gradients(0, interleave)[1] * sample_time;

        for (int gradcount = 1; gradcount < ngrad; gradcount++) {
            trajectory(gradcount, interleave)[0] =
                trajectory(gradcount - 1, interleave)[0] + float(GAMMA) * gradients(gradcount, interleave)[0] * sample_time;
                            trajectory(gradcount, interleave)[1] =
                trajectory(gradcount - 1, interleave)[1] + float(GAMMA) * gradients(gradcount, interleave)[1] * sample_time;
        }
    }

    boost::transform(trajectory, trajectory.begin(), [&](auto t) { return t / (2 * krmax); });
    return trajectory;
}

hoNDArray<floatd2> calculate_vds(double slewmax, double gradmax, double Tgsample, double Tdsample, int Ninterleaves,
                                 double* fov, int numfov, double krmax, int ngmax, int max_nsamples) {

    std::vector<double> x_grad(ngmax, 0.0);
    std::vector<double> y_grad(ngmax, 0.0);

    // rbhpg oversampling of spiral gradient by fixed factor 5
    int oversample_factor = int(8);
    std::vector<double> xgradovs(int32_t(oversample_factor) * int32_t(ngmax), 0.0);
    std::vector<double> ygradovs(int32_t(oversample_factor) * int32_t(ngmax), 0.0);
    double* x_ptr = x_grad.data();
    double* y_ptr = y_grad.data();
    double* x_ptr_ovs = xgradovs.data();
    double* y_ptr_ovs = ygradovs.data();
    int Nints = 0;

    calc_vds(slewmax, gradmax, Tgsample / float(oversample_factor), Tdsample, Ninterleaves, fov, numfov, krmax, ngmax,
             &x_ptr, &y_ptr, &Nints, ngmax * oversample_factor);

    // int i, j;
    // j = 0;
    // for (i = 0; i < Nints; i++) {
    //     if ((i % oversample_factor) == 0) {
    //         x_ptr[j] = xgradovs[i + oversample_factor / 2 - 1];
    //         y_ptr[j] = ygradovs[i + oversample_factor / 2 - 1];
    //         j++;
    //     }
    // }
    //Nints = j;
    GDEBUG_STREAM("Nints (VDS):" << Nints);
    hoNDArray<floatd2> gradient(std::min(Nints, max_nsamples));
    size_t nsamples = gradient.get_number_of_elements();
    std::transform(x_ptr, x_ptr + nsamples, y_ptr, gradient.begin(),
                   [](auto x, auto y) { return -floatd2(float(x), float(y)); });

    // delete[] x_ptr;
    // delete[] y_ptr;
    return gradient;
}

hoNDArray<floatd2> create_rotations(const hoNDArray<floatd2>& trajectories, int nints) {
    auto dims = trajectories.dimensions();
    dims.push_back(nints);
    hoNDArray<floatd2> result(dims);

    int ngrad = trajectories.get_number_of_elements();
    for (int inter = 0; inter < nints; inter++) {
        double rotation = (inter * 2 * PI) / nints;
        for (int gradcount = 0; gradcount < ngrad; gradcount++) {
            floatd2 point = trajectories[gradcount];
            result(gradcount, inter) = floatd2(point[0] * cos(rotation) + point[1] * sin(rotation),
                                               -point[0] * sin(rotation) + point[1] * cos(rotation));
        }
    }
    return result;
}
} // namespace Gadgetron