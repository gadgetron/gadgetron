/** \file   cmr_radial_thickening.h
	\brief  Implement radial thickening method of computing cardiac strain
	\author Angela Gao
	\date   July 29, 2019
*/

#include "cmr_radial_thickening.h"
#include "hoNDArray_reductions.h"
#include "hoNDArray_elemwise.h"
#include "hoNDInterpolator.h"
#include "hoNDBoundaryHandler.h"
#define _USE_MATH_DEFINES
#include <math.h>

namespace Gadgetron {

	//template <typename T>
	//void linspace(const T start, const T end, const int num, hoNDArray<T>& res)
	//{
	//	try
	//	{
	//		res.create(num);
	//		Gadgetron::clear(res);

	//		double stepsize = (end - start) / (num - 1);
	//		for (int s = 0; s < num; s++)
	//		{
	//			res(s) = start + stepsize*s;
	//		}

	//	}
	//	catch (...)
	//	{
	//		GADGET_THROW("Exceptions happened in linspace(...) ... ");
	//	}
	//	template EXPORTCMR void linspace(const float start, const float end, const int num, hoNDArray<float> & res);
	//	template EXPORTCMR void linspace(const double start, const double end, const int num, hoNDArray<double> & res);
	//}


	template <typename T>
	void compute_thickening(const hoNDArray<T>& endo_mask, const hoNDArray<T>& epi_mask, const size_t ref_phase, hoNDArray<T>& edge_endo, hoNDArray<T>& edge_epi, hoNDArray<T>& rad_strain)
	{
		try
		{
			size_t RO = endo_mask.get_size(0);
			size_t E1 = endo_mask.get_size(1);
			size_t PHS = endo_mask.get_size(2);

			typedef hoNDArray<T> ArrayType;

			hoNDArray<T> thetas(RO, E1, PHS), centroidE(PHS), centroidR(PHS);
			Gadgetron::clear(centroidE);
			Gadgetron::clear(centroidR);

			hoNDArray<T> counter(PHS);
			Gadgetron::clear(counter);

			rad_strain.create(RO, E1, PHS);
			Gadgetron::clear(rad_strain);

			edge_epi.create(RO, E1, PHS);
			Gadgetron::clear(edge_epi);

			edge_endo.create(RO, E1, PHS);
			Gadgetron::clear(edge_endo);

			int samples = std::max(RO, E1);

			// find centroid
			for (size_t e1 = 0; e1 < E1; e1++)
			{

				for (size_t ro = 0; ro < RO; ro++)
				{
					for (size_t p = 0; p < PHS; p++)
					{
						if ((epi_mask(ro, e1, p) > 0) & (endo_mask(ro, e1, p) > 0))
						{
							centroidR(p) += ro;
							centroidE(p) += e1;
							counter(p) += 1;
						}
					}
				}
			}

			Gadgetron::divide(centroidR, counter, centroidR);
			Gadgetron::divide(centroidE, counter, centroidE);

			hoNDArray<T> dists(RO, E1, PHS);

			for (size_t p = 0; p < PHS; p++)
			{
				for (size_t e1 = 0; e1 < E1; e1++)
				{
					for (size_t ro = 0; ro < RO; ro++)
					{
						double cr = centroidR(p);
						double ce = centroidE(p);
						dists(ro, e1, p) = std::sqrt(std::pow(ro - cr, 2) + std::pow(e1 - ce, 2));
					}
				}
			}

			for (size_t p = 0; p < PHS; p++)
			{
				ArrayType endo_mask_2D(RO, E1, const_cast<T*>(&endo_mask(0, 0, p)));
				ArrayType epi_mask_2D(RO, E1, const_cast<T*>(&epi_mask(0, 0, p)));

				hoNDBoundaryHandlerBorderValue<ArrayType> bv_epi(epi_mask_2D);
				hoNDInterpolatorLinear<ArrayType> interp_epi(epi_mask_2D, bv_epi);

				hoNDBoundaryHandlerBorderValue<ArrayType> bv_endo(endo_mask_2D);
				hoNDInterpolatorLinear<ArrayType> interp_endo(endo_mask_2D, bv_endo);

				for (size_t e1 = 0; e1 < E1; e1++)
				{
					for (size_t ro = 0; ro < RO; ro++)
					{
						// compute strain for this point
						double x = ro - (double)centroidR(p);
						double y = (double)centroidE(p) - e1;

						double theta = atan(y / (x + FLT_EPSILON)) + M_PI * (x < 0) + M_PI * 2 * (x >= 0) * (y < 0);
						thetas(ro, e1, p) = theta;

						double x_in = ro + 1.0 * cos(theta);
						double x_out = ro - 1.0 * cos(theta);
						double y_in = e1 - 1.0 * sin(theta);
						double y_out = e1 + 1.0 * sin(theta);

						double epi_in = interp_epi(x_in, y_in);
						double epi_out = interp_epi(x_out, y_out);
						double endo_in = interp_endo(x_in, y_in);
						double endo_out = interp_endo(x_out, y_out);

						if (abs(endo_out - endo_in) > 0.8)
						{
							edge_endo(ro, e1, p) = 1;
						}
						if (abs(epi_out - epi_in) > 0.8)
						{
							edge_epi(ro, e1, p) = 1;
						}
					}
				}
			}

			std::cout << Gadgetron::nrm2(edge_epi) << std::endl;
			std::cout << Gadgetron::nrm2(edge_endo) << std::endl;
			int e, r, p;

//#pragma omp parallel for default(none) private(p, e, r) shared(PHS, RO, E1, rad_strain, centroidE, centroidR, thetas, samples, edge_endo, edge_epi)

			for (e = 0; e < E1; e++)
			{
				for (r = 0; r < RO; r++)
				{
					if (edge_epi(r, e, ref_phase) == 1)
					{
						double cr = centroidR(ref_phase);
						double ce = centroidE(ref_phase);
						double theta_pt = thetas(r, e, ref_phase);
						double check_e_ref = e - samples / 8 * sin(theta_pt);
						double check_r_ref = r + samples / 8 * cos(theta_pt);
						
						hoNDArray<T> test_r_ref(samples), test_e_ref(samples);
						double r_stepsize_ref = (check_r_ref - centroidR(ref_phase))/samples;
						double e_stepsize_ref = (check_e_ref - centroidE(ref_phase))/samples;
						for (int s = 0; s < samples; s++)
						{
							test_r_ref(s) = centroidR(ref_phase) + r_stepsize_ref*s;
							test_e_ref(s) = centroidE(ref_phase) + e_stepsize_ref*s;
						}

						ArrayType edge_endo_ref(RO, E1, const_cast<T*>(&edge_endo(0, 0, ref_phase)));
						ArrayType edge_epi_ref(RO, E1, const_cast<T*>(&edge_epi(0, 0, ref_phase)));

						hoNDBoundaryHandlerBorderValue<ArrayType> bv_edge_endo(edge_endo_ref);
						hoNDInterpolatorLinear<ArrayType> interp_endo_ref(edge_endo_ref, bv_edge_endo);

						hoNDBoundaryHandlerBorderValue<ArrayType> bv_edge_epi(edge_epi_ref);
						hoNDInterpolatorLinear<ArrayType> interp_epi_ref(edge_epi_ref, bv_edge_epi);
						
						hoNDArray<T> endo_ref_line(samples), epi_ref_line(samples);
						for (int s = 0; s < samples; s++)
						{
							endo_ref_line(s) = interp_endo_ref(test_r_ref(s), test_e_ref(s));
							epi_ref_line(s) = interp_epi_ref(test_r_ref(s), test_e_ref(s));
						}
						int endo_ind_ref = Gadgetron::amax(endo_ref_line);
						int epi_ind_ref = Gadgetron::amax(epi_ref_line);

						double endo_e_ref = (test_e_ref(endo_ind_ref));
						double endo_r_ref = (test_r_ref(endo_ind_ref));
						double epi_e_ref = (test_e_ref(epi_ind_ref));
						double epi_r_ref = (test_r_ref(epi_ind_ref));
						double myo_dist_ref = std::sqrt(std::pow(epi_e_ref - endo_e_ref, 2) + std::pow(epi_r_ref - endo_r_ref, 2));

						for (p = 0; p < PHS; p++)
						{
							double check_e = centroidE(p) - centroidE(ref_phase) + check_e_ref;
							double check_r = centroidR(p) - centroidR(ref_phase) + check_r_ref;
							hoNDArray<T> test_r(samples), test_e(samples);
							double r_stepsize = (check_r - centroidR(p)) / samples;
							double e_stepsize = (check_e - centroidE(p)) / samples;
							for (int s = 0; s < samples; s++)
							{
								test_r(s) = centroidR(p) + r_stepsize * s;
								test_e(s) = centroidE(p) + e_stepsize * s;
							}

							ArrayType edge_endo_2D(RO, E1, const_cast<T*>(&edge_endo(0, 0, p)));
							ArrayType edge_epi_2D(RO, E1, const_cast<T*>(&edge_epi(0, 0, p)));

							hoNDBoundaryHandlerBorderValue<ArrayType> bv_endo(edge_endo_2D);
							hoNDInterpolatorLinear<ArrayType> interp_endo(edge_endo_2D, bv_endo);

							hoNDBoundaryHandlerBorderValue<ArrayType> bv_epi(edge_epi_2D);
							hoNDInterpolatorLinear<ArrayType> interp_epi(edge_epi_2D, bv_epi);

							hoNDArray<T> endo_line(samples), epi_line(samples);
							for (int s = 0; s < samples; s++)
							{
								endo_line(s) = interp_endo(test_r(s), test_e(s));
								epi_line(s) = interp_epi(test_r(s), test_e(s));
							}
							int endo_ind = Gadgetron::amax(endo_line);
							int epi_ind = Gadgetron::amax(epi_line);

							double endo_e = test_e(endo_ind);
							double endo_r = test_r(endo_ind);
							double epi_e = test_e(epi_ind);
							double epi_r = test_r(epi_ind);
							double myo_dist = sqrt(std::pow(epi_e - endo_e, 2) + std::pow(epi_r - endo_r, 2));

							rad_strain(r, e, p) = (myo_dist - myo_dist_ref) / myo_dist_ref;
						}
					}
				}
			}
		}
		catch (...)
		{
			GADGET_THROW("Exceptions happened in compute_strain(...) ... ");
		}
	}

	template EXPORTCMR void compute_thickening(const hoNDArray<float>& endo_mask, const hoNDArray<float>& epi_mask, const size_t ref_phase, hoNDArray<float>& edge_endo, hoNDArray<float>& edge_epi, hoNDArray<float>& rad_strains);
	template EXPORTCMR void compute_thickening(const hoNDArray<double>& endo_mask, const hoNDArray<double>& epi_mask, const size_t ref_phase, hoNDArray<double>& edge_endo, hoNDArray<double>& edge_epi, hoNDArray<double>& rad_strain);
}
