/** \file   cmr_radial_thickening.h
	\brief  Implement functionalities to handle cardiac time stamps
	\author Angela Gao
	\date   July 9, 2019
*/

#include "cmr_radial_thickening.h"
#include "hoNDArray_reductions.h"
#include "hoNDArray_elemwise.h"
#include "hoNDInterpolator.h"
#include "hoNDBoundaryHandler.h"
#define _USE_MATH_DEFINES
#include <math.h>

namespace Gadgetron {

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

						double x_in = ro + 0.5 * cos(theta);
						double x_out = ro - 0.5 * cos(theta);
						double y_in = e1 - 0.5 * sin(theta);
						double y_out = e1 + 0.5 * sin(theta);

						double epi_in = interp_epi(x_in, y_in);
						double epi_out = interp_epi(x_out, y_out);
						double endo_in = interp_endo(x_in, y_in);
						double endo_out = interp_endo(x_out, y_out);

						if (abs(endo_out - endo_in) > 0.35)
						{
							edge_endo(ro, e1, p) = 1;
						}
						if (abs(epi_out - epi_in) > 0.35)
						{
							edge_epi(ro, e1, p) = 1;
						}
					}
				}
			}
			//std::cout << Gadgetron::nrm2(edge_epi) << std::endl;
			//std::cout << Gadgetron::nrm2(edge_endo) << std::endl;
			int epi_e_ref, epi_r_ref, p;

//#pragma omp parallel for default(none) private(p, epi_e_ref, epi_r_ref) shared(PHS, RO, E1, rad_strain, mask, centroidE, centroidR, thetas)

			for (epi_e_ref = 0; epi_e_ref < E1; epi_e_ref++)
			{
				for (epi_r_ref = 0; epi_r_ref < RO; epi_r_ref++)
				{
					if (edge_epi(epi_r_ref, epi_e_ref, ref_phase) == 1)
					{
						double cr = centroidR(ref_phase);
						double ce = centroidE(ref_phase);
						double theta_pt = thetas(epi_r_ref, epi_e_ref, ref_phase);
						double check_e = epi_e_ref - samples / 8 * sin(theta_pt);
						double check_r = epi_r_ref + samples / 8 * cos(theta_pt);
						
						hoNDArray<T> test_r(samples), test_e(samples);
						double r_stepsize = (check_r - centroidR(ref_phase))/samples;
						double e_stepsize = (check_e - centroidE(ref_phase))/samples;
						for (int s = 0; s < samples; s++)
						{
							test_r(s) = centroidR(ref_phase) + r_stepsize*s;
							test_e(s) = centroidE(ref_phase) + e_stepsize*s;
						}

						ArrayType edge_endo_ref(RO, E1, const_cast<T*>(&edge_endo(0, 0, ref_phase)));

						hoNDBoundaryHandlerBorderValue<ArrayType> bv_edge_endo(edge_endo_ref);
						hoNDInterpolatorLinear<ArrayType> interp_ref(edge_endo_ref, bv_edge_endo);
						
						hoNDArray<T> endo_ref_line(samples);
						for (int s = 0; s < samples; s++)
						{
							endo_ref_line(s) = interp_ref(test_r(s), test_e(s));
						}
						int endo_ind_ref = Gadgetron::amax(endo_ref_line);
						double endo_e_ref = (test_e(endo_ind_ref));
						double endo_r_ref = (test_r(endo_ind_ref));
						double myo_dist_ref = std::sqrt(std::pow(epi_e_ref - endo_e_ref, 2) + std::pow(epi_r_ref - endo_r_ref, 2));

						for (p = 0; p < PHS; p++)
						{

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

							rad_strain(epi_r_ref, epi_e_ref, p) = (myo_dist - myo_dist_ref) / myo_dist_ref;
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
