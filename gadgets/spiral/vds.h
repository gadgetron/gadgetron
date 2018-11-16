#pragma once

#include "gadgetron_spiral_export.h"

namespace Gadgetron{

  void EXPORTGADGETS_SPIRAL 
  calc_vds(double slewmax,double gradmax,double Tgsample,double Tdsample,int Ninterleaves,
	   double* fov, int numfov,double krmax,
	   int ngmax, double** xgrad,double** ygrad,int* numgrad);
  
  void EXPORTGADGETS_SPIRAL 
  calc_traj(double* xgrad, double* ygrad, int ngrad, int Nints, double Tgsamp, double krmax,
	    double** x_trajectory, double** y_trajectory, double** weights);  
}
