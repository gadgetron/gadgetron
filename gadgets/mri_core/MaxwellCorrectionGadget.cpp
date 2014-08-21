#include "MaxwellCorrectionGadget.h"
#include "Gadgetron.h"
#include "GadgetronTimer.h"
#include "Spline.h"
#include "ismrmrd_xml.h"

#include <numeric>
#ifdef USE_OMP
#include <omp.h>
#endif 

namespace Gadgetron{

  #ifdef M_PI
    #undef M_PI
  #endif // M_PI
  #define M_PI 3.14159265358979323846

  MaxwellCorrectionGadget::MaxwellCorrectionGadget()
    : maxwell_coefficients_present_(false)
    , maxwell_coefficients_(4,0)
  {
  }

  MaxwellCorrectionGadget::~MaxwellCorrectionGadget() {}

  int MaxwellCorrectionGadget::process_config(ACE_Message_Block* mb)
  {

    ISMRMRD::IsmrmrdHeader h;
    ISMRMRD::deserialize(mb->rd_ptr(),h);

    if (h.userParameters) {
      for (std::vector<ISMRMRD::UserParameterDouble>::const_iterator i (h.userParameters->userParameterDouble.begin()); 
	   i != h.userParameters->userParameterDouble.end(); i++) {
	if (std::strcmp(i->name.c_str(),"MaxwellCoefficient_0") == 0) {
	  maxwell_coefficients_[0] = i->value;
	} else if (std::strcmp(i->name.c_str(),"MaxwellCoefficient_1") == 0) {
	  maxwell_coefficients_[1] = i->value;
	} else if (std::strcmp(i->name.c_str(),"MaxwellCoefficient_2") == 0) {
	  maxwell_coefficients_[2] = i->value;
	} else if (std::strcmp(i->name.c_str(),"MaxwellCoefficient_3") == 0) {
	  maxwell_coefficients_[3] = i->value;
	} else {
	  GADGET_DEBUG2("WARNING: unused user parameter parameter %s found\n", i->name.c_str());
	}
      }
    } else {
      GADGET_DEBUG1("MaxwellCorrection coefficients are supposed to be in the UserParameters. No user parameter section found\n");
      return GADGET_OK;
    }

    maxwell_coefficients_present_ = true;

    GADGET_DEBUG2("Maxwell Coefficients: %f, %f, %f, %f\n", maxwell_coefficients_[0], maxwell_coefficients_[1], maxwell_coefficients_[2], maxwell_coefficients_[3]);

    return GADGET_OK;
  }

  int MaxwellCorrectionGadget::
  process(GadgetContainerMessage<ISMRMRD::ImageHeader>* m1,
	  GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2)
  {
    if (maxwell_coefficients_present_) {
      //GADGET_DEBUG1("Got coefficients\n");

      int Nx = m2->getObjectPtr()->get_size(0);
      int Ny = m2->getObjectPtr()->get_size(1);
      int Nz = m2->getObjectPtr()->get_size(2);

      float dx = m1->getObjectPtr()->field_of_view[0] / Nx;
      float dy = m1->getObjectPtr()->field_of_view[1] / Ny;
      float dz = m1->getObjectPtr()->field_of_view[2] / Nz;

      /*
      GADGET_DEBUG2("Nx = %d, Ny = %d, Nz = %d\n", Nx, Ny, Nz);
      GADGET_DEBUG2("dx = %f, dy = %f, dz = %f\n", dx, dy, dz);
      GADGET_DEBUG2("img_pos_x = %f, img_pos_y = %f, img_pos_z = %f\n", m1->getObjectPtr()->position[0], m1->getObjectPtr()->position[1], m1->getObjectPtr()->position[2]);
      */

      std::vector<float> dR(3,0);
      std::vector<float> dP(3,0);
      std::vector<float> dS(3,0);
      std::vector<float> p(3,0);

      for (int z = 0; z < Nz; z++) {
	for (int y = 0; y < Ny; y++) {
	  for (int x = 0; x < Nx; x++) {
	   
	    dR[0] = (x-Nx/2+0.5) * dx * m1->getObjectPtr()->read_dir[0];
	    dR[1] = (x-Nx/2+0.5) * dx * m1->getObjectPtr()->read_dir[1];
	    dR[2] = (x-Nx/2+0.5) * dx * m1->getObjectPtr()->read_dir[2];
	    
	    dP[0] = (y-Ny/2+0.5) * dy * m1->getObjectPtr()->phase_dir[0];
	    dP[1] = (y-Ny/2+0.5) * dy * m1->getObjectPtr()->phase_dir[1];
	    dP[2] = (y-Ny/2+0.5) * dy * m1->getObjectPtr()->phase_dir[2];
	    
	    if (Nz > 1) {
	      dS[0] = (z-Nz/2+0.5) * dz * m1->getObjectPtr()->slice_dir[0];
	      dS[1] = (z-Nz/2+0.5) * dz * m1->getObjectPtr()->slice_dir[1];
	      dS[2] = (z-Nz/2+0.5) * dz * m1->getObjectPtr()->slice_dir[2];
	    }

	    p[0] = m1->getObjectPtr()->position[0] + dP[0] + dR[0] + dS[0];
	    p[1] = m1->getObjectPtr()->position[1] + dP[1] + dR[1] + dS[1];
	    p[2] = m1->getObjectPtr()->position[2] + dP[1] + dR[2] + dS[2];

	    //Convert to centimeters
	    p[0] = p[0]/1000.0;
	    p[1] = p[1]/1000.0;
	    p[2] = p[2]/1000.0;

	    float delta_phi = maxwell_coefficients_[0]*p[2]*p[2] +
	      maxwell_coefficients_[1]*(p[0]*p[0] + p[1]*p[1]) + 
	      maxwell_coefficients_[2]*p[0]*p[2] + 
	      maxwell_coefficients_[3]*p[1]*p[2];

	    long index = z*Ny*Nx+y*Nx+x;
	    std::complex<float>* data_ptr = m2->getObjectPtr()->get_data_ptr();

	    std::complex<float> correction = std::polar(1.0f,static_cast<float>(2*M_PI*delta_phi));

	    //data_ptr[index] *= correction;
	  }
	}
      }

    }

    if (this->next()->putq(m1) < 0) {
      GADGET_DEBUG1("Unable to put data on next Gadgets Q\n");
      return GADGET_FAIL;
    }
    return GADGET_OK;
  }

  GADGET_FACTORY_DECLARE(MaxwellCorrectionGadget)
}
