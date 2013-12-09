#include "SpiralToGenericGadget.h"
#include "GadgetIsmrmrdReadWrite.h"
#include "vds.h"

#include <algorithm>
#include <vector>

namespace Gadgetron{

  SpiralToGenericGadget::SpiralToGenericGadget()
    : samples_to_skip_start_(0)
    , samples_to_skip_end_(0)
    , samples_per_interleave_(0)
    , prepared_(false)
  {
  }

  SpiralToGenericGadget::~SpiralToGenericGadget() {}

  int SpiralToGenericGadget::process_config(ACE_Message_Block* mb)
  {
    // Start parsing the ISMRMRD XML header
    //

    boost::shared_ptr<ISMRMRD::ismrmrdHeader> cfg = parseIsmrmrdXMLHeader(std::string(mb->rd_ptr()));

    if( cfg.get() == 0x0 ){
      GADGET_DEBUG1("Unable to parse Ismrmrd header\n");
      return GADGET_FAIL;
    }

    ISMRMRD::ismrmrdHeader::encoding_sequence e_seq = cfg->encoding();

    if (e_seq.size() != 1) {
      GADGET_DEBUG2("Number of encoding spaces: %d\n", e_seq.size());
      GADGET_DEBUG1("This Gadget only supports one encoding space\n");
      return GADGET_FAIL;
    }

    //ISMRMRD::encodingSpaceType e_space = (*e_seq.begin()).encodedSpace();
    ISMRMRD::encodingSpaceType r_space = (*e_seq.begin()).reconSpace();
    ISMRMRD::encodingLimitsType e_limits = (*e_seq.begin()).encodingLimits();

    //
    // Setup the spiral trajectory
    //

    if (!(*e_seq.begin()).trajectoryDescription().present()) {
      GADGET_DEBUG1("Trajectory description needed to calculate trajectory");
      return GADGET_FAIL;
    }

    ISMRMRD::trajectoryDescriptionType traj_desc = (*e_seq.begin()).trajectoryDescription().get();

    if (std::strcmp(traj_desc.identifier().c_str(), "HargreavesVDS2000")) {
      GADGET_DEBUG1("Expected trajectory description identifier 'HargreavesVDS2000', not found.");
      return GADGET_FAIL;
    }

    long interleaves = -1;
    long fov_coefficients = -1;
    long sampling_time_ns = -1;
    double max_grad = -1.0;
    double max_slew = -1.0;
    double fov_coeff = -1.0;
    double kr_max = -1.0;

    for (ISMRMRD::trajectoryDescriptionType::userParameterLong_sequence::iterator i (traj_desc.userParameterLong().begin ()); i != traj_desc.userParameterLong().end(); ++i) {
      if (std::strcmp(i->name().c_str(),"interleaves") == 0) {
	interleaves = i->value();
      } else if (std::strcmp(i->name().c_str(),"fov_coefficients") == 0) {
	fov_coefficients = i->value();
      } else if (std::strcmp(i->name().c_str(),"SamplingTime_ns") == 0) {
	sampling_time_ns = i->value();
      } else {
	GADGET_DEBUG2("WARNING: unused trajectory parameter %s found\n", i->name().c_str());
      }
    }

    for (ISMRMRD::trajectoryDescriptionType::userParameterDouble_sequence::iterator i (traj_desc.userParameterDouble().begin ()); i != traj_desc.userParameterDouble().end(); ++i) {
      if (std::strcmp(i->name().c_str(),"MaxGradient_G_per_cm") == 0) {
	max_grad = i->value();
      } else if (std::strcmp(i->name().c_str(),"MaxSlewRate_G_per_cm_per_s") == 0) {
	max_slew = i->value();
      } else if (std::strcmp(i->name().c_str(),"FOVCoeff_1_cm") == 0) {
	fov_coeff = i->value();
      } else if (std::strcmp(i->name().c_str(),"krmax_per_cm") == 0) {
	kr_max= i->value();
      } else {
	GADGET_DEBUG2("WARNING: unused trajectory parameter %s found\n", i->name().c_str());
      }
    }

    if ((interleaves < 0) || (fov_coefficients < 0) || (sampling_time_ns < 0) || (max_grad < 0) || (max_slew < 0) || (fov_coeff < 0) || (kr_max < 0)) {
      GADGET_DEBUG1("Appropriate parameters for calculating spiral trajectory not found in XML configuration\n");
      return GADGET_FAIL;
    }

    Tsamp_ns_ = sampling_time_ns;
    Nints_ = interleaves;
    interleaves_ = static_cast<int>(Nints_);

    gmax_ = max_grad;
    smax_ = max_slew;
    krmax_ = kr_max;
    fov_ = fov_coeff;

    samples_to_skip_start_  =  0; //n.get<int>(std::string("samplestoskipstart.value"))[0];
    samples_to_skip_end_    = -1; //n.get<int>(std::string("samplestoskipend.value"))[0];

    GADGET_DEBUG2("smax:                    %f\n", smax_);
    GADGET_DEBUG2("gmax:                    %f\n", gmax_);
    GADGET_DEBUG2("Tsamp_ns:                %d\n", Tsamp_ns_);
    GADGET_DEBUG2("Nints:                   %d\n", Nints_);
    GADGET_DEBUG2("fov:                     %f\n", fov_);
    GADGET_DEBUG2("krmax:                   %f\n", krmax_);
    GADGET_DEBUG2("samples_to_skip_start_ : %d\n", samples_to_skip_start_);
    GADGET_DEBUG2("samples_to_skip_end_   : %d\n", samples_to_skip_end_);

    return GADGET_OK;
  }

  int SpiralToGenericGadget::
  process(GadgetContainerMessage<ISMRMRD::AcquisitionHeader> *m1,
	  GadgetContainerMessage< hoNDArray< std::complex<float> > > *m2)
  {
    // Noise should have been consumed by the noise adjust, but just in case...
    //

    bool is_noise = ISMRMRD::FlagBit(ISMRMRD::ACQ_IS_NOISE_MEASUREMENT).isSet(m1->getObjectPtr()->flags);
    if (is_noise) {
      m1->release();
      return GADGET_OK;
    }

    // Compute hoNDArray of trajectory and weights at first pass
    //

    if (!prepared_) {

      int     nfov   = 1;         /*  number of fov coefficients.             */
      int     ngmax  = 1e5;       /*  maximum number of gradient samples      */
      double  *xgrad;             /*  x-component of gradient.                */
      double  *ygrad;             /*  y-component of gradient.                */
      double  *x_trajectory;
      double  *y_trajectory;
      double  *weighting;
      int     ngrad;
      double sample_time = (1.0*Tsamp_ns_) * 1e-9;

      // Calculate gradients 
      calc_vds(smax_,gmax_,sample_time,sample_time,Nints_,&fov_,nfov,krmax_,ngmax,&xgrad,&ygrad,&ngrad);

      samples_per_interleave_ = std::min(ngrad,static_cast<int>(m1->getObjectPtr()->number_of_samples));
      GADGET_DEBUG2("Using %d samples per interleave\n", samples_per_interleave_);

      // Calculate the trajectory and weights
      calc_traj(xgrad, ygrad, samples_per_interleave_, Nints_, sample_time, krmax_, &x_trajectory, &y_trajectory, &weighting);

      std::vector<unsigned long long> trajectory_dimensions;
      trajectory_dimensions.push_back(3);
      trajectory_dimensions.push_back(samples_per_interleave_*Nints_);

      host_traj_ = boost::shared_ptr< hoNDArray<float> >(new hoNDArray<float>(&trajectory_dimensions));

      {
	float* co_ptr = reinterpret_cast<float*>(host_traj_->get_data_ptr());
	
	for (int i = 0; i < (samples_per_interleave_*Nints_); i++) {
	  co_ptr[i*3+0] = -x_trajectory[i]/2;
	  co_ptr[i*3+1] = -y_trajectory[i]/2;
	  co_ptr[i*3+2] = weighting[i];
	}
      }

      delete [] xgrad;
      delete [] ygrad;
      delete [] x_trajectory;
      delete [] y_trajectory;
      delete [] weighting;

      prepared_ = true;
    }

    // Adjustments based in the incoming data
    //

    if (samples_to_skip_end_ == -1) {
      samples_to_skip_end_ = m1->getObjectPtr()->number_of_samples-samples_per_interleave_;
      GADGET_DEBUG2("Adjusting samples_to_skip_end_ = %d\n", samples_to_skip_end_);
    }

    // Define some utility variables
    //

    unsigned int samples_to_copy = m1->getObjectPtr()->number_of_samples-samples_to_skip_end_;
    unsigned int interleave = m1->getObjectPtr()->idx.kspace_encode_step_1;

    // Prepare for a new array continuation for the trajectory/weights of the incoming profile
    //

    std::vector<unsigned long long> trajectory_dimensions;
    trajectory_dimensions.push_back(3);
    trajectory_dimensions.push_back(samples_per_interleave_);
    
    hoNDArray<float> *traj_source = new hoNDArray<float>
      (&trajectory_dimensions, host_traj_->get_data_ptr()+3*samples_per_interleave_*interleave);
    
    // Make a new array as continuation of m1, and pass along
    //

    GadgetContainerMessage< hoNDArray<float> > *cont = new GadgetContainerMessage< hoNDArray<float> >();
    *(cont->getObjectPtr()) = *traj_source;
    m2->cont(cont);
    
    if (this->next()->putq(m1) < 0) {
      GADGET_DEBUG1("Failed to put job on queue.\n");
      return GADGET_FAIL;
    }
    
    return GADGET_OK;
  }
  
  GADGET_FACTORY_DECLARE(SpiralToGenericGadget)
}
