#ifndef MRI_CORE_DATA_H
#define MRI_CORE_DATA_H

#include "ismrmrd/ismrmrd.h"
#include "ismrmrd/meta.h"
#include <vector>
#include <set>
#include "hoNDArray.h"
#include <boost/optional.hpp>

namespace Gadgetron 
{

    /** 
      This is a list of lables of the coordinates described in the ISMRMRD acquisition header.

      It is useful for accumulators and triggers and for labeling the storage used in
      the @IsmrmrdAcquisitionBucket and @IsmrmrdDataBuffered structures. 

   */
    enum IsmrmrdCONDITION {
	KSPACE_ENCODE_STEP_1,
	KSPACE_ENCODE_STEP_2,
	AVERAGE,
	SLICE,
	CONTRAST,
	PHASE,
	REPETITION,
	SET,
	SEGMENT,
	USER_0,
	USER_1,
	USER_2,
	USER_3,
	USER_4,
	USER_5,
	USER_6,
	USER_7,
    N_ACQUISITIONS,
	NONE
      };

    // --------------------------------------------------------------------------
    /// define the calibration mode of ISMRMRD
    // --------------------------------------------------------------------------
    enum ismrmrdCALIBMODE
    {
        ISMRMRD_embedded,
        ISMRMRD_interleaved,
        ISMRMRD_separate,
        ISMRMRD_external,
        ISMRMRD_other,
        ISMRMRD_noacceleration
    };

  class SamplingLimit
  {
  public:
    uint16_t min_;
    uint16_t center_;
    uint16_t max_;

    SamplingLimit()
    {
        min_ = 0;
        center_ = 0;
        max_ = 0;
    }
  };
  
  class SamplingDescription
  {
  public:
    // encoding FOV
    float encoded_FOV_[3];
    // recon FOV
    float recon_FOV_[3];
    
    uint16_t encoded_matrix_[3];
    uint16_t recon_matrix_[3];
    
    // sampled range along RO, E1, E2 (for asymmetric echo and partial fourier)
    // min, max and center
    SamplingLimit sampling_limits_[3];

    SamplingDescription()
    {
        encoded_FOV_[0] = 0;
        encoded_FOV_[1] = 0;
        encoded_FOV_[2] = 0;

        recon_FOV_[0] = 0;
        recon_FOV_[1] = 0;
        recon_FOV_[2] = 0;

        encoded_matrix_[0] = 0;
        encoded_matrix_[1] = 0;
        encoded_matrix_[2] = 0;

        recon_matrix_[0] = 0;
        recon_matrix_[1] = 0;
        recon_matrix_[2] = 0;
    }

  };
  
  struct IsmrmrdDataBuffered
  {
  public:
    //7D, fixed order [E0, E1, E2, CHA, N, S, LOC]
    hoNDArray< std::complex<float> > data_;
    
    //7D, fixed order [TRAJ, E0, E1, E2, N, S, LOC]
    boost::optional<hoNDArray<float>> trajectory_;
    
    //5D, fixed order [E1, E2, N, S, LOC]
    hoNDArray< ISMRMRD::AcquisitionHeader > headers_;
    
    SamplingDescription sampling_;

    // function to check if it's empty
  };
  

  /**
     This class is used to group a sub-unit of the data that would feed into a reconstruction. 
   */
  struct IsmrmrdReconBit
  {
  public:
    IsmrmrdDataBuffered data_;
    boost::optional<IsmrmrdDataBuffered> ref_;
  };

  /**
     This class is used to store a unit of data that would feed into a reconstruction. 
   */
  struct IsmrmrdReconData
  {
  public:
    std::vector<IsmrmrdReconBit> rbit_;
  };

  
  /**
     This class is used to store an array of reconstructed data. 
   */
  class IsmrmrdImageArray
  {
  public:
    //7D, fixed order [X, Y, Z, CHA, N, S, LOC]
    hoNDArray< std::complex<float> > data_;
    
    //3D, fixed order [N, S, LOC]
    hoNDArray< ISMRMRD::ImageHeader > headers_;
    
    //3D, fixed order [N, S, LOC]
    //This element is optional (length is 0 if not present)
    std::vector< ISMRMRD::MetaContainer > meta_;
    
  };

}
#endif //MRI_CORE_DATA_H
