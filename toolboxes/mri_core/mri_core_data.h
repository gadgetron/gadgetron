#ifndef MRI_CORE_DATA_H
#define MRI_CORE_DATA_H

#include "GadgetContainerMessage.h"
#include "ismrmrd.h"

#include <vector>

namespace Gadgetron 
{

  /** 
      This class functions as a storage unit for GadgetContainerMessage pointers
      that point to acquisiton headers, data and trajectories.

      It is the storage used in the @IsmrmrdAcquisitionBucket structure. 

   */
  class IsmrmrdAcquisitionData
  {
  public:
    IsmrmrdAcquisitionData()
      : head_(0)
      , data_(0)
      , traj_(0)
    {

    }

    /**
       Destructor. The memory in the GadgetContainer Messages will be deleted
       when the object is destroyed. 
     */
    ~IsmrmrdAcquisitionData() {
      if (traj_) {
	traj_->release();
	traj_ = 0;
      }
      
      if (data_) {
	data_->release();
	data_ = 0;
      }

      if (head_) {
	head_->release();
	head_ = 0;
      }
    }


    GadgetContainerMessage<ISMRMRD::AcquisitionHeader>* head_;
    GadgetContainerMessage< hoNDArray< std::complex<float> > >* data_;
    GadgetContainerMessage< hoNDArray< float > > * traj_;
  };


  /**

     This class serves as the storage unit for buffered data. 
     The @IsmrmrdAcquisitionData structure contains pointers 
     to the GadgetContainerMessages with the data. 

     Data stored in these buckets will automatically get deleted when the object is
     destroyed. 

   */ 
  class IsmrmrdAcquisitionBucket
  {
    std::vector< IsmrmrdAcquisitionData > data_;
    std::vector< IsmrmrdAcquisitionData > ref_;
  };
  
  
  class SamplingLimit
  {
    uint16_t min_;
    uint16_t max_;
    uint16_t center_;
  };
  
  class SamplingDescription
  {
    
    // encoding FOV
    float encoding_FOV_[3];
    // recon FOV
    float recon_FOV_[3];
    
    uint16_t encoded_matrix_[3];
    uint16_t recon_matrix_[3];
    
    // sampled range along RO, E1, E2 (for asymmetric echo and partial fourier)
    // start and end
    SamplingLimit sampling_limits_[3];
  };
  
  class IsmrmrdDataBuffered
  {
    
    //11D, fixed order [RO, E1, E2, CHA, SLC, N, S]
    hoNDArray< std::complex<float> > data_;
    
    //11D, fixed order [RO, E1, E2, CHA, SLC, PHS, CON, REP, SET, SEG, AVE]
    //This element is optional (length is 0 if not present)
    hoNDArray< float > trajectory_;
    
    //9D, fixed order [E1, E2, SLC, PHS, CON, REP, SET, SEG, AVE]
    hoNDArray< ISMRMRD::AcquisitionHeader > headers_;
    
    SamplingDescription sampling_;
  };
  

  /**
     This class is used to store a unit of data that would feed into a reconstruction. 
   */
  class IsmrmrdReconData
  {
    std::vector<IsmrmrdDataBuffered> data_;
    std::vector<IsmrmrdDataBuffered> ref_;
  };

}
#endif //MRI_CORE_DATA_H
