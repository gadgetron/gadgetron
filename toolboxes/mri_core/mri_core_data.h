#ifndef MRI_CORE_DATA_H_
#define MRI_CORE_DATA_H_

#include "GadgetContainerMessage.h"
#include "ismrmrd/ismrmrd.h"
#include "Gadgetron.h"
#include <vector>
#include <set>

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
	NONE
      };
    
  /** 
      This class functions as a storage unit for statistics related to
      the @IsmrmrdAcquisitionData objects.

   */
  class IsmrmrdAcquisitionBucketStats
  {
    public:
      // Set of labels found in the data or ref part of a bucket
      //11D, fixed order [RO, E1, E2, CHA, SLC, PHS, CON, REP, SET, SEG, AVE]
      std::set<uint16_t> kspace_encode_step_1;
      std::set<uint16_t> kspace_encode_step_2;
      std::set<uint16_t> slice;
      std::set<uint16_t> phase;
      std::set<uint16_t> contrast;
      std::set<uint16_t> repetition;
      std::set<uint16_t> set;
      std::set<uint16_t> segment;
      std::set<uint16_t> average;
  };

  /** 
      This class functions as a storage unit for GadgetContainerMessage pointers
      that point to acquisiton headers, data and trajectories.

      It is the storage used in the @IsmrmrdAcquisitionBucket structure. 

   */
  class IsmrmrdAcquisitionData
  {
  public:
    /**
       Default Constructor
    */
    IsmrmrdAcquisitionData()
      : head_(0)
      , data_(0)
      , traj_(0)
      {

      }
    
    /**
       Constructor
    */
    IsmrmrdAcquisitionData(GadgetContainerMessage<ISMRMRD::AcquisitionHeader>* head,
                           GadgetContainerMessage< hoNDArray< std::complex<float> > >* data,
                           GadgetContainerMessage< hoNDArray< float > >* traj = 0)
    {
      if (head) {
	head_ = head->duplicate();
      } else {
	head_ = 0;
      }

      if (data) {
	data_ = data->duplicate();
      } else {
	data_ = 0;
      }

      if (traj) {
	traj_ = traj->duplicate();
      } else {
	traj_ = 0;
      }
    }

    /** 
	Assignment operator
     */
    IsmrmrdAcquisitionData& operator=(const IsmrmrdAcquisitionData& d)
      {
	if (this != &d) { 
	  if (d.head_) {
	    if (head_) head_->release();
	    head_ = d.head_->duplicate();
	  } else {
	    head_ = 0;
	  }
	  
	  if (d.data_) {
	    if (data_) data_->release();
	    data_ = d.data_->duplicate();
	  } else {
	    data_ = 0;
	  }
	  
	  if (d.traj_) {
	    if (traj_) traj_->release();
	    traj_ = d.traj_->duplicate();
	  } else {
	    traj_ = 0;
	  }
	}
	return *this;
      }

    /**
       Copy constructor
     */
    IsmrmrdAcquisitionData(const IsmrmrdAcquisitionData& d)
      : head_(0)
      , data_(0)
      , traj_(0)
      {
	*this = d;
      }


    /**
       Destructor. The memory in the GadgetContainer Messages will be deleted
       when the object is destroyed. 
     */
    ~IsmrmrdAcquisitionData() {
      if (head_) {
	head_->release();
	head_ = 0;
      }

      if (data_) {
	data_->release();
	data_ = 0;
      }

      if (traj_) {
	traj_->release();
	traj_ = 0;
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
  public:
    std::vector< IsmrmrdAcquisitionData > data_;
    std::vector< IsmrmrdAcquisitionData > ref_;
    std::vector< IsmrmrdAcquisitionBucketStats > datastats_;
    std::vector< IsmrmrdAcquisitionBucketStats > refstats_;
  };
  
  
  class SamplingLimit
  {
  public:
    uint16_t min_;
    uint16_t max_;
    uint16_t center_;
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
  };
  
  class IsmrmrdDataBuffered
  {
  public:
    //7D, fixed order [RO, E1, E2, CHA, SLC, N, S]
    hoNDArray< std::complex<float> > data_;
    
    //11D, fixed order [RO, E1, E2, CHA, SLC, PHS, CON, REP, SET, SEG, AVE]
    //This element is optional (length is 0 if not present)
    hoNDArray< float > trajectory_;
    
    //9D, fixed order [E1, E2, SLC, PHS, CON, REP, SET, SEG, AVE]
    hoNDArray< ISMRMRD::AcquisitionHeader > headers_;
    
    SamplingDescription sampling_;

    // function to check if it's empty
  };
  

  /**
     This class is used to group a sub-unit of the data that would feed into a reconstruction. 
   */
  class IsmrmrdReconBit
  {
  public:
    IsmrmrdDataBuffered data_;
    IsmrmrdDataBuffered ref_;
  };

  /**
     This class is used to store a unit of data that would feed into a reconstruction. 
   */
  class IsmrmrdReconData
  {
  public:
    std::vector<IsmrmrdReconBit> rbit_;
  };
  
}
#endif //MRI_CORE_DATA_H_
