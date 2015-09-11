#ifndef MRI_CORE_ACQUISITION_BUCKET_H
#define MRI_CORE_ACQUISITION_BUCKET_H

#include "GadgetContainerMessage.h"
#include "mri_core_data.h"

namespace Gadgetron 
{

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
  
}
#endif //MRI_CORE_ACQUISITION_BUCKET_H
