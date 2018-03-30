#ifndef GADGETMRIHEADERS_H
#define GADGETMRIHEADERS_H

#include <ace/Basic_Types.h>

//Data flags
/*
#define GADGET_FLAG_ACQ_END                   (1 << 0)
#define GADGET_FLAG_LAST_ACQ_IN_SLICE         (1 << 1)
#define GADGET_FLAG_LAST_ACQ_IN_MEAS          (1 << 2)
#define GADGET_FLAG_LAST_ACQ_IN_CONCAT        (1 << 3)
#define GADGET_FLAG_FIRST_ACQ_IN_SLICE        (1 << 4)
#define GADGET_FLAG_FIRST_ACQ_IN_MEAS         (1 << 5)
#define GADGET_FLAG_FIRST_ACQ_IN_CONCAT       (1 << 6)
#define GADGET_FLAG_IS_NOISE_SCAN             (1 << 7)
#define GADGET_FLAG_IS_PATREF_SCAN            (1 << 8)
#define GADGET_FLAG_IS_PATREFANDIMA_SCAN      (1 << 9)

#define GADGET_FLAG_LAST_IMAGE                (1 << 0)

enum GadgetImageFormats {
	GADGET_IMAGE_COMPLEX_FLOAT = 0,
	GADGET_IMAGE_REAL_FLOAT,
	GADGET_IMAGE_REAL_UNSIGNED_SHORT
};

enum GadgetImageTypes
{
	GADGET_IMAGE_MAGNITUDE = 0,
	GADGET_IMAGE_PHASE,
	GADGET_IMAGE_REAL,
	GADGET_IMAGE_IMAG
};
*/

namespace Gadgetron{

enum GadgetMessageID {
  GADGET_MESSAGE_EXT_ID_MIN                             = 1000,
  GADGET_MESSAGE_ACQUISITION                            = 1001, /**< DEPRECATED */
  GADGET_MESSAGE_NEW_MEASUREMENT                        = 1002, /**< DEPRECATED */
  GADGET_MESSAGE_END_OF_SCAN                            = 1003, /**< DEPRECATED */
  GADGET_MESSAGE_IMAGE_CPLX_FLOAT                       = 1004, /**< DEPRECATED */
  GADGET_MESSAGE_IMAGE_REAL_FLOAT                       = 1005, /**< DEPRECATED */
  GADGET_MESSAGE_IMAGE_REAL_USHORT                      = 1006, /**< DEPRECATED */
  GADGET_MESSAGE_EMPTY                                  = 1007, /**< DEPRECATED */
  GADGET_MESSAGE_ISMRMRD_ACQUISITION                    = 1008,
  GADGET_MESSAGE_ISMRMRD_IMAGE_CPLX_FLOAT               = 1009, /**< DEPRECATED */
  GADGET_MESSAGE_ISMRMRD_IMAGE_REAL_FLOAT               = 1010, /**< DEPRECATED */
  GADGET_MESSAGE_ISMRMRD_IMAGE_REAL_USHORT              = 1011, /**< DEPRECATED */
  GADGET_MESSAGE_DICOM                                  = 1012, /**< DEPRECATED */
  GADGET_MESSAGE_CLOUD_JOB                              = 1013,
  GADGET_MESSAGE_GADGETCLOUD_JOB                        = 1014,
  GADGET_MESSAGE_ISMRMRD_IMAGEWITHATTRIB_CPLX_FLOAT     = 1015, /**< DEPRECATED */
  GADGET_MESSAGE_ISMRMRD_IMAGEWITHATTRIB_REAL_FLOAT     = 1016, /**< DEPRECATED */
  GADGET_MESSAGE_ISMRMRD_IMAGEWITHATTRIB_REAL_USHORT    = 1017, /**< DEPRECATED */
  GADGET_MESSAGE_DICOM_WITHNAME                         = 1018,
  GADGET_MESSAGE_DEPENDENCY_QUERY                       = 1019,
  GADGET_MESSAGE_ISMRMRD_IMAGE_REAL_SHORT               = 1020, /**< DEPRECATED */
  GADGET_MESSAGE_ISMRMRD_IMAGEWITHATTRIB_REAL_SHORT     = 1021, /**< DEPRECATED */
  GADGET_MESSAGE_ISMRMRD_IMAGE                          = 1022,
  GADGET_MESSAGE_RECONDATA                              = 1023,
  GADGET_MESSAGE_ISMRMRD_WAVEFORM                       = 1026,
  GADGET_MESSAGE_EXT_ID_MAX                             = 4096
};


/*
struct ISMRMRD::ImageHeader
{
  ACE_UINT32     flags;
  ACE_UINT16     matrix_size[3];
  ACE_UINT16     channels;
  float          position[3];
  float          quaternion[4];
  float			 table_position;
  ACE_UINT16     slice;
  ACE_UINT16     contrast;
  ACE_UINT16     set;
  ACE_UINT16     phase;
  ACE_UINT16     average;
  ACE_UINT16     repetition;
  ACE_UINT32     time_stamp;
  ACE_UINT32     pmu_time_stamp;
  ACE_UINT16     image_format;
  ACE_UINT16     image_type;
  ACE_UINT16     image_index;
  ACE_UINT16	 image_series_index;

  ACE_UINT16 get_matrix_size(unsigned int index) {
    if (index < 3) {
      return matrix_size[index];
    } else {
      return 0;
    }
  }

  void set_matrix_size(unsigned int index, ACE_UINT16 size) {
    if (index < 3) {
      matrix_size[index] = size;
    }
  }

  float get_position(unsigned int index) {
    if (index < 3) {
      return position[index];
    } else {
      return 0.0f;
    }
  }

  void set_position(unsigned int index, float pos)
  {
    if (index < 3) {
      position[index] = pos;
    }
  }

  float get_quaternion(unsigned int index) {
    if (index < 4) {
      return quaternion[index];
    } else {
      return 0.0f;
    }
  }

  void set_quaternion(unsigned int index, float quar)
  {
    if (index < 4) {
      quaternion[index] = quar;
    }
  }

}; 
*/
}

#endif  //GADGETMRIHEADERS_H
