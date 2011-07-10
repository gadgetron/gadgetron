#ifndef GADGETMRIHEADERS_H
#define GADGETMRIHEADERS_H

#include "gadgetron_export.h"
#include <ace/Basic_Types.h>

//Data flags
#define GADGET_FLAG_ACQ_END                   (1 << 0) /* 0x01 */
#define GADGET_FLAG_LAST_ACQ_IN_SLICE         (1 << 1) /* 0x02 */
#define GADGET_FLAG_LAST_ACQ_IN_MEAS          (1 << 2) /* 0x04 */
#define GADGET_FLAG_LAST_ACQ_IN_CONCAT        (1 << 3) /* 0x08 */
#define GADGET_FLAG_FIRST_ACQ_IN_SLICE        (1 << 4) /* 0x02 */
#define GADGET_FLAG_FIRST_ACQ_IN_MEAS         (1 << 5) /* 0x04 */
#define GADGET_FLAG_FIRST_ACQ_IN_CONCAT       (1 << 6) /* 0x08 */

EXPORTGADGETSCORE enum GadgetMessageID {
  GADGET_MESSAGE_EXT_ID_MIN        = 1000,
  GADGET_MESSAGE_ACQUISITION       = 1001,
  GADGET_MESSAGE_NEW_MEASUREMENT   = 1002,
  GADGET_MESSAGE_END_OF_SCAN       = 1003,
  GADGET_MESSAGE_IMAGE             = 1004,
  GADGET_MESSAGE_EMPTY             = 1005,
  GADGET_MESSAGE_EXT_ID_MAX        = 4096
};
  
struct EXPORTGADGETSCORE LoopCounters {
  ACE_UINT16 line;
  ACE_UINT16 acquisition;
  ACE_UINT16 slice;
  ACE_UINT16 partition;
  ACE_UINT16 echo;
  ACE_UINT16 phase;
  ACE_UINT16 repetition;
  ACE_UINT16 set;
  ACE_UINT16 segment;
  ACE_UINT16 channel;
};

struct EXPORTGADGETSCORE GadgetMessageAcquisition
{
  ACE_UINT32     flags;
  ACE_UINT32     meas_uid;
  ACE_UINT32     scan_counter;
  ACE_UINT32     time_stamp;
  ACE_UINT16     samples;
  ACE_UINT16     channels;
  float          position[3];
  float          quarternion[4];   
  LoopCounters   idx;
  LoopCounters   min_idx;
  LoopCounters   max_idx;
};

struct EXPORTGADGETSCORE GadgetMessageImage
{
  ACE_UINT16     matrix_size[3];
  ACE_UINT16     channels;
  float          position[3];
  float          quarternion[4];
  LoopCounters   data_idx_min;
  LoopCounters   data_idx_max;
  LoopCounters   data_idx_current;
  ACE_UINT32     time_stamp;
}; 

#endif  //GADGETMRIHEADERS_H
