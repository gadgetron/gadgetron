#include <stdint.h>
#include <iostream>
#include <fstream>
#include <map>

#define MDH_NUMBEROFICEPROGRAMPARA 4
#define MDH_FREEHDRPARA            4

enum Trajectory {
  TRAJECTORY_CARTESIAN = 0x01,
  TRAJECTORY_RADIAL    = 0x02,
  TRAJECTORY_SPIRAL    = 0x04,
  TRAJECTORY_BLADE     = 0x08
};

struct mdhLC {
  uint16_t ushLine;
  uint16_t ushAcquisition;
  uint16_t ushSlice;
  uint16_t ushPartition;
  uint16_t ushEcho;
  uint16_t ushPhase;
  uint16_t ushRepetition;
  uint16_t ushSet;
  uint16_t ushSeg;
  uint16_t ushIda;
  uint16_t ushIdb;
  uint16_t ushIdc;
  uint16_t ushIdd;
  uint16_t ushIde;
};

struct mdhCutOff {
  uint16_t ushPre;
  uint16_t ushPost;
};

struct mdhSlicePosVec {
  float flSag;
  float flCor;
  float flTra;
};

struct mdhSD {
  mdhSlicePosVec sSlicePosVec;
};

struct mdhSliceData {
  float aflQuaternion[4];
};


struct sMDH
{
  uint32_t ulFlagsAndDMALength;
  int32_t lMeasUID;
  uint32_t ulScanCounter;
  uint32_t ulTimeStamp;
  uint32_t ulPMUTimeStamp;
  uint32_t aulEvalInfoMask[2];
  uint16_t ushSamplesInScan;
  uint16_t ushUsedChannels;
  mdhLC sLC;
  mdhCutOff sCutOff;

  uint16_t ushKSpaceCentreColumn;
  uint16_t ushCoilSelect;
  float fReadOutOffcentre;
  uint32_t ulTimeSinceLastRF;
  uint16_t ushKSpaceCentreLineNo;
  uint16_t ushKSpaceCentrePartitionNo;
  uint16_t aushIceProgramPara[4];
  uint16_t aushFreePara[4];

  mdhSD sSD;
  mdhSliceData sSliceData;

  uint16_t ushChannelId;
  uint16_t ushPTABPosNeg;
};


typedef struct
{
  sMDH mdh;
  void* previous;
  void* next;
  float* data;
} SiemensMdhNode;

typedef struct
{
  uint32_t matrix_size[3];
  uint32_t pat_matrix_size[3];
  uint32_t base_resolution;
  uint32_t phase_encoding_lines;
  uint32_t partitions;
  uint32_t dimensions;
  float    phase_resolution;
  float    slice_resolution;
  uint32_t dwell_time_us;
  uint32_t acceleration_factor_pe;
  uint32_t acceleration_factor_3d;
} SiemensBaseParameters;

class SiemensRawData
{
public:
  SiemensRawData();
  ~SiemensRawData();

  int ReadRawFile(char* filename);
  long GetNumberOfNodes();
  SiemensMdhNode* GetFirstNode();
  sMDH* GetMinValues();
  sMDH* GetMaxValues();

  SiemensBaseParameters GetBaseParameters() {return m_base_parameters;}
  
  int GetMeasYapsParameter(std::string parameter_name, std::string& value);

protected:
  int ReadMdhNode(std::ifstream* f);
  int DeleteNode(SiemensMdhNode* node);
  int DeleteNodeList();
  int UpdateMinMax();
  int ParseMeasYaps();

  SiemensMdhNode* m_first_node;
  SiemensMdhNode* m_last_node;
  long            m_number_of_nodes;
  sMDH            m_mdh_min;
  sMDH            m_mdh_max;
  bool            m_min_max_is_valid;

  std::map< std::string, std::string >  m_parameter_buffers;
  std::map< std::string, std::string >  m_meas_yaps;
  SiemensBaseParameters  m_base_parameters;
};
