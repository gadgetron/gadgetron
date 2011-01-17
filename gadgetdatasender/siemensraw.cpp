
#include "siemensraw.hpp"

#include <iostream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>

#define SET_MIN(A,B) A = B < A ? B : A
#define SET_MAX(A,B) A = B > A ? B : A

SiemensRawData::SiemensRawData()
  :m_first_node         (NULL),
   m_last_node          (NULL),
   m_number_of_nodes    (0),
   m_min_max_is_valid   (false)
{

}

SiemensRawData::~SiemensRawData()
{
  DeleteNodeList();
}

int SiemensRawData::ReadRawFile(char* filename)
{
  unsigned int header_length = 0;
  unsigned int file_length = 0;
  std::ifstream f(filename, std::ios::in | std::ios::binary);
  
  //Reading the oveall length of the header(s)
  f.read(reinterpret_cast<char*>(&header_length),sizeof(unsigned int));

  //Figure out the number of parameter buffers
  uint32_t number_of_parameter_buffers = 0;
  f.read(reinterpret_cast<char*>(&number_of_parameter_buffers),sizeof(uint32_t));

  std::cout << "Number of parameter buffers: " << number_of_parameter_buffers << std::endl;
  
  for (unsigned int i = 0; i < number_of_parameter_buffers; i++) {
    char buffer_name[1024];
    char* tmp_buffer = 0;
    uint32_t buffer_length;

    f.getline(buffer_name,1024,'\0');
    f.read(reinterpret_cast<char*>(&buffer_length),sizeof(uint32_t));
   
    try {
      tmp_buffer = new char[buffer_length];
    } catch (...) {
      std::cout << "Unable to allocate temporary memory for buffer" << std::endl;
      return -1;
    }
    
    f.read(tmp_buffer,buffer_length);
    m_parameter_buffers[std::string(buffer_name)] = std::string(tmp_buffer,buffer_length);

    if (tmp_buffer) delete [] tmp_buffer;

  }

  f.seekg(0,std::ios_base::end);
  file_length = f.tellg();

  f.seekg(header_length, std::ios_base::beg);
  
  while (ReadMdhNode(&f) > 0 && ((file_length-f.tellg())>(unsigned int)sizeof(sMDH)))
  {
    //...
  }

  f.close();

  UpdateMinMax();

  //std::cout << "Done reading. Total mdh count: " << m_number_of_nodes << std::endl;
  ParseMeasYaps();
  return m_number_of_nodes;
}

int SiemensRawData::ReadMdhNode(std::ifstream* f)
{
  SiemensMdhNode* new_node = NULL;

  if (f->eof())
  {
    return -1;
  }

  try {
    new_node = new SiemensMdhNode;
  } catch (...) {
    std::cout << "SiemensRawData::ReadMdhNode: Failed to allocate new node" << std::endl;
    return -1;
  }

  new_node->data = NULL;
  new_node->previous = NULL;
  new_node->next = NULL;
  
  try
  {
    f->read(reinterpret_cast<char*>(&(new_node->mdh)),sizeof(sMDH));
  } catch (...) { 
    std::cout << "SiemensRawData::ReadMdhNode: Unable to read mdh" << std::endl;
    delete new_node;
    return -1;
  }
  
  //Is this end of acquisition
  if (new_node->mdh.aulEvalInfoMask[0] & 0x0001)
  {
    delete new_node;
    return -1;
  }

  try
  {
    new_node->data = new float[new_node->mdh.ushSamplesInScan * 2];
  } catch (...) {
    std::cout << "SiemensRawData::ReadMdhNode: Unable to allocate memory in data node" << std::endl;
    delete new_node;
    return -1;
  }
 
  try
  {
    f->read(reinterpret_cast<char*>(new_node->data),sizeof(float)*new_node->mdh.ushSamplesInScan*2);
  } catch (...) { 
    std::cout << "SiemensRawData::ReadMdhNode: Unable to read data for data node" << std::endl;
    delete [] new_node->data;
    delete new_node;
    return -1;
  }  
  
  if (m_first_node == NULL) m_first_node = new_node;
  if (m_last_node == NULL)
  {
    m_last_node = new_node;
    m_last_node->previous = NULL;
    m_last_node->next = NULL;
  }
  else
  {
    m_last_node->next = new_node;
    new_node->previous = m_last_node;
    m_last_node = new_node;
  }
  m_number_of_nodes++;
  
  m_min_max_is_valid = false;

  return 1;
}

int SiemensRawData::DeleteNode(SiemensMdhNode* node)
{
  if (node)
  {
    SiemensMdhNode* prev = (SiemensMdhNode*)node->previous;
    SiemensMdhNode* next = (SiemensMdhNode*)node->next;
  
    if (node->data) delete [] node->data;
    delete node;

    if (prev) prev->next = next;
    if (next) next->previous = prev;

  }

  m_min_max_is_valid = false;

  return 1;
}

int SiemensRawData::DeleteNodeList()
{
  SiemensMdhNode* next;
  while(m_first_node)
  {
    next = (SiemensMdhNode*)m_first_node->next;
    DeleteNode(m_first_node);
    m_first_node = next;
  }

  return 1;
}

int SiemensRawData::UpdateMinMax()
{
  SiemensMdhNode* next = m_first_node;

  m_mdh_min = next->mdh;
  m_mdh_max = next->mdh;

  while (next)
  {
    SET_MIN(m_mdh_min.ulScanCounter,         next->mdh.ulScanCounter);
    SET_MAX(m_mdh_max.ulScanCounter,         next->mdh.ulScanCounter);
    SET_MIN(m_mdh_min.ulTimeStamp,           next->mdh.ulTimeStamp);
    SET_MAX(m_mdh_max.ulTimeStamp,           next->mdh.ulTimeStamp);
    SET_MIN(m_mdh_min.ulPMUTimeStamp,        next->mdh.ulPMUTimeStamp);
    SET_MAX(m_mdh_max.ulPMUTimeStamp,        next->mdh.ulPMUTimeStamp);
    SET_MIN(m_mdh_min.ushSamplesInScan,      next->mdh.ushSamplesInScan);
    SET_MAX(m_mdh_max.ushSamplesInScan,      next->mdh.ushSamplesInScan);
    SET_MIN(m_mdh_min.ushUsedChannels,       next->mdh.ushUsedChannels);
    SET_MAX(m_mdh_max.ushUsedChannels,       next->mdh.ushUsedChannels);
    SET_MIN(m_mdh_min.sLC.ushLine,           next->mdh.sLC.ushLine);
    SET_MAX(m_mdh_max.sLC.ushLine,           next->mdh.sLC.ushLine);
    SET_MIN(m_mdh_min.sLC.ushAcquisition,    next->mdh.sLC.ushAcquisition);
    SET_MAX(m_mdh_max.sLC.ushAcquisition,    next->mdh.sLC.ushAcquisition);
    SET_MIN(m_mdh_min.sLC.ushSlice,          next->mdh.sLC.ushSlice);
    SET_MAX(m_mdh_max.sLC.ushSlice,          next->mdh.sLC.ushSlice);
    SET_MIN(m_mdh_min.sLC.ushPartition,      next->mdh.sLC.ushPartition);
    SET_MAX(m_mdh_max.sLC.ushPartition,      next->mdh.sLC.ushPartition);
    SET_MIN(m_mdh_min.sLC.ushEcho,           next->mdh.sLC.ushEcho);
    SET_MAX(m_mdh_max.sLC.ushEcho,           next->mdh.sLC.ushEcho);
    SET_MIN(m_mdh_min.sLC.ushPhase,          next->mdh.sLC.ushPhase);
    SET_MAX(m_mdh_max.sLC.ushPhase,          next->mdh.sLC.ushPhase);
    SET_MIN(m_mdh_min.sLC.ushRepetition,     next->mdh.sLC.ushRepetition);
    SET_MAX(m_mdh_max.sLC.ushRepetition,     next->mdh.sLC.ushRepetition);
    SET_MIN(m_mdh_min.sLC.ushSet,            next->mdh.sLC.ushSet);
    SET_MAX(m_mdh_max.sLC.ushSet,            next->mdh.sLC.ushSet);
    SET_MIN(m_mdh_min.sLC.ushSeg,            next->mdh.sLC.ushSeg);
    SET_MAX(m_mdh_max.sLC.ushSeg,            next->mdh.sLC.ushSeg);
    SET_MIN(m_mdh_min.sLC.ushIda,            next->mdh.sLC.ushIda);
    SET_MAX(m_mdh_max.sLC.ushIda,            next->mdh.sLC.ushIda);
    SET_MIN(m_mdh_min.sLC.ushIdb,            next->mdh.sLC.ushIdb);
    SET_MAX(m_mdh_max.sLC.ushIdb,            next->mdh.sLC.ushIdb);
    SET_MIN(m_mdh_min.sLC.ushIdc,            next->mdh.sLC.ushIdc);
    SET_MAX(m_mdh_max.sLC.ushIdc,            next->mdh.sLC.ushIdc);
    SET_MIN(m_mdh_min.sLC.ushIdd,            next->mdh.sLC.ushIdd);
    SET_MAX(m_mdh_max.sLC.ushIdd,            next->mdh.sLC.ushIdd);
    SET_MIN(m_mdh_min.sLC.ushIde,            next->mdh.sLC.ushIde);
    SET_MAX(m_mdh_max.sLC.ushIde,            next->mdh.sLC.ushIde);
    SET_MIN(m_mdh_min.ushKSpaceCentreColumn, next->mdh.ushKSpaceCentreColumn);
    SET_MAX(m_mdh_max.ushKSpaceCentreColumn, next->mdh.ushKSpaceCentreColumn);
    SET_MIN(m_mdh_min.ushCoilSelect,         next->mdh.ushCoilSelect);
    SET_MAX(m_mdh_max.ushCoilSelect,         next->mdh.ushCoilSelect);
    SET_MIN(m_mdh_min.fReadOutOffcentre,     next->mdh.fReadOutOffcentre);
    SET_MAX(m_mdh_max.fReadOutOffcentre,     next->mdh.fReadOutOffcentre);
    SET_MIN(m_mdh_min.ulTimeSinceLastRF,     next->mdh.ulTimeSinceLastRF);
    SET_MAX(m_mdh_max.ulTimeSinceLastRF,     next->mdh.ulTimeSinceLastRF);
    SET_MIN(m_mdh_min.ushKSpaceCentreLineNo, next->mdh.ushKSpaceCentreLineNo);
    SET_MAX(m_mdh_max.ushKSpaceCentreLineNo, next->mdh.ushKSpaceCentreLineNo);
    SET_MIN(m_mdh_min.ushKSpaceCentrePartitionNo, next->mdh.ushKSpaceCentrePartitionNo);
    SET_MAX(m_mdh_max.ushKSpaceCentrePartitionNo, next->mdh.ushKSpaceCentrePartitionNo);
    for (int i = 0; i < MDH_NUMBEROFICEPROGRAMPARA; i++)
    {
      SET_MIN(m_mdh_min.aushIceProgramPara[i], next->mdh.aushIceProgramPara[i]);
      SET_MAX(m_mdh_max.aushIceProgramPara[i], next->mdh.aushIceProgramPara[i]);
    }
    for (int i = 0; i < MDH_FREEHDRPARA; i++)
    {
      SET_MIN(m_mdh_min.aushFreePara[i], next->mdh.aushFreePara[i]);
      SET_MAX(m_mdh_max.aushFreePara[i], next->mdh.aushFreePara[i]);
    }
    SET_MIN(m_mdh_min.ushChannelId,          next->mdh.ushChannelId);
    SET_MAX(m_mdh_max.ushChannelId,          next->mdh.ushChannelId);

    next = (SiemensMdhNode*)next->next;
  }

  m_min_max_is_valid = true;

  return 1;
}


long SiemensRawData::GetNumberOfNodes()
{
  return m_number_of_nodes;
}

SiemensMdhNode* SiemensRawData::GetFirstNode()
{
  return m_first_node;
}

sMDH* SiemensRawData::GetMinValues()
{
  if (!m_min_max_is_valid) UpdateMinMax();
  return &m_mdh_min;
}


sMDH* SiemensRawData::GetMaxValues()
{
  if (!m_min_max_is_valid) UpdateMinMax();
  return &m_mdh_max;
}

int SiemensRawData::GetMeasYapsParameter(std::string parameter_name, std::string& value)
{
  std::map<std::string, std::string>::iterator it;
  if ((it = m_meas_yaps.find(parameter_name)) == m_meas_yaps.end()) {
    return -1;
  } 

  value = it->second;
  return 0;
}

int SiemensRawData::ParseMeasYaps()
{
  static const char fname[] = "SiemensRawData::ParseMeasYaps";

  if (m_parameter_buffers.find("MeasYaps") == m_parameter_buffers.end()) {
    std::cout << fname << ": Unable to find MeasYaps buffer" << std::endl;
    return -1;
  }

  std::stringstream s(m_parameter_buffers["MeasYaps"],std::ios_base::in | std::ios_base::out);
  std::string line_str;
  
  while (!std::getline(s,line_str).eof()) {
    std::string param;
    std::string value;
    size_t p = line_str.find_first_of(" ");
    param = line_str.substr(0,p);
    p = line_str.find_first_of("=");
    value = line_str.substr(p+2);
    m_meas_yaps[param] = value;
  }

  //Get some parameters
  m_base_parameters.phase_resolution = atof(m_meas_yaps["sKSpace.dPhaseResolution"].c_str());
  m_base_parameters.slice_resolution = atof(m_meas_yaps["sKSpace.dSliceResolution"].c_str());
  m_base_parameters.base_resolution = atoi(m_meas_yaps["sKSpace.lBaseResolution"].c_str());
  m_base_parameters.phase_encoding_lines = atoi(m_meas_yaps["sKSpace.lPhaseEncodingLines"].c_str());
  m_base_parameters.partitions = atoi(m_meas_yaps["sKSpace.lPartitions"].c_str());
  sscanf(m_meas_yaps["sKSpace.ucDimension"].c_str(),"0x%X",&m_base_parameters.dimensions);
  m_base_parameters.acceleration_factor_pe = atoi(m_meas_yaps["sPat.lAccelFactPE"].c_str());
  m_base_parameters.acceleration_factor_3d = atoi(m_meas_yaps["sPat.lAccelFact3D"].c_str());

  m_base_parameters.matrix_size[0] = m_base_parameters.base_resolution;
  m_base_parameters.matrix_size[1] = static_cast<uint32_t>(m_base_parameters.phase_encoding_lines/ m_base_parameters.phase_resolution);

  if (m_base_parameters.dimensions > 2) {
    m_base_parameters.matrix_size[2] = static_cast<uint32_t>(m_base_parameters.partitions/ m_base_parameters.slice_resolution);
  } else {
    m_base_parameters.matrix_size[2] = 1;
  }

  m_base_parameters.pat_matrix_size[0] = m_base_parameters.matrix_size[0];
  m_base_parameters.pat_matrix_size[1] = ((m_base_parameters.matrix_size[1]-1)/m_base_parameters.acceleration_factor_pe + 1)*m_base_parameters.acceleration_factor_pe;
  m_base_parameters.pat_matrix_size[2] = ((m_base_parameters.matrix_size[2]-1)/m_base_parameters.acceleration_factor_3d + 1)*m_base_parameters.acceleration_factor_3d;
  
  return 0;
}
