#ifndef GRAPPACALIBRATIONBUFFER_H
#define GRAPPACALIBRATIONBUFFER_H

#include <vector>
#include <string.h>
#include <memory>
#include <complex>

#include "gadgetrongrappa_export.h"
#include "GadgetMRIHeaders.h"
#include "hoNDArray.h"
#include "GrappaWeights.h"
#include "GrappaWeightsCalculator.h"

class EXPORTGADGETSGRAPPA CalibrationBufferCounter
{

 public:
  CalibrationBufferCounter(unsigned int lines)  {
    lines_sampled_ = std::vector<unsigned int>(lines,0);
    memset(position_, 0, 3*sizeof(float));
    memset(quarternion_, 0, 4*sizeof(float));
  }


  virtual ~CalibrationBufferCounter() {}

  int update_line(unsigned int ky_index, float* position, float* quarternion)
  {
    int ret_val = 0;

    if (!quarterion_equal(quarternion) || !position_equal(position)) {
      for (unsigned int i = 0; i < lines_sampled_.size(); i++) {
	lines_sampled_[i] = 0;
      }
      memcpy(position_,position,3*sizeof(float));
      memcpy(quarternion_,quarternion,4*sizeof(float));
      ret_val = 1;
    }

    if (ky_index >= lines_sampled_.size()) {
      return -1;
    }

    lines_sampled_[ky_index] = 1;

    return ret_val;
  }

  int get_region_of_support(unsigned int& min_ky_index, unsigned int& max_ky_index) {
    
    unsigned int current_start_line = 0;
    min_ky_index = 0;
    max_ky_index = 0;
    while (current_start_line < lines_sampled_.size() ) {
      while ((current_start_line < lines_sampled_.size()) && (lines_sampled_[current_start_line] == 0) ) {	
       	current_start_line++;
      }
      if (current_start_line >= lines_sampled_.size()) continue;

      unsigned int region_start = current_start_line;
      while ((current_start_line < lines_sampled_.size()) && (lines_sampled_[current_start_line] > 0)) {	
       	current_start_line++;
      }
      unsigned int region_end = current_start_line-1;
      if ((region_start < region_end) && ((region_end-region_start) > (max_ky_index-min_ky_index))) {
	min_ky_index = region_start;
	max_ky_index = region_end;
      }
    }
    return 0;
  }

 protected:
  float           position_[3];
  float           quarternion_[4];   

  bool position_equal(float* position) {
    for (unsigned int i = 0; i < 3; i++) {
      if (position_[i] != position[i]) return false;
    }
    return true;
  }

  bool quarterion_equal(float* quarternion) {
    for (unsigned int i = 0; i < 4; i++) {
      if (quarternion_[i] != quarternion[i]) return false;
    }
    return true;
  }

 private:
  std::vector<unsigned int> lines_sampled_;

};

class EXPORTGADGETSGRAPPA GrappaCalibrationBuffer
{

 public:
  GrappaCalibrationBuffer(std::vector<unsigned int> dimensions, 
			  boost::shared_ptr< GrappaWeights<float> > w,
			  GrappaWeightsCalculator<float>* weights_calculator);
  virtual ~GrappaCalibrationBuffer() {}

  int add_data(GadgetMessageAcquisition* m1, hoNDArray< std::complex<float> >* m2);

 private:
  hoNDArray< std::complex<float> > buffer_;
  std::vector<unsigned int> dimensions_;
  boost::shared_ptr< GrappaWeights<float> > weights_;
  GrappaWeightsCalculator<float>* weights_calculator_;
  CalibrationBufferCounter buffer_counter_;

  unsigned int biggest_gap_current_;
  unsigned int acceleration_factor_;
  unsigned int last_line_;
  bool weights_invalid_;
};

#endif
