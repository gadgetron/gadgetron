#pragma once

#define DECLARE_MATRIX_OPERATOR_DEVICE_SUPPORT(COMPONENT)                                  \
                                                                                           \
  virtual void set_device( int device )                                                    \
  {                                                                                        \
    if( device<0 ){                                                                        \
      if( cudaGetDevice( &device_ ) != cudaSuccess ){                                      \
	std::cerr << "##COMPONENT:: unable to get current device." << std::endl;           \
	device_ = 0;                                                                       \
      }                                                                                    \
    }                                                                                      \
    else                                                                                   \
      device_ = device;                                                                    \
  }                                                                                        \
                                                                                           \
protected:                                                                                 \
  virtual int set_device()                                                                 \
  {                                                                                        \
    if( cudaGetDevice( &old_device_ ) != cudaSuccess ){                                    \
      std::cerr << "##COMPONENT:: unable to get current device." << std::endl ;            \
      return -1;                                                                           \
    }                                                                                      \
    if( device_ != old_device_ && cudaSetDevice(device_) != cudaSuccess) {                 \
      std::cerr << "##COMPONENT:: unable to set device " << device_ << std::endl;          \
      return -1;                                                                           \
    }                                                                                      \
    return 0;                                                                              \
  }                                                                                        \
                                                                                           \
  virtual int restore_device()                                                             \
  {                                                                                        \
    if( device_ != old_device_ && cudaSetDevice(old_device_) != cudaSuccess) {             \
      std::cerr << "##COMPONENT:: unable to restore device " << old_device_ << std::endl;  \
      return -1;                                                                           \
    }                                                                                      \
    return 0;                                                                              \
  }                                                                                        \
                                                                                           \
protected:                                                                                 \
  int device_;                                                                             \
                                                                                           \
private:                                                                                   \
  int old_device_; 
