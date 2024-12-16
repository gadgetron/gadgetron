/** \file   EPIReconXObject.h
    \brief  Define the symbols and implement functionality for EPI X reconstruction operator
    \author Souheil Inati
*/

#pragma once

#include "hoNDArray.h"

namespace Gadgetron { namespace EPI {

enum EPIType
{
  FLAT,
  TRAPEZOID,
  SINUSOID,
  ARBITRARY
};

enum EPIReceiverPhaseType
{
  NONE,
  EVEN,
  FULL
};

template <typename T> class EPIReconXObject
{
 public:
  EPIReconXObject();
  virtual ~EPIReconXObject();

  hoNDArray <float> getTrajectoryPos();
  hoNDArray <float> getTrajectoryNeg();

  hoNDArray <float> filterPos_;
  hoNDArray <float> filterNeg_;
  float slicePosition[3];

  virtual int computeTrajectory()=0;

  virtual int apply(mrd::AcquisitionHeader &hdr_in,  hoNDArray <T> &data_in,
		    mrd::AcquisitionHeader &hdr_out, hoNDArray <T> &data_out)=0;
  EPIReceiverPhaseType rcvType_;

 protected:
  hoNDArray <float> trajectoryPos_;
  hoNDArray <float> trajectoryNeg_;

};

template <typename T> EPIReconXObject<T>::EPIReconXObject()
{
}

template <typename T> EPIReconXObject<T>::~EPIReconXObject()
{
}

template <typename T> hoNDArray<float> EPIReconXObject<T>::getTrajectoryPos()
{
  return trajectoryPos_;
}

template <typename T> hoNDArray<float> EPIReconXObject<T>::getTrajectoryNeg()
{
  return trajectoryNeg_;
}

}}
