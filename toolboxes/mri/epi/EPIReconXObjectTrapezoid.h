/** \file   EPIReconXObjectTrapezoid.h
    \brief  Implement functionality for EPI X reconstruction operator for Trapezoidal type
    \author Souheil Inati
*/

#pragma once

#include "EPIReconXObject.h"
#include "hoArmadillo.h"
#include "hoNDArray_elemwise.h"
#include "hoNDArray_linalg.h"
#include "gadgetronmath.h"
#include <complex>

namespace Gadgetron { namespace EPI {

template <typename T> class EPIReconXObjectTrapezoid : public EPIReconXObject<T>
{
 public:
  EPIReconXObjectTrapezoid();
  virtual ~EPIReconXObjectTrapezoid();

  virtual int computeTrajectory();

  virtual int apply(mrd::AcquisitionHeader &hdr_in, hoNDArray <T> &data_in,
		    mrd::AcquisitionHeader &hdr_out, hoNDArray <T> &data_out);

  using EPIReconXObject<T>::filterPos_;
  using EPIReconXObject<T>::filterNeg_;
  using EPIReconXObject<T>::slicePosition;
  using EPIReconXObject<T>::rcvType_;

  bool  balanced_;
  float rampUpTime_;
  float rampDownTime_;
  float flatTopTime_;
  float acqDelayTime_;
  int   numSamples_;
  float dwellTime_;
  int   encodeNx_;
  float encodeFOV_;
  int   reconNx_;
  float reconFOV_;

 protected:
  using EPIReconXObject<T>::trajectoryPos_;
  using EPIReconXObject<T>::trajectoryNeg_;

  hoNDArray <T> Mpos_;
  hoNDArray <T> Mneg_;
  bool operatorComputed_;

  float calcOffCenterDistance(mrd::AcquisitionHeader& hdr_in);
};

template <typename T> EPIReconXObjectTrapezoid<T>::EPIReconXObjectTrapezoid()
{
  rcvType_ = EVEN;
  balanced_ = true;
  rampUpTime_ = 0.0;
  rampDownTime_ = 0.0;
  flatTopTime_ = 0.0;
  acqDelayTime_ = 0.0;
  numSamples_ = 0.0;
  dwellTime_ = 0.0;
  encodeNx_ = 0;
  reconNx_ = 0;
  encodeFOV_ = 0.0;
  reconFOV_ = 0.0;
  operatorComputed_ = false;
}

template <typename T> EPIReconXObjectTrapezoid<T>::~EPIReconXObjectTrapezoid()
{
}

template <typename T> int EPIReconXObjectTrapezoid<T>::computeTrajectory()
{

  // Initialize the k-space trajectory arrays
  trajectoryPos_.create(numSamples_);
  Gadgetron::clear(trajectoryPos_);
  trajectoryNeg_.create(numSamples_);
  Gadgetron::clear(trajectoryNeg_);

  // Temporary trajectory for a symmetric readout
  // first calculate the integral with G = 1;
  int nK = numSamples_;
  hoNDArray <float> k(nK);
  float t;
  int n;

  //GDEBUG_STREAM("Dwell = " << dwellTime_ << "    acqDelayTime = " << acqDelayTime_ << std::endl);
  //GDEBUG_STREAM("rampUpTime = " << rampUpTime_ << "    flatTopTime = " << flatTopTime_ << "    rampDownTime = " << rampDownTime_ << std::endl);

  // Some timings
  float totTime = rampUpTime_ + flatTopTime_ + rampDownTime_;
  float readTime = dwellTime_ * numSamples_;

  // Fix the acqDelayTime for balanced acquisitions
  if (balanced_) {
    acqDelayTime_ = 0.5 * (totTime - readTime);
  }

  // Some Areas
  float totArea = 0.5*rampUpTime_ + flatTopTime_ + 0.5*rampDownTime_;
  float readArea =  0.5*rampUpTime_ + flatTopTime_ + 0.5*rampDownTime_;
  if (rampUpTime_ > 0.0) {
      readArea -= 0.5*(acqDelayTime_)*(acqDelayTime_)/rampUpTime_;
  }
  if (rampDownTime_ > 0.0) {
      readArea -= 0.5*(totTime - (acqDelayTime_+readTime))*(totTime - (acqDelayTime_+readTime))/rampDownTime_;
  }

  // Prephase is set so that k=0 is halfway through the readout time
  float prePhaseArea = 0.5 * totArea;

  // The scale is set so that the read out area corresponds to the number of encoded points
  float scale = encodeNx_ /readArea;

  for (n=0; n<nK; n++)
  {
    t = (n+1.0)*dwellTime_ + acqDelayTime_;  // end of the dwell time
    if (t <= rampUpTime_) {
      // on the ramp up
      k[n] = 0.5 / rampUpTime_ * t*t;
    }
    else if ((t > rampUpTime_) && (t <= (rampUpTime_+flatTopTime_))) {
      // on the flat top
      k[n] = 0.5*rampUpTime_ + (t - rampUpTime_);
    }
    else {
      // on the ramp down
      float v = (rampUpTime_+flatTopTime_+rampDownTime_-t);
      k[n] = 0.5*rampUpTime_ + flatTopTime_ + 0.5*rampDownTime_ - 0.5/rampDownTime_*v*v;
    }
    //GDEBUG_STREAM(n << ":  " << t << "  " << k[n] << " " << std::endl);
  }

  // Fill the positive and negative trajectories
  for (n=0; n<numSamples_; n++)
  {
    trajectoryPos_[n] = scale * (k[n] - prePhaseArea);
    trajectoryNeg_[n] = scale * (-1.0*k[n] + totArea - prePhaseArea);
    //GDEBUG_STREAM(n << ":  " << trajectoryPos_[n] << "  " << trajectoryNeg_[n] << std::endl);
  }

  // reset the operatorComputed_ flag
  operatorComputed_ = false;

  return(0);
}


template <typename T> int EPIReconXObjectTrapezoid<T>::apply(mrd::AcquisitionHeader &hdr_in, hoNDArray <T> &data_in,
		    mrd::AcquisitionHeader &hdr_out, hoNDArray <T> &data_out)
{
  if (!operatorComputed_) {
    // Compute the reconstruction operator
    int Km = std::floor(encodeNx_ / 2.0);
    int Ne = 2*Km + 1;
    int p,q; // counters

    // resize the reconstruction operator
    Mpos_.create(reconNx_,numSamples_);
    Mneg_.create(reconNx_,numSamples_);

    // evenly spaced k-space locations
    arma::vec keven = arma::linspace<arma::vec>(-Km, Km, Ne);
    //keven.print("keven =");

    // image domain locations [-0.5,...,0.5)
    arma::vec x = arma::linspace<arma::vec>(-0.5,(reconNx_-1.)/(2.*reconNx_),reconNx_);
    //x.print("x =");

    // DFT operator
    // Going from k space to image space, we use the IFFT sign convention
    arma::cx_mat F(reconNx_, Ne);
    double fftscale = 1.0 / std::sqrt((double)Ne);
    for (p=0; p<reconNx_; p++) {
      for (q=0; q<Ne; q++) {
	F(p,q) = fftscale * std::exp(std::complex<double>(0.0,1.0*2*M_PI*keven(q)*x(p)));
      }
    }
    //F.print("F =");

    // forward operators
    arma::mat Qp(numSamples_, Ne);
    arma::mat Qn(numSamples_, Ne);
    for (p=0; p<numSamples_; p++) {
      //GDEBUG_STREAM(trajectoryPos_(p) << "    " << trajectoryNeg_(p) << std::endl);
      for (q=0; q<Ne; q++) {
	Qp(p,q) = sinc(trajectoryPos_(p)-keven(q));
	Qn(p,q) = sinc(trajectoryNeg_(p)-keven(q));
      }
    }

    //Qp.print("Qp =");
    //Qn.print("Qn =");

    // recon operators
    arma::cx_mat Mp(reconNx_,numSamples_);
    arma::cx_mat Mn(reconNx_,numSamples_);
    Mp = F * arma::pinv(Qp);
    Mn = F * arma::pinv(Qn);

    /////    Compute the off-center correction:     /////

    // Compute the off-center distance in the RO direction:
    float roOffCenterDistance = calcOffCenterDistance( hdr_in );

    arma::Col<typename realType<T>::Type> my_keven = arma::linspace< arma::Col<typename realType<T>::Type> >(0, numSamples_ -1, numSamples_);
    // find the offset:
    // PV: maybe find not just exactly 0, but a very small number?
    arma::Col<typename realType<T>::Type> trajectoryPosArma = as_arma_col(trajectoryPos_);
    arma::uvec n = find( trajectoryPosArma==0, 1, "first");
    my_keven -= arma::as_scalar(n);
    // Scale it:
    // We have to find the maximum k-trajectory (absolute) increment:
    arma::Col<typename realType<T>::Type> Delta_k = arma::abs( trajectoryPosArma.subvec(1,numSamples_-1) - trajectoryPosArma.subvec(0,numSamples_-2) );
    my_keven *= Delta_k.max();

    // off-center corrections:
    arma::Col<T> myExponent = arma::zeros< arma::Col<T> >(numSamples_);
    myExponent.set_imag( 2*M_PI*roOffCenterDistance/encodeFOV_*(trajectoryPosArma-my_keven) );
    arma::Col<T> offCenterCorrN = arma::exp( myExponent );
    myExponent.set_imag( 2*M_PI*roOffCenterDistance/encodeFOV_*(as_arma_col(trajectoryNeg_)+my_keven) );
    arma::Col<T> offCenterCorrP = arma::exp( myExponent );

    //    GDEBUG_STREAM("roOffCenterDistance_: " << roOffCenterDistance_ << ";       encodeFOV_: " << encodeFOV_);
    //    for (q=0; q<numSamples_; q++) {
    //      GDEBUG_STREAM("keven(" << q << "): " << my_keven(q) << ";       trajectoryPosArma(" << q << "): " << trajectoryPosArma(q) );
    //      GDEBUG_STREAM("offCenterCorrP(" << q << "):" << offCenterCorrP(q) );
    //    }

    // Finally, combine the off-center correction with the recon operator:
    Mp = Mp * diagmat(offCenterCorrP);
    Mn = Mn * diagmat(offCenterCorrN);
    // and save it into the NDArray members:
    for (p=0; p<reconNx_; p++) {
      for (q=0; q<numSamples_; q++) {
        Mpos_(p,q) = Mp(p,q);
        Mneg_(p,q) = Mn(p,q);
      }
    }

    //Mp.print("Mp =");
    //Mn.print("Mn =");

    // set the operator computed flag
    operatorComputed_ = true;
  }

  // convert to armadillo representation of matrices and vectors
  //arma::Mat<typename stdType<T>::Type> adata_in = as_arma_matrix(&data_in);
  //arma::Mat<typename stdType<T>::Type> adata_out = as_arma_matrix(&data_out);

  // Apply it
  if (hdr_in.flags.HasFlags(mrd::AcquisitionFlags::kIsReverse)) {

    // Negative readout
    // adata_out = as_arma_matrix(&Mneg_) * adata_in;
      Gadgetron::gemm(data_out, Mneg_, data_in);
  } else {
    // Forward readout
    // adata_out = as_arma_matrix(&Mpos_) * adata_in;
      Gadgetron::gemm(data_out, Mpos_, data_in);
  }

  // Copy the input header to the output header and set the size and the center sample
  hdr_out = hdr_in;
  hdr_out.center_sample = reconNx_/2;

  return 0;
}

template <typename T> float EPIReconXObjectTrapezoid<T>::calcOffCenterDistance(mrd::AcquisitionHeader& hdr_in)
{
  // armadillo vectors with the position and readout direction:
  arma::fvec pos(3);
  arma::fvec RO_dir(3);

  for (int i=0; i<3; i++) {
    pos(i)    = hdr_in.position[i];
    RO_dir(i) = hdr_in.read_dir[i];
  }

  float roOffCenterDistance = dot(pos, RO_dir);
  GDEBUG_STREAM("roOffCenterDistance: " << roOffCenterDistance );

  return roOffCenterDistance;

}

}}
