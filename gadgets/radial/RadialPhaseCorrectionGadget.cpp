#include "RadialPhaseCorrectionGadget.h"
#include "hoNDArray_elemwise.h"
#include "hoArmadillo.h"
#include "hoNDArray_fileio.h"

#define _USE_MATH_DEFINES
#include <math.h>
#include <cmath>

#ifdef USE_OMP
#include <omp.h>
#endif

namespace Gadgetron{

  RadialPhaseCorrectionGadget::RadialPhaseCorrectionGadget()
    : slices_(-1)
    , sets_(-1)
    , channels_(-1)
    , profiles_counter_(0)
  {
  }

  int RadialPhaseCorrectionGadget::
  process_config(const mrd::Header& header)
  {
    if (header.encoding.size() != 1) {
      GDEBUG("This Gadget only supports one encoding space\n");
      return GADGET_FAIL;
    }

    // Get the encoding space and trajectory description
    mrd::EncodingSpaceType e_space = header.encoding[0].encoded_space;
    mrd::EncodingSpaceType r_space = header.encoding[0].recon_space;
    mrd::EncodingLimitsType e_limits = header.encoding[0].encoding_limits;

    slices_ = e_limits.slice ? e_limits.slice->maximum + 1 : 1;
    sets_ = e_limits.set ? e_limits.set->maximum + 1 : 1;

    if (header.acquisition_system_information) {
      channels_ = header.acquisition_system_information->receiver_channels.value_or(128);
    }

    mode_ = mode.value();
    order_ = order.value();
    profiles_ = profiles.value();

    if( profiles_ < 1 ) {
      GDEBUG("The number of profiles to estimate polynomial fit is too low.\n");
      return GADGET_FAIL;
    }

    fit_calculated = std::vector<bool>(sets_ * slices_, false);
    polyfit = std::vector<double>(channels_ * sets_ * slices_ * (order_ + 1), 0.0);
    profiles_queue = std::map<unsigned int, std::queue<std::unique_ptr<AcquisitionMessage>>>();

    return GADGET_OK;
  }

  int RadialPhaseCorrectionGadget
  ::process( GadgetContainerMessage<mrd::Acquisition> *m1)
  {

    // Pass any noise measurements down the chain
    //

    bool is_noise = m1->getObjectPtr()->head.flags.HasFlags(mrd::AcquisitionFlags::kIsNoiseMeasurement);
    if (is_noise) {
      if (this->next()->putq(m1) < 0) {
        GDEBUG("Failed to pass on noise samples.\n");
        return GADGET_FAIL;
      }
      return GADGET_OK;
    }

    // For now we require that this gadget is inserted before any coil reduction gadgets
    //

    if( channels_ != m1->getObjectPtr()->Coils() ){
      GDEBUG("Unexpected number of coils encountered. Did you insert the phase correction gadget after a coil reduction gadget? In that case invert the order of these gadgets\n");
      return GADGET_FAIL;
    }

    unsigned int slice = m1->getObjectPtr()->head.idx.slice.value_or(0);
    unsigned int set = m1->getObjectPtr()->head.idx.set.value_or(0);
    auto idx = set*slices_+slice;

    if( !fit_calculated[idx] ){

      // Enqueue the first 'profiles_' profiles...
      //

      profiles_queue[idx].push(std::unique_ptr<AcquisitionMessage>(m1));

      // ...before estimating the polynomial fit of order 'order_'
      //

      if( profiles_queue[idx].size() == profiles_ ){

        // Perform polynomial fit,
        // assemble system matrix A.
        //

        arma::mat A( profiles_, order_+1 );

        for( int m=0; m<profiles_; m++ ){

          double angle = get_projection_angle(m);

          for( int n=0; n<order_+1; n++ ){
            A(m,n) = pow( angle, double(n) );
          }
        }

        // Assemble right hand side
        //

        arma::mat b( profiles_, channels_ );
        //double prev_phase[channels_];
        std::vector<double> prev_phase(channels_);

        for( int m=0; m<profiles_; m++ ){

          AcquisitionMessage *mbq = profiles_queue[idx].front().release();
          profiles_queue[idx].pop();

          auto profile = mbq->getObjectPtr()->data;

          // A unique fit for each coil
          //

          for( unsigned int coil=0; coil<channels_; coil++ ){

            // 'arg' returns angles in the interval (-pi;pi)
            // Make sure that no discontinouities arise on the graph as they cannot be fitted
            //

            std::complex<float> sample = profile.get_data_ptr()[coil*profile.get_size(0)+(profile.get_size(0)>>1)];
            double phase = double(std::arg(sample));

            if( m>0 && std::abs(phase-prev_phase[coil])>M_PI ){

              // It appears as if phase wrapping has occurred, make correction...
              //

              if( phase<prev_phase[coil] )
                phase += 2.0*M_PI;
              else
                phase -= 2.0*M_PI;
            }

            b(m,coil) = phase;
            prev_phase[coil] = phase;
          }
        }

        // Linear least squares fit, i.e. solve "A^T A x = b"
        //

        std::vector<size_t> dims; dims.push_back(order_+1); dims.push_back(channels_);
        hoNDArray<double> vec( dims, &polyfit[set*(order_+1)*channels_*slices_+slice*(order_+1)*channels_] );

        arma::mat x = as_arma_matrix(vec);
        x = arma::solve(A.t()*A,A.t()*b);

        // Phase correct buffered profiles
        //

        for( int m=0; m<profiles_; m++ ){

          AcquisitionMessage *acq = profiles_queue[idx].front().release();
          profiles_queue[idx].pop();

          if(!acq) {
            GDEBUG("Unable to interpret data on message queue (3)\n");
            return GADGET_FAIL;
          }

          phase_correct(acq);

          if (this->next()->putq(acq) < 0) {
            GDEBUG("Failed to put data on queue\n");
            return GADGET_FAIL;
          }
        }
        fit_calculated[idx] = true;
      }
    }
    else{

      // Phase correct profile
      //

      phase_correct(m1);

      if (this->next()->putq(m1) < 0) {
        GDEBUG("Failed to put data on queue\n");
        return GADGET_FAIL;
      }
    }

    return GADGET_OK;
  }


  double RadialPhaseCorrectionGadget
  ::get_projection_angle( unsigned int idx )
  {
    if(!(mode_ == 2 || mode_ == 3 )){
      throw std::runtime_error("RadialPhaseCorrectionGadget: currently only trajectory modes 2 and 3 are supported (golden ratio)");;
    }

    double angle_step;
    if( mode_ == 2 )
      angle_step = M_PI/((std::sqrt(5.0)+1.0)*0.5); // GR_ORIGINAL
    else if( mode_ == 3 ){
      angle_step = M_PI*(3.0-std::sqrt(5.0))*0.5;   // GR_SMALLEST
    }
    return fmod(idx*angle_step, 2.0*M_PI);
  }

  void RadialPhaseCorrectionGadget
  ::phase_correct( GadgetContainerMessage<mrd::Acquisition> *m1 )
  {
    unsigned int slice = m1->getObjectPtr()->head.idx.slice.value_or(0);
    unsigned int set = m1->getObjectPtr()->head.idx.set.value_or(0);
    double angle = get_projection_angle(profiles_counter_);

    for( unsigned int coil=0; coil<channels_; coil++ ){

      double estimated_phase = 0.0;

      for( unsigned int i=0; i<order_+1; i++ ){

        double weight = polyfit[set*(order_+1)*channels_*slices_ +
                                slice*(order_+1)*channels_ +
                                coil*(order_+1) +
                                i ];

        double power = std::pow(angle, double(i));

        estimated_phase += (weight*power);
      }

      auto profile = m1->getObjectPtr()->data;
#ifdef USE_OMP
#pragma omp parallel for
#endif
      for( int i=0; i<profile.get_size(0); i++ ){
        std::complex<float> sample = profile.get_data_ptr()[coil*profile.get_size(0)+i];
        float phase = std::arg(sample);
        float mag = std::abs(sample);
        profile.get_data_ptr()[coil*profile.get_size(0)+i] = std::polar( mag, phase-float(estimated_phase) );
      }
    }
    profiles_counter_++;
  }

  GADGET_FACTORY_DECLARE(RadialPhaseCorrectionGadget)

} // namespace Gadgetron
