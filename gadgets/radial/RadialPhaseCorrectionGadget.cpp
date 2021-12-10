#include "RadialPhaseCorrectionGadget.h"
#include "hoNDArray_elemwise.h"
#include "hoArmadillo.h"
#include "hoNDArray_fileio.h"
#include "ismrmrd/xml.h"

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
  process_config( ACE_Message_Block *mb )
  {
    ISMRMRD::IsmrmrdHeader h;
    ISMRMRD::deserialize(mb->rd_ptr(),h);
    
    
    if (h.encoding.size() != 1) {
      GDEBUG("This Gadget only supports one encoding space\n");
      return GADGET_FAIL;
    }
    
    // Get the encoding space and trajectory description
    ISMRMRD::EncodingSpace e_space = h.encoding[0].encodedSpace;
    ISMRMRD::EncodingSpace r_space = h.encoding[0].reconSpace;
    ISMRMRD::EncodingLimits e_limits = h.encoding[0].encodingLimits;

    slices_ = e_limits.slice ? e_limits.slice->maximum + 1 : 1;
    sets_ = e_limits.set ? e_limits.set->maximum + 1 : 1;

    if (h.acquisitionSystemInformation) {
      channels_ = h.acquisitionSystemInformation->receiverChannels ? *h.acquisitionSystemInformation->receiverChannels : 128;
    }

    mode_ = mode.value();
    order_ = order.value();
    profiles_ = profiles.value();

    if( profiles_ < 1 ) {
      GDEBUG("The number of profiles to estimate polynomial fit is too low.\n");
      return GADGET_FAIL;
    }

    fit_calculated_ = boost::shared_array<bool>(new bool[sets_*slices_]);
    polyfit_ = boost::shared_array<double>(new double[(order_+1)*channels_*sets_*slices_]);   
    profiles_queue_ = boost::shared_array< ACE_Message_Queue<ACE_MT_SYNCH> >(new ACE_Message_Queue<ACE_MT_SYNCH>[slices_*sets_]);

    size_t bsize = sizeof(GadgetContainerMessage<ISMRMRD::AcquisitionHeader>)*profiles_*10;
    for( unsigned int i=0; i<slices_*sets_; i++ ){
      fit_calculated_[i] = false;
      profiles_queue_[i].high_water_mark(bsize);
      profiles_queue_[i].low_water_mark(bsize);
    }    
    
    return GADGET_OK;
  }
  
  int RadialPhaseCorrectionGadget
  ::process( GadgetContainerMessage<ISMRMRD::AcquisitionHeader> *m1,
             GadgetContainerMessage< hoNDArray< std::complex<float> > > *m2 )
  {

    // Pass any noise measurements down the chain
    //
    
    bool is_noise = m1->getObjectPtr()->isFlagSet(ISMRMRD::ISMRMRD_ACQ_IS_NOISE_MEASUREMENT);
    if (is_noise) { 
      if (this->next()->putq(m1) < 0) {
        GDEBUG("Failed to pass on noise samples.\n");
        return GADGET_FAIL;
      }
      return GADGET_OK;
    }

    // For now we require that this gadget is inserted before any coil reduction gadgets
    //

    if( channels_ != m1->getObjectPtr()->active_channels ){
      GDEBUG("Unexpected number of coils encountered. Did you insert the phase correction gadget after a coil reduction gadget? In that case invert the order of these gadgets\n");
      return GADGET_FAIL;
    }

    unsigned int slice = m1->getObjectPtr()->idx.slice;
    unsigned int set = m1->getObjectPtr()->idx.set;
    int idx = set*slices_+slice;

    if( !fit_calculated_[idx] ){

      // Enqueue the first 'profiles_' profiles...
      //
      
      profiles_queue_[idx].enqueue_tail(m1);

      // ...before estimating the polynomial fit of order 'order_'
      //

      if( profiles_queue_[idx].message_count() == profiles_ ){

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
        ACE_Message_Queue<ACE_MT_SYNCH>::ITERATOR iter(profiles_queue_[idx]);
        
        for( int m=0; m<profiles_; m++ ){                     
          
          ACE_Message_Block* mbq = 0x0;
          iter.next( mbq );
          iter.advance();
          
          if(!mbq) {
            GDEBUG("Unable to interpret data on message queue (1)\n");
            return GADGET_FAIL;
          }
          
          GadgetContainerMessage< hoNDArray< std::complex<float> > > *_profile = 
            AsContainerMessage< hoNDArray< std::complex<float> > >(mbq->cont());
        
          if(!_profile) {
            GDEBUG("Unable to interpret data on message queue (2)\n");
            return GADGET_FAIL;
          }
          
          hoNDArray< std::complex<float> > *profile = _profile->getObjectPtr();

          // A unique fit for each coil
          //

          for( unsigned int coil=0; coil<channels_; coil++ ){
            
            // 'arg' returns angles in the interval (-pi;pi)
            // Make sure that no discontinouities arise on the graph as they cannot be fitted
            //
            
            std::complex<float> sample = profile->get_data_ptr()[coil*profile->get_size(0)+(profile->get_size(0)>>1)];
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
        hoNDArray<double> vec( &dims, &polyfit_[set*(order_+1)*channels_*slices_+slice*(order_+1)*channels_] );

        arma::mat x = as_arma_matrix(vec);
        x = arma::solve(A.t()*A,A.t()*b);

        /*
        static int counter = 0;
        char filename[256];
        sprintf((char*)filename, "_polyfit_%d.real", counter);
        write_nd_array<double>( &vec, filename );
        */
        
        // Phase correct buffered profiles
        //

        for( int m=0; m<profiles_; m++ ){          

          ACE_Message_Block *mbq;
          if( profiles_queue_[idx].dequeue_head(mbq) < 0 ){
            GDEBUG("Message dequeue failed\n");
            GADGET_FAIL;
          }

          GadgetContainerMessage<ISMRMRD::AcquisitionHeader> *header = 
            AsContainerMessage<ISMRMRD::AcquisitionHeader>(mbq);
          
          if(!header) {
            GDEBUG("Unable to interpret data on message queue (3)\n");
            return GADGET_FAIL;
          }

          phase_correct(header);

          if (this->next()->putq(header) < 0) {
            GDEBUG("Failed to put data on queue\n");
            return GADGET_FAIL;
          }          
        }
        fit_calculated_[idx] = true;
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
  ::phase_correct( GadgetContainerMessage<ISMRMRD::AcquisitionHeader> *m1 )
  {
    unsigned int slice = m1->getObjectPtr()->idx.slice;
    unsigned int set = m1->getObjectPtr()->idx.set;
    double angle = get_projection_angle(profiles_counter_);

    for( unsigned int coil=0; coil<channels_; coil++ ){

      double estimated_phase = 0.0;

      for( unsigned int i=0; i<order_+1; i++ ){

        double weight = polyfit_[set*(order_+1)*channels_*slices_ +
                                 slice*(order_+1)*channels_ +
                                 coil*(order_+1) + 
                                 i ];

        double power = std::pow(angle, double(i));

        estimated_phase += (weight*power);
      }
      
      GadgetContainerMessage< hoNDArray< std::complex<float> > > *_profile = 
        AsContainerMessage<hoNDArray< std::complex<float> > >(m1->cont());
      
      if(!_profile) {
        GDEBUG("Unable to phase correct profile\n");
        return;
      }

      hoNDArray< std::complex<float> > *profile = _profile->getObjectPtr();      
#ifdef USE_OMP
#pragma omp parallel for
#endif
      for( int i=0; i<profile->get_size(0); i++ ){
        std::complex<float> sample = profile->get_data_ptr()[coil*profile->get_size(0)+i];
        float phase = std::arg(sample);
        float mag = std::abs(sample);
        profile->get_data_ptr()[coil*profile->get_size(0)+i] = std::polar( mag, phase-float(estimated_phase) );
      }
    }
    profiles_counter_++;
  }

  GADGET_FACTORY_DECLARE(RadialPhaseCorrectionGadget)

} // namespace Gadgetron
