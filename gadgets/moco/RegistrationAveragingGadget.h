#ifndef RegistrationAveragingGadget_H
#define RegistrationAveragingGadget_H

#include "Gadget.h"
#include "GadgetMRIHeaders.h"
#include "hoNDArray.h"
#include "complext.h"
#include "PhysioInterpolationGadget.h"
#include "GadgetronTimer.h"
#include "gadgetron_moco_export.h"
#include "hoNDArray_fileio.h"

#ifdef USE_CUDA
#include "cuNDArray_reductions.h"
#endif // USE_CUDA

#include <ismrmrd/ismrmrd.h>
#include <complex>

#include <queue>

namespace Gadgetron{  

  template<class ARRAY_TYPE, unsigned int D> class opticalFlowSolver;
  
  /**
     This is an abstract gadget class and consequently should not be included in any xml configuration file.
     "Instantiate" instead the cpuRegistrationAveragingGadget or gpuRegistrationAveragingGadget.
  */
  template<class ARRAY_TYPE, unsigned int D> class EXPORTGADGETS_MOCO RegistrationAveragingGadget 
    : public Gadget2<ISMRMRD::ImageHeader, hoNDArray< typename ARRAY_TYPE::element_type > > // se note below
  {
    //
    // We use hoNDArray to interface the gadget chain, even if ARRAY_TYPE is a cuNDArray
    // Instead of hard coding the interface to use single precision (float), 
    // "typename ARRAY_TYPE::element_type" could in principle denote a double precision type (double) as well.
    // Registration of complex images is however not supported currently...
    //
    
  public:
    
    RegistrationAveragingGadget() {
      this->of_solver_ = 0x0;
      this->number_of_phases_ = 0; // This is a property queried from the PhysioInterpolationGadget
    }

    virtual ~RegistrationAveragingGadget() {
      if( this->of_solver_ ) delete this->of_solver_;
    }

  protected:
    GADGET_PROPERTY(alpha, float, "Alpha regularization parameters", 0.05);
    GADGET_PROPERTY(beta,  float, "Beta regularization parameters", 1.0);
    GADGET_PROPERTY(limit, float, "Iteration limit", 0.01);
    GADGET_PROPERTY(num_multiresolution_levels, int, "Number of multiresolution levels", 3);
    GADGET_PROPERTY(max_iterations_per_level, int, "Maximum number of iterations per level", 500);
    GADGET_PROPERTY(output_convergence, bool, "Output convergence", false);
    GADGET_PROPERTY(phases, unsigned short, "Number of cardiac phases", 30);


    virtual int process_config(ACE_Message_Block *mb)
    {
      this->alpha_ = (typename ARRAY_TYPE::element_type)alpha.value();
      this->beta_  = (typename ARRAY_TYPE::element_type)beta.value();
      this->limit_ = (typename ARRAY_TYPE::element_type)limit.value();
      this->output_convergence_ = output_convergence.value();
      this->num_multires_levels_ = num_multiresolution_levels.value();
      this->max_iterations_per_level_ = max_iterations_per_level.value();

      this->number_of_phases_ = phases.value();
      
      GDEBUG("Configured for %d phases\n", this->number_of_phases_); 
      return GADGET_OK;
    }
    
    virtual int process( GadgetContainerMessage< ISMRMRD::ImageHeader > *m1,
                         GadgetContainerMessage< hoNDArray< typename ARRAY_TYPE::element_type > > *m2 )
    {

      //GDEBUG("\nSERIES: %d, PHASE: %d", m1->getObjectPtr()->image_series_index, m1->getObjectPtr()->phase );

      // If this image header corresponds to series 0, it is not part of the sorted phases.
      // Just pass those images along...
      //

      if( m1->getObjectPtr()->image_series_index < 9 ){
        return this->next()->putq(m1);
      }
      
      // At first pass allocate the image buffer array.
      //
      
      if( this->phase_images_.empty() ){
      
        this->image_dimensions_ = *m2->getObjectPtr()->get_dimensions();
        this->phase_images_ =  decltype(this->phase_images_){};
	
        size_t bsize = sizeof(GadgetContainerMessage<ISMRMRD::ImageHeader>)*100*this->number_of_phases_;
	

        // Setup the optical flow solver
        //
	
        if( this->setup_solver() != GADGET_OK ){
          GDEBUG("Failed to set up optical flow solver\n");
          return GADGET_FAIL;
        }
      }
      
      //
      // Put the incoming images on the appropriate queue (based on the phase index).
      // 
      unsigned int phase = m1->getObjectPtr()->phase;

      this->phase_images_[phase].emplace(m1);

      return GADGET_OK;
    }

    // All the work is done here in the close method
    //

    virtual int close(unsigned long flags)
    {
      if( !this->phase_images_.empty()){
      
        GDEBUG("RegistrationAveragingGadget::close (performing registration and averaging images)\n");
      
        // Make sure we have the same number of images on all phase queues
        // (It doesn't really matter, but if not the case something probably went wrong upstream)
        //

        unsigned int num_images = this->phase_images_[0].size();

        GDEBUG("Number of images for phase 0: %d", num_images );
        
        for( unsigned int phase = 0; phase< this->number_of_phases_; phase++ ){

          unsigned int num_images_phase = this->phase_images_[phase].size();
          GDEBUG("Number of images for phase %d: %d", phase, num_images_phase );

          if( num_images != num_images_phase ){
            GDEBUG("Failed to set up registration, a different number of images received for each phase\n");
            return Gadget::close(flags);
          }
        }
      
        if( num_images == 0 ){
          GDEBUG("No images to register\n");
          return Gadget::close(flags);
        }

        for( unsigned int phase=0; phase < this->number_of_phases_; phase++ ){
	
          unsigned int num_image_elements = this->image_dimensions_[0]*image_dimensions_[1];
          std::vector<size_t> moving_dims = this->image_dimensions_;
          moving_dims.push_back(num_images-1);
	
          GadgetContainerMessage<ISMRMRD::ImageHeader> *header;

          ARRAY_TYPE fixed_image;
          ARRAY_TYPE moving_image(moving_dims);
	
          for( unsigned int image=0; image<num_images; image++ ){
	  
            auto m1 = this->phase_images_[phase].front().release();

            this->phase_images_[phase].pop();

            GadgetContainerMessage< hoNDArray<typename ARRAY_TYPE::element_type> > *m2 = 
              AsContainerMessage< hoNDArray<typename ARRAY_TYPE::element_type> >(m1->cont());
	  
            if( m2 == 0x0 ) {
              GDEBUG("Unexpected continuation on queue\n");
              m1->release();
              return Gadget::close(flags);
            }
	  
            if( image == 0 ){

              // Setup the fixed image.
              // If ARRAY_TYPE is an cuNDArray the following assignment uploads the array to the device,
              // for an 'hoNDArray' it merely copies the array.
              fixed_image = *m2->getObjectPtr();

              // We are going to pass on the averaged image using this header
              header = m1; 

              // The continuation will be a new array (set after registration).
              // No registration is however performed if we received only one image. 
              // In the latter case keep the current continuation.
              if( num_images > 1 ){	      
                m1->cont(0x0); 
                m2->release();
              }
            }
            else{

              // Assign this image as the 'image-1'th frame in the moving image
              ARRAY_TYPE tmp_moving(image_dimensions_, moving_image.get_data_ptr()+(image-1)*num_image_elements);
              tmp_moving = *m2->getObjectPtr(); // Copy as for the fixed image
              m1->release();	    
            }
          }
	
          if( num_images > 1 ){
	  
            // Perform registration for the current phase
            //
	  
            boost::shared_ptr<ARRAY_TYPE> deformations;
            {
              GadgetronTimer timer("Running registration");
              deformations = this->of_solver_->solve( &fixed_image, &moving_image );
            }

            // Deform moving images based on the registration
            //

            boost::shared_ptr<ARRAY_TYPE> deformed_moving;
            {
              GadgetronTimer timer("Applying deformation");
              deformed_moving = this->of_solver_->deform( &moving_image, deformations );
            }

            // Accumulate the deformed moving images (into one image) and add this image to the fixed image. 
            // Then divide by the number of images to get the average.
            //	  

            fixed_image += ((deformed_moving->get_number_of_dimensions() == 3) ? *sum(deformed_moving.get(), 2) : *deformed_moving);
            fixed_image /= ((typename ARRAY_TYPE::element_type)num_images);

            // Pass along averaged image
            //
	  
            if( set_continuation( header, &fixed_image ) < 0 ) {
              GDEBUG("Failed to set continuation\n");
              header->release();
              return Gadget::close(flags);
            }
          }

          if( this->next()->putq(header) < 0 ) {
            GDEBUG("Failed to put registrered image on queue\n");
            header->release();
            return Gadget::close(flags);
          }
        }
      }
    
      return Gadget::close(flags);
    }

    virtual int setup_solver() = 0;
    virtual int set_continuation( GadgetContainerMessage<ISMRMRD::ImageHeader>* m1, ARRAY_TYPE *continuation ) = 0;

  protected:
    opticalFlowSolver<ARRAY_TYPE,D> *of_solver_;
    typename ARRAY_TYPE::element_type alpha_;
    typename ARRAY_TYPE::element_type beta_;
    typename ARRAY_TYPE::element_type limit_;
    bool output_convergence_;
    unsigned int num_multires_levels_;
    unsigned int max_iterations_per_level_;

  private:
    std::map<int, std::queue<std::unique_ptr<GadgetContainerMessage<ISMRMRD::ImageHeader>>>> phase_images_;
    std::vector<size_t> image_dimensions_;
    unsigned short number_of_phases_;    
  };
}

#endif //RegistrationAveragingGadget_H
