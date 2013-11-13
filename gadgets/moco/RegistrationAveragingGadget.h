#ifndef RegistrationAveragingGadget_H
#define RegistrationAveragingGadget_H

#include "Gadget.h"
#include "GadgetMRIHeaders.h"
#include "hoNDArray.h"
#include "complext.h"
#include "PhysioInterpolationGadget.h"
#include "GadgetStreamController.h"
#include "GadgetronTimer.h"
#include "gadgetron_moco_export.h"
#include "hoNDArray_fileio.h"

#include <ismrmrd.h>
#include <complex>
#include <boost/shared_ptr.hpp>
#include <boost/shared_array.hpp>

namespace Gadgetron{  

  template<class ARRAY_TYPE, unsigned int D> class opticalFlowSolver;
  
  /**
     This is an abstract gadget class and consequently should not be included in any xml configuration file.
     "Instatiate" instead the cpuRegistrationAveragingGadget or gpuRegistrationAveragingGadget.
  */
  template<class ARRAY_TYPE, unsigned int D> class EXPORTGADGETS_MOCO RegistrationAveragingGadget 
    : public Gadget2<ISMRMRD::ImageHeader, hoNDArray< typename ARRAY_TYPE::element_type > > // se note below
  {
    //
    // We use hoNDArray to interface the gadget chain, even if ARRAY_TYPE is a cuNDArray
    // Instead of hard coding the interface to use single precision (float), 
    // "typename ARRAY_TYPE::element_type" could in principle denote a double precison type (double) as well.
    // Registration of complex images is however not supported currently...
    //
    
  public:
    
    RegistrationAveragingGadget() {
      this->of_solver_ = 0x0;
      this->number_of_phases_ = 0; // This is a property queried from the PhysioInterpolationGadget
      this->set_parameter(std::string("alpha").c_str(), "0.05");
      this->set_parameter(std::string("beta").c_str(), "1.0");
      this->set_parameter(std::string("limit").c_str(), "0.01");
      this->set_parameter(std::string("num_multiresolution_levels").c_str(), "3");
      this->set_parameter(std::string("max_iterations_per_level").c_str(), "500");    
      this->set_parameter(std::string("output_convergence").c_str(), "false");
    }

    virtual ~RegistrationAveragingGadget() {
      if( this->of_solver_ ) delete this->of_solver_;
    }

  protected:

    virtual int process_config(ACE_Message_Block *mb)
    {
      this->alpha_ = (typename ARRAY_TYPE::element_type)this->get_double_value("alpha");
      this->beta_  = (typename ARRAY_TYPE::element_type)this->get_double_value("beta");
      this->limit_ = (typename ARRAY_TYPE::element_type)this->get_double_value("limit");
      this->output_convergence_ = this->get_bool_value(std::string("output_convergence").c_str());
      this->num_multires_levels_ = this->get_int_value(std::string("num_multiresolution_levels").c_str());
      this->max_iterations_per_level_ = this->get_int_value(std::string("max_iterations_per_level").c_str());
      
      // Fow now we require the existence of a gadget named "PhysioInterpolationGadget" upstream,
      // to determine the number of incoming phases.
      //
      
      GadgetStreamController *controller = this->get_controller();
    
      if( controller == 0x0 ){
        GADGET_DEBUG1("Failed to get controller\n");
        return GADGET_FAIL;
      }
      
      PhysioInterpolationGadget *physio = 
        dynamic_cast<PhysioInterpolationGadget*>( controller->find_gadget(std::string("PhysioInterpolationGadget")) );
      
      if( physio == 0x0 ){
        GADGET_DEBUG1("Could not find (or cast) PhysioInterpolationGadget in gadget stream\n");
        return GADGET_FAIL;
      }
      
      this->number_of_phases_ = physio->get_number_of_phases();      
      
      GADGET_DEBUG2("Configured for %d phases\n", this->number_of_phases_); 
      return GADGET_OK;
    }
    
    virtual int process( GadgetContainerMessage< ISMRMRD::ImageHeader > *m1,
                         GadgetContainerMessage< hoNDArray< typename ARRAY_TYPE::element_type > > *m2 )
    {

      //GADGET_DEBUG2("\nSERIES: %d, PHASE: %d", m1->getObjectPtr()->image_series_index, m1->getObjectPtr()->phase );

      // If this image header corresponds to series 0, it is not part of the sorted phases.
      // Just pass those images along...
      //

      if( m1->getObjectPtr()->image_series_index < 9 ){
        return this->next()->putq(m1);
      }
      
      // At first pass allocate the image buffer array.
      //
      
      if( this->phase_images_.get() == 0x0 ){
      
        this->image_dimensions_ = *m2->getObjectPtr()->get_dimensions();
        this->phase_images_ = boost::shared_array< ACE_Message_Queue<ACE_MT_SYNCH> >
          (new ACE_Message_Queue<ACE_MT_SYNCH>[this->number_of_phases_]);      
	
        size_t bsize = sizeof(GadgetContainerMessage<ISMRMRD::ImageHeader>)*100*this->number_of_phases_;
	
        for( unsigned int i=0; i<this->number_of_phases_; i++ ){
          this->phase_images_[i].high_water_mark(bsize);
          this->phase_images_[i].low_water_mark(bsize);      
        }
	
        // Setup the optical flow solver
        //
	
        if( this->setup_solver() != GADGET_OK ){
          GADGET_DEBUG1("Failed to set up optical flow solver\n");
          return GADGET_FAIL;
        }
      }
      
      //
      // Put the incoming images on the appropriate queue (based on the phase index).
      // 
      
      unsigned int phase = m1->getObjectPtr()->phase;
      
      if( this->phase_images_[phase].enqueue_tail(m1) < 0 ) {
        GADGET_DEBUG1("Failed to add image to buffer\n");
        return GADGET_FAIL;
      }

      return GADGET_OK;
    }

    // All the work is done here in the close method
    //

    virtual int close(unsigned long flags)
    {
      if( this->phase_images_.get() ){
      
        GADGET_DEBUG1("RegistrationAveragingGadget::close (performing registration and averaging images)\n");
      
        // Make sure we have the same number of images on all phase queues
        // (It doesn't really matter, but if not the case something probably went wrong upstream)
        //

        unsigned int num_images = this->phase_images_[0].message_count();

        GADGET_DEBUG2("Number of images for phase 0: %d", num_images );
        
        for( unsigned int phase = 0; phase< this->number_of_phases_; phase++ ){

          unsigned int num_images_phase = this->phase_images_[phase].message_count();
          GADGET_DEBUG2("Number of images for phase %d: %d", phase, num_images_phase );

          if( num_images != num_images_phase ){
            GADGET_DEBUG1("Failed to set up registration, a different number of images received for each phase\n");
            return Gadget::close(flags);
          }
        }
      
        if( num_images == 0 ){
          GADGET_DEBUG1("No images to register\n");
          return Gadget::close(flags);
        }

        for( unsigned int phase=0; phase < this->number_of_phases_; phase++ ){
	
          unsigned int num_image_elements = this->image_dimensions_[0]*image_dimensions_[1];
          std::vector<unsigned int> moving_dims = this->image_dimensions_;
          moving_dims.push_back(num_images-1);
	
          GadgetContainerMessage<ISMRMRD::ImageHeader> *header;

          ARRAY_TYPE fixed_image;
          ARRAY_TYPE moving_image(&moving_dims);
	
          for( unsigned int image=0; image<num_images; image++ ){
	  
            ACE_Message_Block *mbq;
	  
            if( this->phase_images_[phase].dequeue_head(mbq) < 0 ) {
              GADGET_DEBUG1("Image header dequeue failed\n");
              return Gadget::close(flags);
            }
	  
            GadgetContainerMessage<ISMRMRD::ImageHeader> *m1 = AsContainerMessage<ISMRMRD::ImageHeader>(mbq);
	  
            if( m1 == 0x0 ) {
              GADGET_DEBUG1("Unexpected image type on queue\n");
              return Gadget::close(flags);
            }
	  
            GadgetContainerMessage< hoNDArray<typename ARRAY_TYPE::element_type> > *m2 = 
              AsContainerMessage< hoNDArray<typename ARRAY_TYPE::element_type> >(m1->cont());
	  
            if( m2 == 0x0 ) {
              GADGET_DEBUG1("Unexpected continuation on queue\n");
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
              ARRAY_TYPE tmp_moving(&image_dimensions_, moving_image.get_data_ptr()+(image-1)*num_image_elements);
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
	  
            /*{
            // The debug code below only compiles for cuNDArrays.
            // To use (temporarily) comment out
            // list(APPEND CPU_GADGETS cpuRegistrationAveragingGadget.cpp)
            // in the CMakeList.txt
            //
            char filename[256];
            sprintf((char*)filename, "fixed_%d.real", phase);
            write_nd_array<float>( fixed_image.to_host().get(), filename );
            sprintf((char*)filename, "moving_%d.real", phase);
            write_nd_array<float>( moving_image.to_host().get(), filename );
            sprintf((char*)filename, "deformed_moving_%d.real", phase);
            write_nd_array<float>( deformed_moving->to_host().get(), filename );
            sprintf((char*)filename, "deformation_%d.real", phase);
            write_nd_array<float>( deformations->to_host().get(), filename );
            } */

	 
            // Accumulate the deformed moving images (into one image) and add this image to the fixed image. 
            // Then divide by the number of images to get the average.
            //	  
	  
            fixed_image += ((deformed_moving->get_number_of_dimensions() == 3) ? *sum(deformed_moving.get(), 2) : *deformed_moving);
            fixed_image /= ((typename ARRAY_TYPE::element_type)num_images);
	  
            // Pass along averaged image
            //
	  
            if( set_continuation( header, &fixed_image ) < 0 ) {
              GADGET_DEBUG1("Failed to set continuation\n");
              header->release();
              return Gadget::close(flags);
            }
          }

          if( this->next()->putq(header) < 0 ) {
            GADGET_DEBUG1("Failed to put registrered image on queue\n");
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
    boost::shared_array< ACE_Message_Queue<ACE_MT_SYNCH> > phase_images_;
    std::vector<unsigned int> image_dimensions_;
    unsigned short number_of_phases_;    
  };
}

#endif //RegistrationAveragingGadget_H
