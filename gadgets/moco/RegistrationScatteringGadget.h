#ifndef RegistrationScatteringGadget_H
#define RegistrationScatteringGadget_H

#include "Gadget.h"
#include "GadgetMRIHeaders.h"
#include "hoNDArray.h"
#include "complext.h"
#include "PhysioInterpolationGadget.h"
#include "GadgetronTimer.h"
#include "gadgetron_moco_export.h"
#include "hoNDArray_fileio.h"

#include <ismrmrd/ismrmrd.h>
#include <complex>
#include <queue>

namespace Gadgetron{  

  template<class ARRAY_TYPE, unsigned int D> class opticalFlowSolver;
  
  /**
     This is an abstract gadget class and consequently should not be included in any xml configuration file.
     Use instead the gpuRegistrationScatteringGadget.
  */
  template<class ARRAY_TYPE, unsigned int D> class EXPORTGADGETS_MOCO RegistrationScatteringGadget 
    : public Gadget2<ISMRMRD::ImageHeader, hoNDArray< typename ARRAY_TYPE::element_type > > // se note below
  {
    //
    // We use hoNDArray to interface the gadget chain, even if ARRAY_TYPE is a cuNDArray
    // Instead of hard coding the interface to use single precision (float), 
    // "typename ARRAY_TYPE::element_type" could in principle denote a double precision type (double) as well.
    // Registration of complex images is however not supported currently.
    //
    
  public:
    
    RegistrationScatteringGadget() {
      this->of_solver_ = 0x0;
      this->number_of_phases_ = 0; // This is a property queried from the PhysioInterpolationGadget
    }

    virtual ~RegistrationScatteringGadget() {
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
      
      if( this->image_dimensions_.empty() ){
      
        this->image_dimensions_ = *m2->getObjectPtr()->get_dimensions();
        this->phase_images_ = decltype(this->phase_images_){};

        
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
      if( !this->phase_images_.empty() ){
      
        GDEBUG("RegistrationScatteringGadget::close (performing registration and scattering images)\n");
      
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

        // These are the dimensions of the vector field written out
        // - just a plain 'write_nd_array' below for now...
        //

        std::vector<size_t> reg_dims = this->image_dimensions_; // x,y
        reg_dims.push_back(num_images-1); // this many registrations 
        reg_dims.push_back(2); // 2d flow vectors
        reg_dims.push_back(this->number_of_phases_);
        ARRAY_TYPE reg_field(&reg_dims);
        unsigned int num_reg_elements_phase = reg_dims[0]*reg_dims[1]*reg_dims[2]*reg_dims[3];

        for( unsigned int phase=0; phase < this->number_of_phases_; phase++ ){
	
          unsigned int num_image_elements = this->image_dimensions_[0]*image_dimensions_[1];
          std::vector<size_t> fixed_dims = this->image_dimensions_;
          fixed_dims.push_back(num_images-1);
	
          std::vector< GadgetContainerMessage<ISMRMRD::ImageHeader>*> headers;

          ARRAY_TYPE fixed_image(&fixed_dims);
          ARRAY_TYPE moving_image;
	
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

              // Setup the moving image.
              // If ARRAY_TYPE is an cuNDArray the following assignment uploads the array to the device,
              // for an 'hoNDArray' it merely copies the array.
              //

              moving_image = *m2->getObjectPtr();
              headers.push_back(m1);
            }
            else{

              // Assign this image as the 'image-1'th frame in the moving image
              //

              ARRAY_TYPE tmp_fixed(&image_dimensions_, fixed_image.get_data_ptr()+(image-1)*num_image_elements);
              tmp_fixed = *m2->getObjectPtr(); // Copy as for the moving image
              headers.push_back(m1);

              // The continuation will be a new array (set after registration).
              // No registration is however performed if we received only one image. 
              // In the latter case keep the current continuation.
              //

              if( num_images > 1 ){
                m1->cont(0x0);
                m2->release();
              }             
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

            // Copy displacement field to vector field array
            //
            
            {              
              std::vector<size_t> phase_reg_dims = reg_dims; phase_reg_dims.pop_back();
              ARRAY_TYPE tmp_in( &phase_reg_dims, deformations->get_data_ptr() ); // the vector field has an extra dimension for CK (to be discarded)
              ARRAY_TYPE tmp_out( &phase_reg_dims, reg_field.get_data_ptr()+phase*num_reg_elements_phase );
              tmp_out = tmp_in;
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
            // list(APPEND CPU_GADGETS cpuRegistrationScatteringGadget.cpp)
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


            // Pass along the deformed moving images
            //	  
	  
            for( unsigned int i=0; i<headers.size(); i++ ){
              
              if( i==0 ){
                GDEBUG("Putting image %d image on queue\n", i);
                
                if( this->next()->putq(headers[i]) < 0 ) {
                  GDEBUG("Failed to put registrered image on queue\n");
                  headers[i]->release();
                  return Gadget::close(flags);
                }
              }
              else{                
                std::vector<size_t> moving_dims = *moving_image.get_dimensions();
                cuNDArray<float> subimage( &moving_dims, deformed_moving->get_data_ptr()+(i-1)*num_image_elements);
                
                if( set_continuation( headers[i], &subimage ) < 0 ) {
                  GDEBUG("Failed to set continuation\n");
                  headers[i]->release();
                  return Gadget::close(flags);
                }
                
                GDEBUG("Putting image %d image on queue\n", i);
                
                if( this->next()->putq(headers[i]) < 0 ) {
                  GDEBUG("Failed to put registrered image on queue\n");
                  headers[i]->release();
                  return Gadget::close(flags);
                }
              }
            }
          }
        }
        
        // Write out the result after permutation to the data order
        // - to be betetr suited for a subsequent reconstruction pass
        //
        
        std::vector<size_t> order;
        order.push_back(0); 
        order.push_back(1);
        order.push_back(4);
        order.push_back(2);
        order.push_back(3);
        
        GDEBUG("Writing out displacement field with dimensions: %d %d %d %d %d\n", order[0], order[1], order[2], order[3], order[4]);
        auto displacement_field = permute(reg_field,order);
        write_displacement_field( &displacement_field );

      }
      
      return Gadget::close(flags);
    }
    
    virtual int setup_solver() = 0;
    virtual int set_continuation( GadgetContainerMessage<ISMRMRD::ImageHeader>* m1, ARRAY_TYPE *continuation ) = 0;
    virtual int write_displacement_field( ARRAY_TYPE *vec_field ) = 0;
    
  protected:
    opticalFlowSolver<ARRAY_TYPE,D> *of_solver_;
    typename ARRAY_TYPE::element_type alpha_;
    typename ARRAY_TYPE::element_type beta_;
    typename ARRAY_TYPE::element_type limit_;
    bool output_convergence_;
    unsigned int num_multires_levels_;
    unsigned int max_iterations_per_level_;
    
  private:
    std::map<int,std::queue<std::unique_ptr<GadgetContainerMessage<ISMRMRD::ImageHeader>>>> phase_images_;
    std::vector<size_t> image_dimensions_;
    unsigned short number_of_phases_;    
  };
}

#endif //RegistrationScatteringGadget_H
