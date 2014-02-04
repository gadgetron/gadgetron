#include "CMRT3DGadget.h"
#include "GadgetIsmrmrdReadWrite.h"
#include "hoNDFFT.h"
#include "hoNDArray_utils.h"
#include "hoNDArray_fileio.h"
#include "cuNDArray_utils.h"
#include "cuNDArray_elemwise.h"
#include "cuNDArray_operators.h"
#include "radial_utilities.h"
#include "vector_td_operators.h"

static const float alpha_ = 2.0f; // oversampling for radial NFFT. If to be changed, also change setup arguments
static const float W_ = 5.5f; // Kaiser-Bessel Window size for the radial NFFT
static const float readout_oversampling_factor_ = 1.0f; // There is no "readout" oversampling for the radial NFFT

namespace Gadgetron {

  /**
   *   Expects ISMRMRD XML configuration
   *
   */

  int CMRT3DGadget::process_config(ACE_Message_Block* mb)
  {
    boost::shared_ptr<ISMRMRD::ismrmrdHeader> cfg = parseIsmrmrdXMLHeader(std::string(mb->rd_ptr()));
  
    ISMRMRD::ismrmrdHeader::encoding_sequence enc = cfg->encoding();
  
    if (enc.size() != 1) {
      GADGET_DEBUG2("Number of encoding spaces: %d\n", enc.size());
      GADGET_DEBUG1("This simple CMRT3DGadget only supports one encoding space\n");
      return GADGET_FAIL;
    }

    ISMRMRD::encodingSpaceType e_space = (*enc.begin()).encodedSpace();
    ISMRMRD::encodingSpaceType r_space = (*enc.begin()).reconSpace();
    ISMRMRD::encodingLimitsType e_limits = (*enc.begin()).encodingLimits();
  
    // Matrix size x is the oversampled readout size : size FH*2 -- not for the hyperplarization?
    // Matrix size y is the AP/RL size

    image_space_dimensions_3D_.push_back(e_space.matrixSize().y());
    image_space_dimensions_3D_.push_back(e_space.matrixSize().y());
    image_space_dimensions_3D_.push_back(e_space.matrixSize().x()/*/2*/);

    GADGET_DEBUG2("Matrix size: %d, %d, %d\n", 
                  image_space_dimensions_3D_[0],
                  image_space_dimensions_3D_[1],
                  image_space_dimensions_3D_[2] );
  
    num_projections_expected_ = get_int_value(std::string("projections_per_recon").c_str());
    projections_percentage_ = get_int_value(std::string("projections_percentage").c_str());
    num_projections_to_use_ = num_projections_expected_/(100/projections_percentage_);

    GADGET_DEBUG2("Number of projections (expected/utilization percentage): %d/%d\n", num_projections_expected_, projections_percentage_ );
    GADGET_DEBUG2("I.e. using %d projections for the reconstruction\n", num_projections_to_use_ );
  
    std::vector<size_t> dims;
    dims.push_back(image_space_dimensions_3D_[0]); // number of samples per radial profile
    dims.push_back(num_projections_to_use_);       // number of radial profiles
    dims.push_back(image_space_dimensions_3D_[2]); // number of slices

    buffer_ = boost::shared_ptr< cuNDArray< complext<float> > >( new cuNDArray< complext<float> >(&dims) );

    // Calculate trajectories and dcw for the radial NFFTs
    //
  
    boost::shared_ptr< cuNDArray<floatd2> > traj = calculate_trajectory();
    boost::shared_ptr< cuNDArray<float> > dcw = calculate_density_compensation();
  

    if( !traj.get() || !dcw.get() ){
      GADGET_DEBUG1("Failed to initialize radial trajecotory/dcw\n");
      return GADGET_FAIL;
    }

    // Setup radial NFFT encoding operator
    // 

    E_ = boost::shared_ptr< cuNFFTOperator<float,2> >( new cuNFFTOperator<float,2>() );

    E_->set_dcw( dcw );

    E_->setup( uint64d2(image_space_dimensions_3D_[0], image_space_dimensions_3D_[1]),
               uint64d2(image_space_dimensions_3D_[0], image_space_dimensions_3D_[1])<<1, // !! <-- alpha_
               W_ );
    
    E_->preprocess( traj.get() );

    return GADGET_OK;
  }

  int CMRT3DGadget::process(GadgetContainerMessage<ISMRMRD::ImageHeader>* m1,
                            GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2)
  {
    // Check if we should ignore this image
    // - to simulate undersampling in the number of slices
  
    if( (images_received_%(100/projections_percentage_)) != 0 ){
      // Ignore this image
      images_received_++;
      return GADGET_OK;
    }

    // We will not pass along m2, so we can modify its array safely
    //

    hoNDArray< std::complex<float> > *host_image = m2->getObjectPtr();
  
    // Some validity checks
    //
  
    if( !(host_image->get_number_of_dimensions() == 2 || ( host_image->get_number_of_dimensions()==3 && host_image->get_size(2)==1 ))) {
      GADGET_DEBUG2("The input image has an unexpected number of dimensions\n", host_image->get_number_of_dimensions());
      return GADGET_FAIL;
    }
    
    if( host_image->get_size(0) != image_space_dimensions_3D_[0] || 
        host_image->get_size(1) != image_space_dimensions_3D_[2] ){
      GADGET_DEBUG2("The input image has unexpected dimensionality: %d %d %d %d\n", host_image->get_size(0), host_image->get_size(1), image_space_dimensions_3D_[0], image_space_dimensions_3D_[1] );
      return GADGET_FAIL;
    }

    // Perform batched 1D FFTs along the phase encoding direction of the input image 
    // I.e. permute and then perform FFT

    std::vector<size_t> order;
    order.push_back(1); order.push_back(0);
  
    *host_image = *permute( host_image, &order );
    hoNDFFT<float>::instance()->fft( host_image, 0 );

    // Next copy each line into the buffer 
    //

    GADGET_DEBUG2("Received image #%d\n", images_received_);    

    for( size_t row=0;row<host_image->get_size(1); row++ ){

      size_t offset_in = 
        row*host_image->get_size(0);
    
      size_t offset_out = 
        row*host_image->get_size(0)*num_projections_to_use_+
        images_used_*host_image->get_size(0);
    
      if( cudaMemcpy( buffer_->get_data_ptr()+offset_out, 
                      host_image->get_data_ptr()+offset_in, 
                      host_image->get_size(0)*sizeof(complext<float>),
                      cudaMemcpyHostToDevice ) != cudaSuccess ){
        GADGET_DEBUG2("Upload to device for line %d failed\n", row);
        return GADGET_FAIL;
      }
    }
  
    // Another image has been received and uploaded...
    //

    images_received_++;
    images_used_++;

    // When we are ready to perform reconstruction, do it...
    //

    if( images_used_ == num_projections_to_use_ ){

      GADGET_DEBUG1("\n\nPerforming reconstruction\n");
    
      std::vector<size_t> dims;
      dims.push_back(image_space_dimensions_3D_[0]);
      dims.push_back(image_space_dimensions_3D_[1]);
      dims.push_back(image_space_dimensions_3D_[2]);
    
      cuNDArray< complext<float> > result(&dims); 
    
      E_->mult_MH( buffer_.get(), &result );

      /*
        boost::shared_ptr< hoNDArray<complext<float> > > host_result = result.to_host();
        write_nd_array<complext<float> >(host_result.get(), "result.cplx");   

        boost::shared_ptr< hoNDArray<float> > host_norm = abs(&result)->to_host();
        write_nd_array<float>( host_norm.get(), "result.real" );*/
      
      // Create new image header/image to pass along
      //

      GadgetContainerMessage<ISMRMRD::ImageHeader> *m = 
        new GadgetContainerMessage<ISMRMRD::ImageHeader>();
    
      GadgetContainerMessage< hoNDArray< std::complex<float> > > *cm = 
        new GadgetContainerMessage< hoNDArray< std::complex<float> > >();      
    
      *m->getObjectPtr() = *m1->getObjectPtr();
      m->cont(cm);
    
      // std::complex<float> and Gadgetron::complext<float> are binary compatible
      boost::shared_ptr< hoNDArray< complext<float> > > host_result = result.to_host();
      *cm->getObjectPtr() = *((hoNDArray< std::complex<float> >*) host_result.get());
    
      m->getObjectPtr()->matrix_size[0] = dims[0];
      m->getObjectPtr()->matrix_size[1] = dims[1];
      m->getObjectPtr()->matrix_size[2] = dims[2];
      m->getObjectPtr()->channels       = 1;
      m->getObjectPtr()->image_index    = 1;
    
      if (this->next()->putq(m) < 0) {
        GADGET_DEBUG1("Failed to put result image on to queue\n");
        m->release();
        return GADGET_FAIL;
      }
    }
    
    m1->release();
    return GADGET_OK;
  }

  boost::shared_ptr< cuNDArray<floatd2> >
  CMRT3DGadget::calculate_trajectory()
  {
    // Define trajectories
    boost::shared_ptr< cuNDArray<floatd2> > traj = 
      compute_radial_trajectory_fixed_angle_2d<float>
      ( image_space_dimensions_3D_[0], num_projections_to_use_, 1 /*number of frames*/ );
  
    if (!traj.get()) {
      GADGET_DEBUG1("Failed to compute radial trajectory");
      return boost::shared_ptr< cuNDArray<floatd2> >();
    }
  
    return traj;
  }

  boost::shared_ptr< cuNDArray<float> >
  CMRT3DGadget::calculate_density_compensation()
  { 
    // Compute density compensation weights
    boost::shared_ptr< cuNDArray<float> > dcw = 
      compute_radial_dcw_fixed_angle_2d
      ( image_space_dimensions_3D_[0], num_projections_to_use_, alpha_, 1.0f/readout_oversampling_factor_ );
  
    if (!dcw.get()) {
      GADGET_DEBUG1("Failed to compute density compensation weights\n");
      return boost::shared_ptr< cuNDArray<float> >();
    }
  
    return dcw;
  }

  GADGET_FACTORY_DECLARE(CMRT3DGadget)
}
