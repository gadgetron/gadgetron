#include "Gadgetron.h"
#include "SpiralGadget.h"
#include "GadgetXml.h"
#include "hoNDArray_fileio.h"
#include "cuNDArray.h"
#include "vector_td.h"
#include "vector_td_operators.h"
#include "check_CUDA.h"

#include <vector>

void calc_vds(double slewmax,double gradmax,double Tgsample,double Tdsample,int Ninterleaves,
			  double* fov, int numfov,double krmax,
			  int ngmax, double** xgrad,double** ygrad,int* numgrad);

void calc_traj(double* xgrad, double* ygrad, int ngrad, int Nints, double Tgsamp, double krmax,
			  double** x_trajectory, double** y_trajectory,
              double** weights);


SpiralGadget::SpiralGadget()
  : samples_to_skip_start_(0)
  , samples_to_skip_end_(0)
  , samples_per_adc_(0)
  , adcs_per_interleave_(0)
  , host_data_buffer_(0)
{

}
 
SpiralGadget::~SpiralGadget()
{
  if (host_data_buffer_) delete [] host_data_buffer_;
}

int SpiralGadget::process_config(ACE_Message_Block* mb)
{

  TiXmlDocument doc;
  doc.Parse(mb->rd_ptr());

  //GADGET_DEBUG1("Calculating trajectory\n");

  double  smax = GetDoubleParameterValueFromXML(&doc, "spiral", "MaxSlewRate_Gcms");   /*	Maximum slew rate, G/cm/s		*/
  double  gmax = GetDoubleParameterValueFromXML(&doc, "spiral", "MaxGradient_Gcm");   /* 	maximum gradient amplitude, G/cm	*/
  double  Tdsamp = GetIntParameterValueFromXML(&doc, "spiral", "SamplingTime_ns")/(1.0e9);     /*	Data Sample period (s)			*/
  //double  Tgsamp = 1e-5;     /*	Data Sample period (s)			*/
  int     Nints = GetIntParameterValueFromXML(&doc, "spiral", "Interleaves");  	/*	Number of interleaves			*/
  double  fov = GetDoubleParameterValueFromXML(&doc, "spiral", "FOVCoeff_1");	/*	FOV coefficients		*/
  int     Nfov = 1;       /*	Number of FOV coefficients		*/
  double  krmax = GetDoubleParameterValueFromXML(&doc, "spiral", "krmax_cm");  	/*	Maximum k-space extent (/cm)		*/
  int     ngmax = 1e5;      /*	Maximum number of gradient samples	*/
  double  *xgrad; 	/* 	X-component of gradient.	*/
  double  *ygrad;     /*	Y-component of gradient.	*/
  double  *x_trajectory;
  double  *y_trajectory;
  double  *weighting;
  int     ngrad;	
  //int     count;

  samples_to_skip_start_  = GetIntParameterValueFromXML(&doc, "spiral", "SamplesToSkipStart");
  samples_to_skip_end_    = GetIntParameterValueFromXML(&doc, "spiral", "SamplesToSkipEnd");
  samples_per_adc_        = GetIntParameterValueFromXML(&doc, "spiral", "SamplesPerADC");
  adcs_per_interleave_    = GetIntParameterValueFromXML(&doc, "spiral", "ADCsPerInterleave");

  image_dimensions_.push_back(GetIntParameterValueFromXML(&doc, "encoding", "matrix_x"));
  image_dimensions_.push_back(GetIntParameterValueFromXML(&doc, "encoding", "matrix_y"));

  /*
  GADGET_DEBUG2("smax:                    %f\n", smax);
  GADGET_DEBUG2("gmax:                    %f\n", gmax);
  GADGET_DEBUG2("Tdsamp:                  %f\n", Tdsamp);
  GADGET_DEBUG2("Nints:                   %d\n", Nints);
  GADGET_DEBUG2("fov:                     %f\n", fov);
  GADGET_DEBUG2("krmax:                   %f\n", krmax);
  GADGET_DEBUG2("samples_to_skip_start_ : %d\n", samples_to_skip_start_);
  GADGET_DEBUG2("samples_to_skip_end_   : %d\n", samples_to_skip_end_);
  GADGET_DEBUG2("samples_per_adc_       : %d\n", samples_per_adc_);
  GADGET_DEBUG2("adcs_per_interleave_   : %d\n", adcs_per_interleave_);
  GADGET_DEBUG2("matrix_size_x          : %d\n", image_dimensions_[0]);
  GADGET_DEBUG2("matrix_size_y          : %d\n", image_dimensions_[1]);
  */


  /*	Call C-function Here to calculate gradients */
  calc_vds(smax,gmax,Tdsamp,Tdsamp,Nints,&fov,Nfov,krmax,ngmax,&xgrad,&ygrad,&ngrad);

  //GADGET_DEBUG2("ngrad (before adjust)   : %d\n", ngrad);
  if ((adcs_per_interleave_*samples_per_adc_-ngrad) != samples_to_skip_end_) {
    ngrad = (adcs_per_interleave_*samples_per_adc_-samples_to_skip_end_) < ngrad ? adcs_per_interleave_*samples_per_adc_-samples_to_skip_end_ : ngrad;
    samples_to_skip_end_ = (adcs_per_interleave_*samples_per_adc_-ngrad);
  }
  //GADGET_DEBUG2("ngrad (after adjust)   : %d\n", ngrad);


  /* Calcualte the trajectory and weights*/
  calc_traj(xgrad, ygrad, ngrad, Nints, Tdsamp, krmax, &x_trajectory, &y_trajectory, &weighting);
  
  host_traj_ = boost::shared_ptr< hoNDArray<floatd2::Type> >(new hoNDArray<floatd2::Type>);
  host_weights_ = boost::shared_ptr< hoNDArray<float> >(new hoNDArray<float>);

  std::vector<unsigned int> trajectory_dimensions;
  trajectory_dimensions.push_back(ngrad*Nints);

  if (!host_traj_->create(&trajectory_dimensions)) {
    GADGET_DEBUG1("Unable to allocate memory for trajectory\n");
    return GADGET_FAIL;
  }
  
  if (!host_weights_->create(&trajectory_dimensions)) {
    GADGET_DEBUG1("Unable to allocate memory for weights\n");
    return GADGET_FAIL;
  }
  

  float* co_ptr = reinterpret_cast<float*>(host_traj_->get_data_ptr());
  float* we_ptr =  reinterpret_cast<float*>(host_weights_->get_data_ptr());

  for (int i = 0; i < (ngrad*Nints); i++) {
    co_ptr[i*2]   = -x_trajectory[i]/2;
    co_ptr[i*2+1] = -y_trajectory[i]/2;
    we_ptr[i] = weighting[i];
  }

  delete [] xgrad;
  delete [] ygrad;
  delete [] x_trajectory;
  delete [] y_trajectory;
  delete [] weighting;

  unsigned int slices = GetIntParameterValueFromXML(&doc, "encoding", "slices");

  std::vector<unsigned int> data_dimensions;
  data_dimensions.push_back(ngrad*Nints);
  data_dimensions.push_back(GetIntParameterValueFromXML(&doc, "encoding", "channels"));

  host_data_buffer_ = new hoNDArray<float_complext::Type>[slices];
  if (!host_data_buffer_) {
    GADGET_DEBUG1("Unable to allocate array for host data buffer\n");
    return GADGET_FAIL;
  }

  for (unsigned int i = 0; i < slices; i++) {
    if (!host_data_buffer_[i].create(&data_dimensions)) {
      GADGET_DEBUG1("Unable to allocate memory for data buffer\n");
      return GADGET_FAIL;
    }
  }


  //Make NFFT plan
  // Matrix sizes
  uintd2 matrix_size = uintd2(image_dimensions_[0],image_dimensions_[1]);
  uintd2 matrix_size_os = uintd2(image_dimensions_[0]*2,image_dimensions_[1]*2);
  
  // Kernel width
  float W = 5.5f;
  
  // Upload host arrays to device arrays
  cuNDArray<floatd2::Type> traj(host_traj_.get());
  gpu_weights_ = cuNDArray<float>(host_weights_.get());
    
  // Initialize plan
  // NFFT_plan<float, 2> plan( matrix_size, matrix_size_os, W );
  plan_ = NFFT_plan<float, 2>( matrix_size, matrix_size_os, W );

  // Preprocess
  bool success = plan_.preprocess( &traj, NFFT_plan<float,2>::NFFT_PREP_ALL );
  
  if (!success) {
    GADGET_DEBUG1("NFFT preprocess failed\n");
    return GADGET_FAIL;
  }

  return GADGET_OK;
}

int SpiralGadget::
process(GadgetContainerMessage<GadgetMessageAcquisition>* m1,
	GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2)
{
 
  unsigned int samples_to_copy = m1->getObjectPtr()->samples-samples_to_skip_end_;
  unsigned int interleave = m1->getObjectPtr()->idx.line;
  unsigned int slice = m1->getObjectPtr()->idx.slice;

  unsigned int samples_per_channel =  host_data_buffer_->get_size(0);

  std::complex<float>* data_ptr    = reinterpret_cast< std::complex<float>* >(host_data_buffer_[slice].get_data_ptr());
  std::complex<float>* profile_ptr = m2->getObjectPtr()->get_data_ptr();

  for (unsigned int c = 0; c < m1->getObjectPtr()->channels; c++) {
    memcpy(data_ptr+c*samples_per_channel+interleave*samples_to_copy,
	   profile_ptr+c*m1->getObjectPtr()->samples, samples_to_copy*sizeof(std::complex<float>));
  }

  if (m1->getObjectPtr()->flags & GADGET_FLAG_LAST_ACQ_IN_SLICE) {
    //GADGET_DEBUG1("Las scan in slice\n");

    unsigned int num_batches = m1->getObjectPtr()->channels;
    
    cuNDArray<float_complext::Type> data(&host_data_buffer_[slice]);
 
    // Setup result image
    std::vector<unsigned int> image_dims;
    image_dims.push_back(image_dimensions_[0]);
    image_dims.push_back(image_dimensions_[1]); 
    image_dims.push_back(num_batches);
    cuNDArray<float_complext::Type> image; image.create(&image_dims);
   
    bool  success = plan_.compute( &data, &image, &gpu_weights_, NFFT_plan<float,2>::NFFT_BACKWARDS );
    if (!success) {
      GADGET_DEBUG1("NFFT compute failed\n");
      return GADGET_FAIL;
    }

    boost::shared_ptr< hoNDArray<float_complext::Type> > image_host = image.to_host();

    GadgetContainerMessage< hoNDArray< std::complex<float> > >* m4 = 
      new GadgetContainerMessage< hoNDArray< std::complex<float> > >();

    if (!m4->getObjectPtr()->create(&image_dimensions_)) {
      GADGET_DEBUG1("Unable to allocate memory for combined image\n"); 
      m4->release();
      return GADGET_FAIL;
    }
    
    unsigned int npixels = image_dimensions_[0]*image_dimensions_[1];
    std::complex<float>* recon_ptr    = reinterpret_cast< std::complex<float>* >(image_host->get_data_ptr());
    std::complex<float>* comb_ptr     = reinterpret_cast< std::complex<float>* >(m4->getObjectPtr()->get_data_ptr());

    for (unsigned int i = 0; i < npixels; i++) {
      float mag = 0.0;
      float phase = 0.0;
      for (unsigned int c = 0; c < num_batches; c++) {
	float mag_tmp = norm(recon_ptr[c*npixels+i]);
	phase += mag_tmp*arg(recon_ptr[c*npixels+i]);
	mag += mag_tmp;
      }
      comb_ptr[i] = std::polar((float)sqrt(mag),phase)*std::complex<float>(npixels,0.0);
    }


    GadgetContainerMessage<GadgetMessageImage>* m3 = 
       new GadgetContainerMessage<GadgetMessageImage>();
    
    m3->cont(m4);
    
    m3->getObjectPtr()->matrix_size[0] = image_dimensions_[0];
    m3->getObjectPtr()->matrix_size[1] = image_dimensions_[1];
    m3->getObjectPtr()->matrix_size[2] = 1;
    m3->getObjectPtr()->channels       = 1;
    m3->getObjectPtr()->data_idx_min       = m1->getObjectPtr()->min_idx;
    m3->getObjectPtr()->data_idx_max       = m1->getObjectPtr()->max_idx;
    m3->getObjectPtr()->data_idx_current   = m1->getObjectPtr()->idx;	

    memcpy(m3->getObjectPtr()->position,m1->getObjectPtr()->position,
	   sizeof(float)*3);

    memcpy(m3->getObjectPtr()->quarternion,m1->getObjectPtr()->quarternion,
	   sizeof(float)*4);
 
    if (this->next()->putq(m3) < 0) {
      m3->release();
      return GADGET_FAIL;
    }

  }

  m1->release();
  return GADGET_OK;
}


GADGET_FACTORY_DECLARE(SpiralGadget)
