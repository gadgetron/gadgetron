#ifndef FATWATERGADGET_H
#define FATWATERGADGET_H

#include <fatwater.h>
#include "Gadget.h"
#include "hoNDArray.h"
#include "gadgetron_fatwater_export.h"

#include "mri_core_data.h"
#include "vector_td_io.h"
namespace Gadgetron{

  class EXPORTGADGETFATWATER FatWaterGadget : 
  public Gadget1<IsmrmrdImageArray>
    {
    public:
      GADGET_DECLARE(FatWaterGadget)
      FatWaterGadget();

      using value_range = vector_td<float,2>;

//      GADGET_PROPERTY(range_frequency_offset,value_range , "Range of field map values in Hz", value_range(-500,500));
      GADGET_PROPERTY(number_of_frequency_offsets,unsigned int, "Number of field map strengths", 200);
      GADGET_PROPERTY(range_r2star,value_range,"Range of R2* values in Hz",value_range(5,500));
      GADGET_PROPERTY(number_of_r2stars,unsigned int, "Number of R2* value to use during graph-cut",5);
      GADGET_PROPERTY(number_of_r2stars_fine,unsigned int,"Number of R2* values used for refinement after graph-cut",200);
      GADGET_PROPERTY(graph_cut_iterations,unsigned int, "Number of graph cut iterations to run",40);
      GADGET_PROPERTY(regularization_lambda,float,"Strength of the spatial regularization",0.02);
      GADGET_PROPERTY(regularization_offset,float, "Fixed value to add to the regularization for increased smoothness in low signal areas",0.01);
      GADGET_PROPERTY(do_gradient_descent, bool, "Use gradient descent after graph-cut",true);
      GADGET_PROPERTY(downsample_data, unsigned int, "Number of times to downsample data before calculating the field map",0);
      GADGET_PROPERTY(save_field_map , bool, "Save the field map",false);
      GADGET_PROPERTY(save_r2star_map, bool, "Save the R2* map",false);
      GADGET_PROPERTY(sample_time_us, float, "Sample time in microseconds for frequency offset correction. Set to 0 for disabled",0);


	
    protected:
      virtual int process(GadgetContainerMessage<IsmrmrdImageArray>* m1);
      virtual int process_config(ACE_Message_Block* mb);


      

    private:
      std::vector<float> echoTimes_;
      float fieldStrength_;
      FatWater::Config config;


      GadgetContainerMessage <IsmrmrdImageArray> *
      FatWaterImageArray(const FatWater::Parameters &parameters, hoNDArray <std::complex<float>> &&wfimages,
                               const hoNDArray <ISMRMRD::ImageHeader> &headers,
                               const std::vector<ISMRMRD::MetaContainer> &metadata) const;

      GadgetContainerMessage <ISMRMRD::ImageHeader> *
      MakeImageMessage(FatWater::Parameters parameters, hoNDArray<float> &&array,
                       const hoNDArray<ISMRMRD::ImageHeader> &header, uint16_t image_index) const;
  };
}
#endif //FATWATERGADGET_H
