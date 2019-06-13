
#ifndef MFIOperator_H
#define MFIOperator_H

#include "hoNDArray.h"
#include "cuNDArray.h"

#ifdef USE_OMP
    #include "omp.h"
#endif // USE_OMP

#include "cuNDArray_utils.h"
#include "hoNDArray_utils.h"
#include "cuNFFT.h"

namespace Gadgetron{

  class MFIOperator{

  public:

      //Constructor
      MFIOperator();

      //Deconstructor
      ~MFIOperator();

      //Alternate constructor calls prepare
      MFIOperator(boost::shared_ptr<cuNFFT_plan<float,2>> plan_, cuNDArray<floatd2>& gpu_traj, std::vector<size_t> &data_d,std::vector<size_t> &image_d, uint64d2 &image_d_os, double &st,float dTE);


<<<<<<< HEAD
      bool prepare(boost::shared_ptr<cuNFFT_plan<float, 2>> plan_, cuNDArray<floatd2>& gpu_traj,
          std::vector<size_t>& data_d, std::vector<size_t>& image_d, uint64d2& image_d_os, double& st, float dTE);

      hoNDArray<complext<float>> MFI_apply(hoNDArray<complext<float>>& ho_image, hoNDArray<float> B0_map);
=======
      bool Prepare(boost::shared_ptr<cuNFFT_plan<float,2>> plan_, cuNDArray<floatd2>& gpu_traj, std::vector<size_t> &data_d,std::vector<size_t> &image_d, uint64d2 &image_d_os, double &st,float dTE);

      hoNDArray<complext<float>> MFI_Apply(hoNDArray<complext<float>> &ho_image,hoNDArray<float> B0_map);
>>>>>>> 78428447f8c051b458fbf61261aaddb589182734

  private:

      bool is_prepared;

      float sample_time;
      float delTE;

      int fmax;
      size_t L;

      std::vector<size_t> data_dims;
      std::vector<size_t> image_dims;

      hoNDArray<complext<float>> MFI_C;

      cuNDArray<complext<float>> cu_phase_mask;
      cuNDArray<complext<float>> cu_kspace_filter;

      boost::shared_ptr<cuNFFT_plan<float,2>> nfft_plan_;

<<<<<<< HEAD
      void calc_MFI_coeff();

      void calc_phase_mask();

      void calc_kspace_filter(std::vector<size_t>& matrix_size);
=======
      void Calc_MFI_Coeff();

      void Calc_Phase_Mask();

      void Calc_Kspace_Filter(std::vector<size_t> &matrix_size);
>>>>>>> 78428447f8c051b458fbf61261aaddb589182734

  };
}

#endif
