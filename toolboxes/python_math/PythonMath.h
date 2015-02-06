#ifndef GADGETRON_PYTHON_MATH_H
#define GADGETRON_PYTHON_MATH_H

#include "python_math_export.h"
#include <boost/thread/mutex.hpp>
#include "hoNDArray.h"

namespace Gadgetron
{

  class EXPORTPYTHONMATH PythonMath
  {

  public:
    static PythonMath* instance();

    /**
       Wrapper for ismrmrdtools.grappa.calculate_grappa_unmixing(source_data, acc_factor, kernel_size=(4,5), data_mask=None, csm=None, regularization_factor=0.001, target_data=None):
     */
    void calculate_grappa_unmixing(hoNDArray< std::complex<float> >* source_data, unsigned int acc_factor, hoNDArray< std::complex<float> >* unmix_out,
				  unsigned int* kernel_size = 0, hoNDArray< unsigned int >* data_mask = 0, hoNDArray< std::complex<float> >* csm = 0,
				  float regularization_factor = 0.001, hoNDArray< std::complex<float> >* target_data = 0, 
				  hoNDArray< float >* gmap_out = 0);

  protected:
    ///Protected constructor. 
    PythonMath();

    static PythonMath* instance_;
    boost::mutex mtx_;
  };


}

#endif
