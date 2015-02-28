#include "log.h"
#include "python_toolbox.h"
#include "hoNDArray_fileio.h"

using namespace Gadgetron;

int main(int argc, char** argv)
{
  GINFO("This is the Python test application\n");

  if (argc < 2) {
    GERROR("You must supply an input file\n");
    return -1;
  }

  boost::shared_ptr< hoNDArray< std::complex<float> > > source_data = read_nd_array< std::complex<float> >(argv[1]);

  size_t coils = source_data->get_size(2);
  size_t ny = source_data->get_size(1);
  size_t nx = source_data->get_size(0);
  GINFO("Array dimensions [%d, %d, %d]\n", nx, ny, coils);

  hoNDArray< std::complex<float> > unmix;
  hoNDArray<float> gmap;

  PythonFunction<hoNDArray<std::complex<float> >, hoNDArray<float> > calculate_grappa_unmixing("ismrmrdtools.grappa", "calculate_grappa_unmixing");

  std::tie(unmix, gmap) = calculate_grappa_unmixing(*source_data.get(), 3);

  write_nd_array(&unmix, "unmix.cplx");

  return 0;
}
