#include "python_toolbox.h"
#include "ismrmrd/ismrmrd.h"

using namespace Gadgetron;

int main(int argc, char** argv)
{
    int a = -42;
    float b = 3.141592;
    std::string c("hello, world");
    unsigned int d(117);
    std::complex<double> e(2.12894, -1.103103);

    std::vector<size_t> dims;
    dims.push_back(4);
    dims.push_back(4);
    dims.push_back(4);
    hoNDArray<std::complex<float> > arr(dims);

    // Call a function with no return value (print all arguments)
    PythonFunction<> foo("__builtin__", "print");
    foo(a, b, c, d, e, arr);

    // Call a function with a single return value
    PythonFunction<float> atan2("math", "atan2");
    int x = 7, y = 4;
    float atan = atan2(x, y);
    std::cout << atan << std::endl;

    // Call a function that returns a tuple
    PythonFunction<float,float> divmod("__builtin__", "divmod");
    float w = 6.89;
    float z = 4.12;
    float fsum = 0, fdiff = 0;
    std::tie(fsum, fdiff) = divmod(w, z);
    std::cout << fsum << ", " << fdiff << std::endl;

    // Call a function that expects an iterable argument (tuple)
    PythonFunction<int> tuplen("__builtin__", "len");
    int l = tuplen(std::make_tuple(-7, 0, 7));
    std::cout << "tuple length: " << l << std::endl;

    // Generate an hoNDArray of even #s using numpy
    PythonFunction<hoNDArray<float>> arange("numpy", "arange");
    hoNDArray<float> evens = arange(0, 100, 2, "f64");
    std::cout << "number of even numbers between 0 and 100: " <<
            evens.get_number_of_elements() << std::endl;

    {
        GILLock gl;     // this is needed
        boost::python::object main(boost::python::import("__main__"));
        boost::python::object global(main.attr("__dict__"));
        boost::python::exec("def modify(head): head.version = 42; return head",
                global, global);
    }

    ISMRMRD::ImageHeader img_head, img_head2;
    img_head.version = 0;
    std::cout << "version before: " << img_head.version << std::endl;
    PythonFunction<ISMRMRD::ImageHeader> modify_img_header("__main__", "modify");
    img_head2 = modify_img_header(img_head);
    std::cout << "version after: " << img_head2.version << std::endl;

    ISMRMRD::AcquisitionHeader acq_head, acq_head2;
    acq_head.version = 0;
    std::cout << "version before: " << img_head.version << std::endl;
    PythonFunction<ISMRMRD::AcquisitionHeader> modify_acq_header("__main__", "modify");
    acq_head2 = modify_acq_header(acq_head);
    std::cout << "version after: " << acq_head2.version << std::endl;

    {
        GILLock gl;     // this is needed
        boost::python::object main(boost::python::import("__main__"));
        boost::python::object global(main.attr("__dict__"));
        boost::python::exec("from numpy.random import random\n"
                "def rand_cplx_array(length): \n"
                "    return random(length) + 1j * random(length)\n",
                global, global);
    }

    std::vector<std::complex<double> > vec;
    PythonFunction<std::vector<std::complex<double> > > make_vec("__main__", "rand_cplx_array");
    vec = make_vec(32);
    std::cout << vec[16] << std::endl;

    return 0;
}
