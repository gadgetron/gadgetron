#include "python_toolbox.h"

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
    std::cout << evens.get_number_of_elements() << std::endl;

    return 0;
}
