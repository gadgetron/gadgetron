
/** \file   mri_core_python.cpp
    \brief  Python binding for mri_core functinalities
    \author Hui Xue
*/

#include "mri_core_grappa_python.h"
#include "python_toolbox.h"
#include <boost/python.hpp>

#ifdef USE_OMP
    #include "omp.h"
#endif // USE_OMP

using namespace boost::python;

const std::string hello() {
    return std::string("hello, gadgetron");
}

BOOST_PYTHON_MODULE(gadgetron_toolbox_mri_core_python)
{
    // for test purpose
    def("hello", hello);

    class_<Gadgetron::grappa2D>("grappa2D")
        .def("initialize", &Gadgetron::grappa2D::initialize)
        .def("calib", &Gadgetron::grappa2D::calib)
        .def("recon", &Gadgetron::grappa2D::recon)
        .def("recon_data_matrix", &Gadgetron::grappa2D::recon_data_matrix)
        .def("recon_fill_back_kspace", &Gadgetron::grappa2D::recon_fill_back_kspace)
        .def("get_A", &Gadgetron::grappa2D::get_A)
        .def("get_B", &Gadgetron::grappa2D::get_B)
        .def("get_ker", &Gadgetron::grappa2D::get_ker)
        .def("get_data_A", &Gadgetron::grappa2D::get_data_A)
        .def("get_data_A_index", &Gadgetron::grappa2D::get_data_A_index)
        .def("help", &Gadgetron::grappa2D::help)
        .def("status", &Gadgetron::grappa2D::status)
        ;
}
