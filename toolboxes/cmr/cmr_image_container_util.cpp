/** \file       cmr_image_container_util.cpp
    \brief      Functionalities to process image container
    \author     Hui Xue
*/

#include "cmr_image_container_util.h"
#include <fstream>

namespace Gadgetron
{
    template <typename T, unsigned int D>
    void resample_image_container(const hoNDImageContainer2D<hoMRImage<T, D> >& input, const std::vector<size_t>& dim_out, hoNDImageContainer2D<hoMRImage<T, D> >& output)
    {
        try
        {
            typedef hoMRImage<T, D> ImageType;

            std::vector<size_t> cols = input.cols();

            size_t SLC = cols.size();
            size_t PHS = cols[0];

            size_t slc;
            long phs;

            for (slc = 0; slc < SLC; slc++)
            {
#pragma omp parallel for default(none) shared(PHS, input, dim_out, output) private(phs)
                for (phs = 0; phs < PHS; phs++)
                {
                    ImageType& a_image = const_cast<ImageType&>(input(slc, phs));

                    if (Gadgetron::nrm2(a_image) > 0.1 && a_image.get_number_of_elements() > 0)
                    {
                        hoNDBoundaryHandlerBorderValue< hoMRImage<T, D> > bhBorderValue(a_image);

                        hoNDInterpolatorBSpline< hoMRImage<T, D>, D > interp(5);
                        interp.setArray(a_image);

                        hoMRImage<T, D> input_image(a_image);
                        hoMRImage<T, D> output_image;

                        Gadgetron::resampleImage(input_image, interp, dim_out, output_image);
                        output(slc, phs) = output_image;

                        output(slc, phs).header_ = input(slc, phs).header_;
                        output(slc, phs).attrib_ = input(slc, phs).attrib_;
                    }
                    else
                    {
                        output(slc, phs).create(dim_out);
                        Gadgetron::clear(output(slc, phs));
                        output(slc, phs).header_ = input(slc, phs).header_;
                        output(slc, phs).attrib_ = input(slc, phs).attrib_;
                    }

                    output(slc, phs).header_.matrix_size[0] = dim_out[0];
                    if (D > 1) output(slc, phs).header_.matrix_size[1] = dim_out[1];

                    if (D > 2)
                    {
                        output(slc, phs).header_.matrix_size[2] = dim_out[2];
                    }
                }
            }
        }
        catch (...)
        {
            GADGET_THROW("Error happened in resample_image_container(...) ... ");
        }
    }

    template EXPORTCMR void resample_image_container(const hoNDImageContainer2D<hoMRImage<float, 2> > & input, const std::vector<size_t>& dim_out, hoNDImageContainer2D<hoMRImage<float, 2> > & output);
    template EXPORTCMR void resample_image_container(const hoNDImageContainer2D<hoMRImage<double, 2> > & input, const std::vector<size_t>& dim_out, hoNDImageContainer2D<hoMRImage<double, 2> > & output);

    template EXPORTCMR void resample_image_container(const hoNDImageContainer2D<hoMRImage<float, 3> > & input, const std::vector<size_t>& dim_out, hoNDImageContainer2D<hoMRImage<float, 3> > & output);
    template EXPORTCMR void resample_image_container(const hoNDImageContainer2D<hoMRImage<double, 3> > & input, const std::vector<size_t>& dim_out, hoNDImageContainer2D<hoMRImage<double, 3> > & output);

    template EXPORTCMR void resample_image_container(const hoNDImageContainer2D<hoMRImage<std::complex<float>, 2> > & input, const std::vector<size_t>& dim_out, hoNDImageContainer2D<hoMRImage<std::complex<float>, 2> > & output);
    template EXPORTCMR void resample_image_container(const hoNDImageContainer2D<hoMRImage<std::complex<double>, 2> > & input, const std::vector<size_t>& dim_out, hoNDImageContainer2D<hoMRImage<std::complex<double>, 2> > & output);

    template EXPORTCMR void resample_image_container(const hoNDImageContainer2D<hoMRImage<std::complex<float>, 3> > & input, const std::vector<size_t>& dim_out, hoNDImageContainer2D<hoMRImage<std::complex<float>, 3> > & output);
    template EXPORTCMR void resample_image_container(const hoNDImageContainer2D<hoMRImage<std::complex<double>, 3> > & input, const std::vector<size_t>& dim_out, hoNDImageContainer2D<hoMRImage<std::complex<double>, 3> > & output);

    // -------------------------------------------------------------------------

    template <typename T, unsigned int D>
    void resample_image_container(const hoNDImageContainer2D<hoMRImage<T, D> >& input, double scale_ratio, hoNDImageContainer2D<hoMRImage<T, D> >& output)
    {
        try
        {
            typedef hoMRImage<T, D> ImageType;

            std::vector<size_t> cols = input.cols();

            size_t SLC = cols.size();
            size_t PHS = cols[0];

            size_t slc;
            long phs;

            for (slc = 0; slc < SLC; slc++)
            {
#pragma omp parallel for default(none) shared(PHS, input, output) private(phs)
                for (phs = 0; phs < PHS; phs++)
                {
                    ImageType& a_image = const_cast<ImageType&>(input(slc, phs));

                    std::vector<size_t> dim_out;
                    a_image.get_dimensions(dim_out);
                    for (size_t d = 0; d < dim_out.size(); d++) dim_out[d] *= scale_ratio;

                    if (Gadgetron::nrm2(a_image) > 0.1 && a_image.get_number_of_elements() > 0)
                    {
                        hoNDBoundaryHandlerBorderValue< hoMRImage<T, D> > bhBorderValue(a_image);

                        hoNDInterpolatorBSpline< hoMRImage<T, D>, D > interp(5);
                        interp.setArray(a_image);

                        hoMRImage<T, D> input_image(a_image);
                        hoMRImage<T, D> output_image;

                        Gadgetron::resampleImage(input_image, interp, dim_out, output_image);
                        output(slc, phs) = output_image;
                    }
                    else
                    {
                        output(slc, phs).create(dim_out);
                        Gadgetron::clear(output(slc, phs));
                    }

                    output(slc, phs).header_ = input(slc, phs).header_;
                    output(slc, phs).attrib_ = input(slc, phs).attrib_;

                    output(slc, phs).header_.matrix_size[0] = dim_out[0];
                    if (D > 1) output(slc, phs).header_.matrix_size[1] = dim_out[1];

                    if (D > 2)
                    {
                        output(slc, phs).header_.matrix_size[2] = dim_out[2];
                    }
                }
            }
        }
        catch (...)
        {
            GADGET_THROW("Error happened in resample_image_container(...) ... ");
        }
    }

    template EXPORTCMR void resample_image_container(const hoNDImageContainer2D<hoMRImage<float, 2> > & input, double scale_ratio, hoNDImageContainer2D<hoMRImage<float, 2> > & output);
    template EXPORTCMR void resample_image_container(const hoNDImageContainer2D<hoMRImage<double, 2> > & input, double scale_ratio, hoNDImageContainer2D<hoMRImage<double, 2> > & output);

    template EXPORTCMR void resample_image_container(const hoNDImageContainer2D<hoMRImage<float, 3> > & input, double scale_ratio, hoNDImageContainer2D<hoMRImage<float, 3> > & output);
    template EXPORTCMR void resample_image_container(const hoNDImageContainer2D<hoMRImage<double, 3> > & input, double scale_ratio, hoNDImageContainer2D<hoMRImage<double, 3> > & output);

    template EXPORTCMR void resample_image_container(const hoNDImageContainer2D<hoMRImage<std::complex<float>, 2> > & input, double scale_ratio, hoNDImageContainer2D<hoMRImage<std::complex<float>, 2> > & output);
    template EXPORTCMR void resample_image_container(const hoNDImageContainer2D<hoMRImage<std::complex<double>, 2> > & input, double scale_ratio, hoNDImageContainer2D<hoMRImage<std::complex<double>, 2> > & output);

    template EXPORTCMR void resample_image_container(const hoNDImageContainer2D<hoMRImage<std::complex<float>, 3> > & input, double scale_ratio, hoNDImageContainer2D<hoMRImage<std::complex<float>, 3> > & output);
    template EXPORTCMR void resample_image_container(const hoNDImageContainer2D<hoMRImage<std::complex<double>, 3> > & input, double scale_ratio, hoNDImageContainer2D<hoMRImage<std::complex<double>, 3> > & output);

    // -------------------------------------------------------------------------

    template <typename T>
    void convert_4D_array_to_container(const hoNDArray<T>& input, hoNDImageContainer2D < hoMRImage<T, 2> > & output)
    {
        try
        {
            size_t RO = input.get_size(0);
            size_t E1 = input.get_size(1);
            size_t PHS = input.get_size(2);
            size_t SLC = input.get_size(3);

            std::vector<size_t> cols(SLC, PHS);
            output.create(cols);

            std::vector<size_t> dims(2);
            dims[0] = RO;
            dims[1] = E1;

            size_t slc, phs;

            for (slc = 0; slc < SLC; slc++)
            {
                for (phs = 0; phs < PHS; phs++)
                {
                    output(slc, phs).create(dims);

                    memcpy(output(slc, phs).begin(), &input(0, 0, phs, slc), sizeof(T) * RO * E1);
                }
            }
        }
        catch (...)
        {
            GADGET_THROW("Error happened in convert_4D_array_to_container(...) ... ");
        }
    }

    template EXPORTCMR void convert_4D_array_to_container(const hoNDArray<float>& input, hoNDImageContainer2D < hoMRImage<float, 2> > & output);
    template EXPORTCMR void convert_4D_array_to_container(const hoNDArray<double>& input, hoNDImageContainer2D < hoMRImage<double, 2> > & output);
    template EXPORTCMR void convert_4D_array_to_container(const hoNDArray<std::complex<float> >& input, hoNDImageContainer2D < hoMRImage<std::complex<float>, 2> > & output);
    template EXPORTCMR void convert_4D_array_to_container(const hoNDArray<std::complex<double> >& input, hoNDImageContainer2D < hoMRImage<std::complex<double>, 2> > & output);

    // -------------------------------------------------------------------------

    template <typename T>
    void convert_container_to_4D_array(const hoNDImageContainer2D < hoMRImage<T, 2> > & input, hoNDArray<T>& output)
    {
        try
        {
            std::vector<size_t> cols = input.cols();

            size_t SLC = cols.size();
            size_t PHS = cols[0];

            size_t RO = input(0, 0).get_size(0);
            size_t E1 = input(0, 0).get_size(1);

            output.create(RO, E1, PHS, SLC);
            Gadgetron::clear(output);

            size_t phs, slc;

            for (slc = 0; slc < SLC; slc++)
            {
                for (phs = 0; phs < PHS; phs++)
                {
                    memcpy(&output(0, 0, phs, slc), input(slc, phs).begin(), sizeof(T) * RO * E1);
                }
            }
        }
        catch (...)
        {
            GADGET_THROW("Error happened in convert_container_to_4D_array(...) ... ");
        }
    }

    template EXPORTCMR void convert_container_to_4D_array(const hoNDImageContainer2D < hoMRImage<float, 2> > & input, hoNDArray<float>& output);
    template EXPORTCMR void convert_container_to_4D_array(const hoNDImageContainer2D < hoMRImage<double, 2> > & input, hoNDArray<double>& output);
    template EXPORTCMR void convert_container_to_4D_array(const hoNDImageContainer2D < hoMRImage<std::complex<float>, 2> > & input, hoNDArray<std::complex<float> >& output);
    template EXPORTCMR void convert_container_to_4D_array(const hoNDImageContainer2D < hoMRImage<std::complex<double>, 2> > & input, hoNDArray<std::complex<double> >& output);

    // -------------------------------------------------------------------------

    template <typename T>
    void sort_4D_array(const hoNDArray<T>& input, const std::vector<size_t> basal_to_apex_slices, hoNDArray<T>& output)
    {
        try
        {
            output = input;

            size_t RO = input.get_size(0);
            size_t E1 = input.get_size(1);
            size_t PHS = input.get_size(2);
            size_t SLC = input.get_size(3);

            size_t slc;
            for (slc=0; slc<SLC; slc++)
            {
                memcpy(&output(0, 0, 0, basal_to_apex_slices[slc]), &input(0, 0, 0, slc), sizeof(T) * RO * E1 * PHS);
            }
        }
        catch (...)
        {
            GADGET_THROW("Error happened in sort_4D_array(...) ... ");
        }
    }

    template EXPORTCMR void sort_4D_array(const hoNDArray<float>& input, const std::vector<size_t> basal_to_apex_slices, hoNDArray<float>& output);
    template EXPORTCMR void sort_4D_array(const hoNDArray<double>& input, const std::vector<size_t> basal_to_apex_slices, hoNDArray<double>& output);
    template EXPORTCMR void sort_4D_array(const hoNDArray<std::complex<float> >& input, const std::vector<size_t> basal_to_apex_slices, hoNDArray<std::complex<float> >& output);
    template EXPORTCMR void sort_4D_array(const hoNDArray<std::complex<double> >& input, const std::vector<size_t> basal_to_apex_slices, hoNDArray<std::complex<double> >& output);

    // -------------------------------------------------------------------------

    template <typename T>
    void sort_image_container(const hoNDImageContainer2D < hoMRImage<T, 2> > & input, const std::vector<size_t> basal_to_apex_slices, hoNDImageContainer2D < hoMRImage<T, 2> > & output)
    {
        try
        {
            output.copyFrom(input);

            std::vector<size_t> cols = input.cols();
            size_t SLC = cols.size();

            size_t phs, slc;
            for (slc = 0; slc < SLC; slc++)
            {
                for (phs = 0; phs < cols[slc]; phs++)
                {
                    output(basal_to_apex_slices[slc], phs) = input(slc, phs);
                }
            }
        }
        catch (...)
        {
            GADGET_THROW("Error happened in sort_image_container(...) ... ");
        }
    }

    template EXPORTCMR void sort_image_container(const hoNDImageContainer2D < hoMRImage<float, 2> > & input, const std::vector<size_t> basal_to_apex_slices, hoNDImageContainer2D < hoMRImage<float, 2> > & output);
    template EXPORTCMR void sort_image_container(const hoNDImageContainer2D < hoMRImage<double, 2> > & input, const std::vector<size_t> basal_to_apex_slices, hoNDImageContainer2D < hoMRImage<double, 2> > & output);
    template EXPORTCMR void sort_image_container(const hoNDImageContainer2D < hoMRImage<std::complex<float>, 2> > & input, const std::vector<size_t> basal_to_apex_slices, hoNDImageContainer2D < hoMRImage<std::complex<float>, 2> > & output);
    template EXPORTCMR void sort_image_container(const hoNDImageContainer2D < hoMRImage<std::complex<double>, 2> > & input, const std::vector<size_t> basal_to_apex_slices, hoNDImageContainer2D < hoMRImage<std::complex<double>, 2> > & output);

    // -------------------------------------------------------------------------

    template <typename ImageType>
    void serialize_contrainer_meta_attributes(hoNDImageContainer2D<ImageType>& input, std::vector<std::vector<std::string> >& attribs)
    {
        try
        {
            std::vector<size_t> cols = input.cols();

            size_t R = cols.size();
            attribs.resize(R);

            size_t r, c;

            for (r = 0; r < R; r++)
            {
                for (c = 0; c < cols[r]; c++)
                {
                    std::stringstream xml_ss;
                    ISMRMRD::serialize(input(r, c).attrib_, xml_ss);
                    std::string xml_str = xml_ss.str();

                    attribs[r].push_back(xml_str);
                }
            }
        }
        catch (...)
        {
            GERROR_STREAM("Exceptions happened in serialize_contrainer_meta_attributes(...) ...");
        }
    }

    template EXPORTCMR void serialize_contrainer_meta_attributes(hoNDImageContainer2D<hoMRImage<float, 2> >& input, std::vector<std::vector<std::string> >& attribs);
    template EXPORTCMR void serialize_contrainer_meta_attributes(hoNDImageContainer2D<hoMRImage<double, 2> > & input, std::vector<std::vector<std::string> >& attribs);
    template EXPORTCMR void serialize_contrainer_meta_attributes(hoNDImageContainer2D<hoMRImage<short, 2> > & input, std::vector<std::vector<std::string> >& attribs);
    template EXPORTCMR void serialize_contrainer_meta_attributes(hoNDImageContainer2D<hoMRImage<unsigned short, 2> > & input, std::vector<std::vector<std::string> >& attribs);
    template EXPORTCMR void serialize_contrainer_meta_attributes(hoNDImageContainer2D<hoMRImage<std::complex<float>, 2> > & input, std::vector<std::vector<std::string> >& attribs);
    template EXPORTCMR void serialize_contrainer_meta_attributes(hoNDImageContainer2D<hoMRImage<std::complex<double>, 2> > & input, std::vector<std::vector<std::string> >& attribs);

    template EXPORTCMR void serialize_contrainer_meta_attributes(hoNDImageContainer2D<hoMRImage<float, 3> > & input, std::vector<std::vector<std::string> >& attribs);
    template EXPORTCMR void serialize_contrainer_meta_attributes(hoNDImageContainer2D<hoMRImage<double, 3> > & input, std::vector<std::vector<std::string> >& attribs);
    template EXPORTCMR void serialize_contrainer_meta_attributes(hoNDImageContainer2D<hoMRImage<short, 3> > & input, std::vector<std::vector<std::string> >& attribs);
    template EXPORTCMR void serialize_contrainer_meta_attributes(hoNDImageContainer2D<hoMRImage<unsigned short, 3> > & input, std::vector<std::vector<std::string> >& attribs);
    template EXPORTCMR void serialize_contrainer_meta_attributes(hoNDImageContainer2D<hoMRImage<std::complex<float>, 3> > & input, std::vector<std::vector<std::string> >& attribs);
    template EXPORTCMR void serialize_contrainer_meta_attributes(hoNDImageContainer2D<hoMRImage<std::complex<double>, 3> > & input, std::vector<std::vector<std::string> >& attribs);

    // -------------------------------------------------------------------------

    void write_contrainer_meta_attributes(const std::vector<std::vector<std::string> >& attribs, const std::string& prefix)
    {
        try
        {
            size_t R = attribs.size();

            size_t r, c;

            for (r = 0; r < R; r++)
            {
                for (c = 0; c < attribs[r].size(); c++)
                {
                    std::stringstream filename;
                    filename << prefix << "_row_" << r << "_col_" << c << ".attrib";

                    std::string xml_str = attribs[r][c];
                    uint32_t xml_length = static_cast<uint32_t>(xml_str.size());

                    std::ofstream outfile;
                    outfile.open(filename.str().c_str(), std::ios::out | std::ios::binary);

                    if (outfile.good())
                    {
                        outfile.write(xml_str.c_str(), xml_length);
                        outfile.close();
                    }
                }
            }
        }
        catch (...)
        {
            GERROR_STREAM("Exceptions happened in write_contrainer_meta_attributes(...) ...");
        }
    }

    // -------------------------------------------------------------------------
}
