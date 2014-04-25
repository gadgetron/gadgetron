#ifndef IMAGEWRITER_H
#define IMAGEWRITER_H

#include <fstream>

#include "GadgetImageMessageReader.h"
#include "GadgetronCommon.h"
#include "gtPlusIOAnalyze.h"
#include "gtPlusIONifti.h"
#include "GtPlusDefinition.h"

namespace Gadgetron
{
    template <typename T> class ImageWriter : public GadgetImageMessageReader<T>
    {
    public:
        ImageWriter()
            : number_of_calls_(0)
        {}

        virtual ~ImageWriter() {};

        virtual ACE_Message_Block* read(ACE_SOCK_Stream* socket) 
        {
            // Invoke parent's read
            ACE_Message_Block* mb = GadgetImageMessageReader<T>::read(socket);

            if (!mb) {
                GADGET_DEBUG1("Read failed in parent\n");
                return 0;
            }

            GadgetContainerMessage<ISMRMRD::ImageHeader> * img_head_mb =
                dynamic_cast<GadgetContainerMessage<ISMRMRD::ImageHeader> *>(mb);

            if (!img_head_mb) {
                GADGET_DEBUG1("Failed in dynamic cast\n");
                mb->release();
                return 0;
            }

            //GADGET_DEBUG2("Received image with %d channels\n", img_head_mb->getObjectPtr()->channels);

            GadgetContainerMessage<hoNDArray< T > > * img_data_mb =
                dynamic_cast<GadgetContainerMessage<hoNDArray< T > > *>(img_head_mb->cont());

            if (!img_data_mb) {
                GADGET_DEBUG1("Failed in dynamic cast\n");
                mb->release();
                return 0;
            }

            if (this->process_image(img_head_mb->getObjectPtr(), img_data_mb->getObjectPtr()) < 0) {
                GADGET_DEBUG1("Failed to process image\n");
                mb->release();
                return 0;
            }

            return mb;
        }

        virtual int process_image(ISMRMRD::ImageHeader* img_head,
            hoNDArray< T >* data)
        {
            ACE_DEBUG( (LM_DEBUG, ACE_TEXT("Image Writer writing image\n")) );

            char filename[1024];

            switch (sizeof(T)) {

            case (8): //Complex float
                sprintf(filename, "out_%05d.cplx", (int)number_of_calls_);
                break;
            case (4): //Real floats
                sprintf(filename, "out_%05d.real", (int)number_of_calls_);
                break;
            case (2): //Unsigned short
                sprintf(filename, "out_%05d.short", (int)number_of_calls_);
                break;
            default:
                sprintf(filename, "out_%05d.cplx", (int)number_of_calls_);
                break;
            }

            std::ofstream outfile;
            outfile.open (filename, std::ios::out|std::ios::binary);

            if (outfile.good()) {
                int ndim = 4;
                int dims[4];
                size_t elements = 1;
                dims[0] = img_head->matrix_size[0]; elements*=dims[0];
                dims[1] = img_head->matrix_size[1]; elements*=dims[1];
                dims[2] = img_head->matrix_size[2]; elements*=dims[2];
                dims[3] = img_head->channels; elements*=dims[3];

                outfile.write((char*)&ndim,sizeof(int));
                outfile.write((char*)dims,sizeof(int)*4);
                outfile.write((char*)data->get_data_ptr(),sizeof(T)*elements);
                outfile.close();
                number_of_calls_++;
            } else {
                GADGET_DEBUG1("File is not good for writing\n");
                return GADGET_FAIL;
            }

            return GADGET_OK;
        }

    protected:
        size_t number_of_calls_;
    };

    // for images with meta attributes
    template <typename T> class ImageAttribWriter : public GadgetImageAttribMessageReader<T>
    {
    public:

        typedef GadgetImageAttribMessageReader<T> BaseClass;
        typedef typename BaseClass::size_t_type size_t_type;

        ImageAttribWriter() : number_of_calls_(0) {}
        virtual ~ImageAttribWriter() {};

        virtual ACE_Message_Block* read(ACE_SOCK_Stream* socket) 
        {
            // Invoke parent's read
            ACE_Message_Block* mb = BaseClass::read(socket);

            if (!mb)
            {
                GADGET_DEBUG1("Read failed in parent\n");
                return 0;
            }

            GadgetContainerMessage<ISMRMRD::ImageHeader> * img_head_mb = dynamic_cast<GadgetContainerMessage<ISMRMRD::ImageHeader> *>(mb);

            if (!img_head_mb)
            {
                GADGET_DEBUG1("Failed in dynamic cast\n");
                mb->release();
                return 0;
            }

            GadgetContainerMessage<hoNDArray< T > > * img_data_mb = dynamic_cast<GadgetContainerMessage<hoNDArray< T > > *>(img_head_mb->cont());

            if (!img_data_mb)
            {
                GADGET_DEBUG1("Failed in dynamic cast\n");
                mb->release();
                return 0;
            }

            GadgetContainerMessage<GtImageAttribType> * img_attrib_mb = dynamic_cast<GadgetContainerMessage<GtImageAttribType> *>(img_data_mb->cont());

            if (!img_head_mb)
            {
                GADGET_DEBUG1("Failed in dynamic cast\n");
                mb->release();
                return 0;
            }

            if (this->process_image(img_head_mb->getObjectPtr(), img_data_mb->getObjectPtr(), img_attrib_mb->getObjectPtr()) < 0)
            {
                GADGET_DEBUG1("Failed to process image\n");
                mb->release();
                return 0;
            }

            return mb;
        }

        virtual int process_image(ISMRMRD::ImageHeader* img_head, hoNDArray< T >* data, GtImageAttribType* img_attrib)
        {
            ACE_DEBUG( (LM_DEBUG, ACE_TEXT("ImageAttrib Writer writing image\n")) );

            char filename[1024], filenameMeta[1024];

            switch (sizeof(T))
            {
            case (8): //Complex float
                sprintf(filename, "out_%05d.cplx", (int)number_of_calls_);
                break;
            case (4): //Real floats
                sprintf(filename, "out_%05d.real", (int)number_of_calls_);
                break;
            case (2): //Unsigned short
                sprintf(filename, "out_%05d.short", (int)number_of_calls_);
                break;
            default:
                sprintf(filename, "out_%05d.cplx", (int)number_of_calls_);
                break;
            }

            sprintf(filenameMeta, "out_%05d.xml", (int)number_of_calls_);

            std::ofstream outfile;
            outfile.open (filename, std::ios::out|std::ios::binary);

            if (outfile.good())
            {
                int ndim = 4;
                int dims[4];
                size_t elements = 1;
                dims[0] = img_head->matrix_size[0]; elements*=dims[0];
                dims[1] = img_head->matrix_size[1]; elements*=dims[1];
                dims[2] = img_head->matrix_size[2]; elements*=dims[2];
                dims[3] = img_head->channels; elements*=dims[3];

                outfile.write((char*)&ndim,sizeof(int));
                outfile.write((char*)dims,sizeof(int)*4);
                outfile.write((char*)data->get_data_ptr(),sizeof(T)*elements);
                outfile.close();
                number_of_calls_++;
            }
            else
            {
                GADGET_DEBUG1("File is not good for writing\n");
                return GADGET_FAIL;
            }

            std::ofstream outfile_meta;
            outfile_meta.open (filename, std::ios::out|std::ios::binary);

            if ( outfile_meta.good() )
            {
                char* buf = NULL;
                size_t_type len;
                if ( !img_attrib->serialize(buf, len) )
                {
                    GADGET_DEBUG1("Image meta attributes serialization failed\n");
                    return GADGET_FAIL;
                }

                outfile_meta.write(buf, len);
                outfile_meta.close();
                delete [] buf;
            }
            else
            {
                GADGET_DEBUG1("File of xml header is not good for writing\n");
                return GADGET_FAIL;
            }

            return GADGET_OK;
        }

    protected:
        size_t number_of_calls_;
    };

    template <typename T> class AnalyzeImageAttribWriter : public ImageAttribWriter<T>
    {
    public:

        typedef ImageAttribWriter<T> BaseClass;
        typedef typename BaseClass::size_t_type size_t_type;

        AnalyzeImageAttribWriter() : BaseClass()
        {
            dst_folder_ = "./";
        }

        AnalyzeImageAttribWriter(const std::string& prefix, const std::string& dst_folder="./") : BaseClass()
        {
            prefix_ = prefix;
            dst_folder_ = dst_folder;
        }

        virtual ~AnalyzeImageAttribWriter() {}

        Gadgetron::gtPlus::gtPlusIOAnalyze gtplus_exporter_;

        virtual void write_image(hoNDArray< T >* data, const std::string& filename)
        {
            GADGET_EXPORT_ARRAY(dst_folder_, gtplus_exporter_, *data, filename);
        }

        virtual int process_image(ISMRMRD::ImageHeader* img_head, hoNDArray< T >* data, GtImageAttribType* img_attrib)
        {
            // image data role
            std::vector<std::string> dataRole;
            if ( !img_attrib->attribute4_.get(GTPLUS_DATA_ROLE, dataRole) )
            {
                dataRole.push_back("Image");
            }

            long long imageNumber;
            img_attrib->attribute1_.get(GTPLUS_IMAGENUMBER, 0, imageNumber);

            long long cha, slc, e2, con, phs, rep, set;
            img_attrib->attribute1_.get(GTPLUS_CHA,        0, cha);
            img_attrib->attribute1_.get(GTPLUS_SLC,        0, slc);
            img_attrib->attribute1_.get(GTPLUS_E2,         0, e2);
            img_attrib->attribute1_.get(GTPLUS_CONTRAST,   0, con);
            img_attrib->attribute1_.get(GTPLUS_PHASE,      0, phs);
            img_attrib->attribute1_.get(GTPLUS_REP,        0, rep);
            img_attrib->attribute1_.get(GTPLUS_SET,        0, set);

            std::ostringstream ostr;

            if ( !prefix_.empty() )
            {
                ostr << prefix_ << "_";
            }

            size_t n;
            for ( n=0; n<dataRole.size(); n++ )
            {
                ostr << dataRole[n] << "_";
            }

            ostr << "SLC" << slc << "_"
                 << "E2" << e2 << "_"
                 << "CON" << con << "_"
                 << "PHS" << phs << "_"
                 << "REP" << rep << "_"
                 << "SET" << set;

            std::string filename = ostr.str();

            ACE_DEBUG( (LM_DEBUG, ACE_TEXT("Writing image %s\n"), filename.c_str()) );

            this->write_image(data, filename);

            std::string filename_attrib(dst_folder_);
            filename_attrib.append("/");
            filename_attrib.append(filename);
            filename_attrib.append("_attrib.xml");

            std::ofstream outfile_attrib;
            outfile_attrib.open (filename_attrib.c_str(), std::ios::out|std::ios::binary);

            if ( outfile_attrib.good() )
            {
                char* buf = NULL;
                size_t_type len;
                if ( !img_attrib->serialize(buf, len) )
                {
                    GADGET_DEBUG1("Image meta attributes serialization failed\n");
                    return GADGET_FAIL;
                }

                outfile_attrib.write(buf, len);
                outfile_attrib.close();
                delete [] buf;
            }
            else
            {
                GADGET_DEBUG1("File of xml header is not good for writing\n");
                return GADGET_FAIL;
            }

            number_of_calls_++;

            return GADGET_OK;
        }

    protected:

        using BaseClass::number_of_calls_;

        std::string prefix_;
        std::string dst_folder_;
    };

    template <typename T> class AnalyzeComplexImageAttribWriter : public AnalyzeImageAttribWriter<T>
    {
    public:

        typedef AnalyzeImageAttribWriter<T> BaseClass;
        typedef typename BaseClass::size_t_type size_t_type;

        AnalyzeComplexImageAttribWriter() : BaseClass()
        {
        }

        AnalyzeComplexImageAttribWriter(const std::string& prefix, const std::string& dst_folder="./") : BaseClass(prefix, dst_folder)
        {
        }

        virtual ~AnalyzeComplexImageAttribWriter() {}

        using BaseClass::gtplus_exporter_;

        virtual void write_image(hoNDArray< T >* data, const std::string& filename)
        {
            GADGET_EXPORT_ARRAY_COMPLEX_REAL_IMAG(dst_folder_, gtplus_exporter_, *data, filename);
        }

    protected:

        using BaseClass::number_of_calls_;
        using BaseClass::dst_folder_;
    };

    template <typename T> class NiftiImageAttribWriter : public AnalyzeImageAttribWriter<T>
    {
    public:

        typedef AnalyzeImageAttribWriter<T> BaseClass;
        typedef typename BaseClass::size_t_type size_t_type;

        NiftiImageAttribWriter() : BaseClass()
        {
        }

        NiftiImageAttribWriter(const std::string& prefix, const std::string& dst_folder="./") : BaseClass(prefix, dst_folder)
        {
        }

        virtual ~NiftiImageAttribWriter() {}

        Gadgetron::gtPlus::gtPlusIONifti gtplus_nifiti_exporter_;

        virtual void write_image(hoNDArray< T >* data, const std::string& filename)
        {
            GADGET_EXPORT_ARRAY(dst_folder_, gtplus_nifiti_exporter_, *data, filename);
        }

    protected:

        using BaseClass::number_of_calls_;
        using BaseClass::dst_folder_;
    };

    template <typename T> class NiftiComplexImageAttribWriter : public AnalyzeImageAttribWriter<T>
    {
    public:

        typedef AnalyzeImageAttribWriter<T> BaseClass;
        typedef typename BaseClass::size_t_type size_t_type;

        NiftiComplexImageAttribWriter() : BaseClass()
        {
        }

        NiftiComplexImageAttribWriter(const std::string& prefix, const std::string& dst_folder="./") : BaseClass(prefix, dst_folder)
        {
        }

        virtual ~NiftiComplexImageAttribWriter() {}

        Gadgetron::gtPlus::gtPlusIONifti gtplus_nifiti_exporter_;

        virtual void write_image(hoNDArray< T >* data, const std::string& filename)
        {
            GADGET_EXPORT_ARRAY_COMPLEX(dst_folder_, gtplus_nifiti_exporter_, *data, filename);
        }

    protected:

        using BaseClass::number_of_calls_;
        using BaseClass::dst_folder_;
    };
}
#endif //IMAGE_WRITER
