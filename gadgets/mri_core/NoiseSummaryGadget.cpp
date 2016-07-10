#include "GadgetIsmrmrdReadWrite.h"
#include "NoiseSummaryGadget.h"

#include <boost/filesystem.hpp>
#include <limits>
#include <fstream>

using namespace boost::filesystem;

namespace Gadgetron
{
    NoiseSummaryGadget::NoiseSummaryGadget()
    {
        processed_in_close_ = false;
    }

    NoiseSummaryGadget::~NoiseSummaryGadget()
    {
    }

    int NoiseSummaryGadget::close(unsigned long flags)
    {
        if ( BaseClass::close(flags) != GADGET_OK ) return GADGET_FAIL;

        if ( !processed_in_close_ ) {
            processed_in_close_ = true;
            
            if ( !workingDirectory.value().empty() ) {
                noise_dependency_folder_ = workingDirectory.value();
            } else {
                GERROR("Unable to determin noise dependecy folder\n");
                return GADGET_FAIL;
            }

            GDEBUG_STREAM("Noise dependency folder is " << noise_dependency_folder_);


            // list the content in the noise dependency folder
            path p = path(noise_dependency_folder_) / path(noise_file.value());

            // declear the attributes
            Gadgetron::GadgetContainerMessage<ISMRMRD::MetaContainer>* m1 = new Gadgetron::GadgetContainerMessage<ISMRMRD::MetaContainer>();
            
            hoNDArray< std::complex<float> > noise_covariance_matrix;
            float noise_dwell_time;

            if ( !boost::filesystem::exists( p) ) {
                m1->getObjectPtr()->append("status", "failed");                
            } else {

                std::ifstream infile;
                infile.open (p.string(), std::ios::in|std::ios::binary);

                if (infile.good()) {
                    //Read the XML header of the noise scan
                    uint32_t xml_length;
                    infile.read( reinterpret_cast<char*>(&xml_length), 4);
                    std::string xml_str(xml_length,'\0');
                    infile.read(const_cast<char*>(xml_str.c_str()), xml_length);
	
                    infile.read( reinterpret_cast<char*>(&noise_dwell_time), sizeof(float));
                    
                    size_t len;
                    infile.read( reinterpret_cast<char*>(&len), sizeof(size_t));
                    
                    char* buf = new char[len];
                    if ( buf == NULL ) return false;
                    
                    infile.read(buf, len);

                    if ( !noise_covariance_matrix.deserialize(buf, len) ) {
                        delete [] buf;
                        GERROR("Unable to deserialize matrix\n");
                        return GADGET_FAIL;
                    }

                    delete [] buf;
                    infile.close();
                } else {
                    GDEBUG("Noise covariance matrix file is not found. Error\n");
                    return GADGET_FAIL;
                }

                size_t coils = noise_covariance_matrix.get_size(0);
                
                //Collect stats
                float mean_sigma = 0.0;
                float max_sigma = 0.0;
                float min_sigma = std::numeric_limits<float>::max();
                
                for (size_t c = 0; c < coils; c++) {
                    float sigma = std::sqrt(std::real(noise_covariance_matrix[c*coils+c]));
                    mean_sigma += sigma;
                    if (sigma > max_sigma) max_sigma = sigma;
                    if (sigma < min_sigma) min_sigma = sigma;
                }

                GDEBUG("Min Sigma: %f\n", min_sigma);
                GDEBUG("Max Sigma: %f\n", max_sigma);
                GDEBUG("Mean Sigma: %f\n", mean_sigma);

                m1->getObjectPtr()->append("noise_dwell_time_us",noise_dwell_time);
                m1->getObjectPtr()->append("min_sigma",min_sigma);
                m1->getObjectPtr()->append("max_sigma",max_sigma);
                m1->getObjectPtr()->append("mean_sigma",mean_sigma);
                m1->getObjectPtr()->append("channels", static_cast<long>(coils));
                m1->getObjectPtr()->append("status", "success");
            }
            
            // send the found dependencies
            GadgetContainerMessage<GadgetMessageIdentifier>* mb = new GadgetContainerMessage<GadgetMessageIdentifier>();
            mb->getObjectPtr()->id = GADGET_MESSAGE_DEPENDENCY_QUERY;
            mb->cont(m1);

            int ret =  this->controller_->output_ready(mb);
        }
        
        return GADGET_OK;
    }

    GADGET_FACTORY_DECLARE(NoiseSummaryGadget)

} // namespace Gadgetron
