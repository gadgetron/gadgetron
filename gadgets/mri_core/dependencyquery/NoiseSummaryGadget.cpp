#include "NoiseSummaryGadget.h"

#include <boost/filesystem.hpp>
#include <limits>
#include <fstream>
#include "Dependency.h"

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
            auto message = new Gadgetron::GadgetContainerMessage<DependencyQuery::Dependency>();
            auto& dependencies = message->getObjectPtr()->dependencies;

            
            hoNDArray< std::complex<float> > noise_covariance_matrix;
            float noise_dwell_time;

            if ( !boost::filesystem::exists( p) ) {
                dependencies.append("status", "failed");
            } else {

                std::ifstream infile(p.string(), std::ios::in|std::ios::binary);

                if (infile.good()) {
                    //Read the XML header of the noise scan
                    uint32_t xml_length;
                    infile.read( reinterpret_cast<char*>(&xml_length), 4);
                    std::string xml_str(xml_length,'\0');
                    infile.read(const_cast<char*>(xml_str.c_str()), xml_length);
	
                    infile.read( reinterpret_cast<char*>(&noise_dwell_time), sizeof(float));
                    
                    size_t len;
                    infile.read( reinterpret_cast<char*>(&len), sizeof(size_t));

                    auto buf = std::make_unique<char[]>(len);

                    infile.read(buf.get(), len);

                    if ( !noise_covariance_matrix.deserialize(buf.get(), len) ) {
                        GERROR("Unable to deserialize matrix\n");
                        return GADGET_FAIL;
                    }
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

                dependencies.append("noise_dwell_time_us",noise_dwell_time);
                dependencies.append("min_sigma",min_sigma);
                dependencies.append("max_sigma",max_sigma);
                dependencies.append("mean_sigma",mean_sigma);
                dependencies.append("channels", static_cast<long>(coils));
                dependencies.append("status", "success");
            }

            this->next()->putq(message);
            
        }
        
        return GADGET_OK;
    }

    GADGET_FACTORY_DECLARE(NoiseSummaryGadget)

} // namespace Gadgetron
