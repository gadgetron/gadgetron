#include "AutoScaleGadget.h"

using namespace Gadgetron;
using namespace Gadgetron::Core;

namespace {
	template<class T>
	Image<T> autoscale(Image<T> &image, float max_value, unsigned int histogram_bins) {
		auto &header = std::get<ISMRMRD::ImageHeader>(image);
		auto &data = std::get<hoNDArray<T>>(image);
		auto &meta = std::get<optional<ISMRMRD::MetaContainer>>(image);	
		if (header.image_type == ISMRMRD::ISMRMRD_IMTYPE_MAGNITUDE) { //Only scale magnitude images for now
			auto max = *std::max_element(data.begin(), data.end());
			auto current_scale_ = 1.0;
			auto histogram_ = std::vector<size_t>(histogram_bins);    
			std::fill(histogram_.begin(), histogram_.end(), 0);
			for (auto& d : data) {
				size_t bin = static_cast<size_t>(std::floor((d/max)*histogram_bins));
				if (bin >= histogram_bins) {
					bin = histogram_bins-1;
				}
				histogram_[bin]++;
			}
			//Find 99th percentile
			long long cumsum = 0;
			size_t counter = 0;
			while (cumsum < (0.99*data.get_number_of_elements())) {
				cumsum += (long long)(histogram_[counter++]);
			}
			max = (counter+1)*(max/histogram_bins);
			GDEBUG("Max: %f\n",max);
			current_scale_ = max_value/max;
			for (auto& d : data){
				d *= current_scale_;
			}
		}
		return Image<T>(header,data,meta);
	}

    template<class T>
    Image<std::complex<T>> autoscale(Image<std::complex<T>> &image, float max_value, unsigned int histogram_bins) {
        GERROR("Autoscaling image is not well defined for complex images.");
		return image;
    }
}

namespace Gadgetron {
	AnyImage AutoScaleGadget::process_function(AnyImage image) const {
		return visit([&](auto &image) -> AnyImage { return autoscale(image, max_value, histogram_bins); }, image);
	}
	GADGETRON_GADGET_EXPORT(AutoScaleGadget);
}
