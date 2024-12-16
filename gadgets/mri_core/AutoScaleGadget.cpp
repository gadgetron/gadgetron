#include "AutoScaleGadget.h"

using namespace Gadgetron;

namespace {
	template<class T>
	mrd::Image<T> autoscale(mrd::Image<T> &image, float max_value, unsigned int histogram_bins) {
		// Only scale magnitude images for now
		if (image.head.image_type == mrd::ImageType::kMagnitude) {
			auto& data = image.data;
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
			// Find 99th percentile
			long long cumsum = 0;
			size_t counter = 0;
			while (cumsum < (0.99*data.size())) {
				cumsum += (long long)(histogram_[counter++]);
			}
			max = (counter+1)*(max/histogram_bins);
			current_scale_ = max_value/max;
			for (auto& d : data){
				d *= current_scale_;
			}
		}
		return image;
	}

    template<class T>
    mrd::Image<std::complex<T>> autoscale(mrd::Image<std::complex<T>> &image, float max_value, unsigned int histogram_bins) {
        GERROR("Autoscaling image is not well defined for complex images.");
		return image;
    }
}

namespace Gadgetron {
	mrd::AnyImage AutoScaleGadget::process_function(mrd::AnyImage image) const {
		return visit([&](auto &image) -> mrd::AnyImage { return autoscale(image, max_value, histogram_bins); }, image);
	}
	GADGETRON_GADGET_EXPORT(AutoScaleGadget);
}
