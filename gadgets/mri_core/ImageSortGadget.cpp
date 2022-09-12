#include "ImageSortGadget.h"
#include "Node.h"
#include "Types.h"
#include "log.h"

using namespace Gadgetron::Core;

namespace Gadgetron {

    ImageSortGadget::ImageSortGadget(const Context &context, const GadgetProperties &properties)
        : ChannelGadget(context, properties) {}

    void ImageSortGadget::process(InputChannel<Core::AnyImage> &input, OutputChannel &output) {
        
        // Lambda, gets index from correct field in ImageHeader based on parameterized sorting dimension
        auto getImageIndex = [&](const auto& i){   
          auto &header = std::get<ISMRMRD::ImageHeader>(i);
          if (sorting_dimension.size() == 0) {
              return uint16_t(-1);
            } else if (sorting_dimension.compare("average") == 0) {
              return header.average;
            } else if (sorting_dimension.compare("slice") == 0) {
              return header.slice;
            } else if (sorting_dimension.compare("contrast") == 0) {
              return header.contrast;
            } else if (sorting_dimension.compare("phase") == 0) {
              return header.phase;
            } else if (sorting_dimension.compare("repetition") == 0) {
              return header.repetition;
            } else if (sorting_dimension.compare("set") == 0) {
              return header.set;
            } else {
              return uint16_t(-1);
            }
          return uint16_t(-1);
        };

        // Lambda, adds image and correct index to vector of ImageEntries 
        auto addToImages = [&](auto image){ 
          if (getImageIndex(image) < 0) {
            output.push(image);
          }
          images_.push_back(ImageEntry{image, getImageIndex(image)});
        };

        // Lambda, comparison method for ImageEntry indices
        auto image_entry_compare = [](const auto& i, const auto& j) {return i.index_<j.index_; };

        // Add all the images from input channel to vector of ImageEntries
        for (auto image : input) {
          visit(addToImages, image);
        }

        // If vector of ImageEntries isn't empty, sort them and write to output channel 
        if (images_.size()) {
          std::sort(images_.begin(),images_.end(), image_entry_compare);
          for (auto it = images_.begin(); it != images_.end(); it++) {
            visit([&](auto image){output.push(image);}, it->image_);
          }
          images_.clear();
        }
    }
    GADGETRON_GADGET_EXPORT(ImageSortGadget);
}