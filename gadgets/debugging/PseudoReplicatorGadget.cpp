#include "PseudoReplicatorGadget.h"

namespace Gadgetron {

    void PseudoReplicatorGadget::process(Core::InputChannel<IsmrmrdReconData>& input, Core::OutputChannel& out) {
        for (IsmrmrdReconData reconData : input) {
			std::mt19937 engine(seed);
			std::normal_distribution<float> distribution;

			auto reconDataCopy = reconData;

			// First just send the normal data to obtain standard image
			out.push(std::move(reconData));

			// Now for the noisy projections
			for (int i =0; i < repetitions; i++){
				auto reconDataTemp = reconDataCopy;
				auto & datasets = reconDataTemp.rbit_;
				for (auto & buffer : datasets){
					auto & data = buffer.data_.data_;
					auto dataptr = data.get_data_ptr();
					for (size_t k =0; k <  data.get_number_of_elements(); k++){
						dataptr[k] += std::complex<float>(distribution(engine),distribution(engine));
					}
				}
				out.push(std::move(reconDataTemp));
			}
        }
    }
    GADGETRON_GADGET_EXPORT(PseudoReplicatorGadget);

}