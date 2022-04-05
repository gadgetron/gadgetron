#include "PseudoReplicatorGadget.h"
#include <random>

namespace Gadgetron {

    void PseudoReplicatorGadget::process(Core::InputChannel<IsmrmrdReconData>& input, Core::OutputChannel& out) {

        for (IsmrmrdReconData reconData : input) {
			std::mt19937 engine(5489UL);
			std::normal_distribution<float> distribution;

			auto m_copy = reconData;

			//First just send the normal data to obtain standard image
			out.push(std::move(reconData));

			//Now for the noisy projections
			for (int i =0; i < repetitions; i++){

				auto cm = new GadgetContainerMessage<IsmrmrdReconData>();
				*cm->getObjectPtr() = m_copy;
				auto & datasets = cm->getObjectPtr()->rbit_;

				for (auto & buffer : datasets){
					auto & data = buffer.data_.data_;
					auto dataptr = data.get_data_ptr();
					for (size_t k =0; k <  data.get_number_of_elements(); k++){
						dataptr[k] += std::complex<float>(distribution(engine),distribution(engine));
					}
				}
				GDEBUG("Sending out Pseudoreplica\n");
				out.push(std::move(cm));
			}
        }
    }
    GADGETRON_GADGET_EXPORT(PseudoReplicatorGadget);

}
