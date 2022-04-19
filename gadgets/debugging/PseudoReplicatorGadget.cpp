/*
 * PseudoReplicator.cpp
 *
 *  Created on: Jun 23, 2015
 *      Author: u051747
 */

#include "PseudoReplicatorGadget.h"
#include <random>
namespace Gadgetron {

PseudoReplicatorGadget::PseudoReplicatorGadget() : Gadget1<IsmrmrdReconData>() {
	// TODO Auto-generated constructor stub

}

PseudoReplicatorGadget::~PseudoReplicatorGadget() {
	// TODO Auto-generated destructor stub
}

int PseudoReplicatorGadget::process_config(ACE_Message_Block*) {

	repetitions_ = repetitions.value();
	return GADGET_OK;
}

int PseudoReplicatorGadget::process(GadgetContainerMessage<IsmrmrdReconData>* m) {

	std::mt19937 engine(5489UL);
	std::normal_distribution<float> distribution;

	auto m_copy = *m->getObjectPtr();
	//First just send the normal data to obtain standard image
	if (this->next()->putq(m) == GADGET_FAIL)
			return GADGET_FAIL;

	//Now for the noisy projections
	for (int i =0; i < repetitions_; i++){

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

		if (this->next()->putq(cm) == GADGET_FAIL)
			return GADGET_FAIL;

	}
	return GADGET_OK;

}

GADGET_FACTORY_DECLARE(PseudoReplicatorGadget)

} /* namespace Gadgetron */
