#include "Gadgetron.h"
#include "FlowPhaseSubtractionGadget.h"
#include "GadgetIsmrmrdReadWrite.h"

FlowPhaseSubtractionGadget::FlowPhaseSubtractionGadget()
{
}

FlowPhaseSubtractionGadget::~FlowPhaseSubtractionGadget()
{
}

int FlowPhaseSubtractionGadget::process_config(ACE_Message_Block* mb)
{
	boost::shared_ptr<ISMRMRD::ismrmrdHeader> cfg = parseIsmrmrdXMLHeader(std::string(mb->rd_ptr()));

	std::vector<long> dims;
	ISMRMRD::ismrmrdHeader::encoding_sequence e_seq = cfg->encoding();
	if (e_seq.size() != 1) {
		GADGET_DEBUG2("Number of encoding spaces: %d\n", e_seq.size());
		GADGET_DEBUG1("This Gadget only supports one encoding space\n");
		return GADGET_FAIL;
	}

	ISMRMRD::encodingSpaceType e_space = (*e_seq.begin()).encodedSpace();
	ISMRMRD::encodingSpaceType r_space = (*e_seq.begin()).reconSpace();
	ISMRMRD::encodingLimitsType e_limits = (*e_seq.begin()).encodingLimits();

	sets_ = e_limits.set().present() ? e_limits.set().get().maximum() + 1 : 1;

	if (sets_ > 2) {
		GADGET_DEBUG1("Phase subtraction only implemented for two sets for now\n");
		GADGET_DEBUG2("Number of sets detected: %d, bailing out.\n", sets_);
		return GADGET_FAIL;
	}

	return GADGET_OK;
}

int FlowPhaseSubtractionGadget::
process(GadgetContainerMessage<GadgetMessageImage>* m1,
		GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2)
{

	if (sets_ < 2) {
		if (this->next()->putq(m1) < 0) {
			m1->release();
			return GADGET_FAIL;
		}
		return GADGET_OK;
	}

	if (buffer_.message_count() < (sets_-1)) {
		buffer_.enqueue_tail(m1);
		return GADGET_OK;
	} else {
		ACE_Message_Block* mbq;
		if (buffer_.dequeue_head(mbq) < 0) {
			GADGET_DEBUG1("Message dequeue failed\n");
			return GADGET_FAIL;
		}

		GadgetContainerMessage<GadgetMessageImage>* pm1 = AsContainerMessage<GadgetMessageImage>(mbq);
		GadgetContainerMessage< hoNDArray< std::complex<float> > >* pm2 = AsContainerMessage<hoNDArray< std::complex<float> > >(mbq->cont());

		if (pm1->getObjectPtr()->set == m1->getObjectPtr()->set) {
			GADGET_DEBUG2("Data with same set number (%d,%d) detected. Error.\n", pm1->getObjectPtr()->set, m1->getObjectPtr()->set);
			return GADGET_FAIL;
		}

		if (m2->getObjectPtr()->get_number_of_elements() != pm2->getObjectPtr()->get_number_of_elements()) {
			GADGET_DEBUG1("Mismatch in number of elements detected. Bailing out.\n");
			pm1->release();
			m1->reset();
			return GADGET_FAIL;
		}

		std::complex<float>* p0 = pm2->getObjectPtr()->get_data_ptr();
		std::complex<float>* p1 = m2->getObjectPtr()->get_data_ptr();
		for (unsigned long i = 0; i < m2->getObjectPtr()->get_number_of_elements(); i++) {
			std::complex<float> tmp = std::polar((std::abs(p0[i])+std::abs(p1[i]))/2, std::arg(p1[i])-std::arg(p0[i]));
			p1[i] = tmp;
		}

		pm1->release();

		m1->getObjectPtr()->set = 0;

		if (this->next()->putq(m1) < 0) {
			m1->release();
			return GADGET_FAIL;
		}
		return GADGET_OK;
	}

	m1->release();
	return GADGET_OK;
}


GADGET_FACTORY_DECLARE(FlowPhaseSubtractionGadget)
