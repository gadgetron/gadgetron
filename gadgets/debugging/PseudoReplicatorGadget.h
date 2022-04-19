/*
 * PseudoReplicator.h

 *
 *  Created on: Jun 23, 2015
 *      Author: David Hansen
 */
#pragma once
#include "Gadget.h"
#include "mri_core_data.h"
#include "gadgetron_debugging_export.h"

namespace Gadgetron {

class EXPORTGADGETSDEBUGGING PseudoReplicatorGadget : public Gadget1<IsmrmrdReconData>{
public:
	GADGET_PROPERTY(repetitions,int,"Number of pseudoreplicas to produce",10);
	PseudoReplicatorGadget()  ;
	virtual ~PseudoReplicatorGadget();

	virtual int process_config(ACE_Message_Block *);
	virtual int process(GadgetContainerMessage<IsmrmrdReconData>*);

private:
	int repetitions_;
};

} /* namespace Gadgetron */
