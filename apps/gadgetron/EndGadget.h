/*
 * EndGadget.h
 *
 *  Created on: Nov 3, 2011
 *      Author: hansenms
 */

#ifndef ENDGADGET_H_
#define ENDGADGET_H_

#include "Gadget.h"
#include "GadgetMessageInterface.h"

namespace Gadgetron{
class EndGadget : public Gadget
{
	virtual int close(unsigned long flags)
	{
		GDEBUG("Close called in EndGadget with flags %d\n", flags);

		GadgetContainerMessage<GadgetMessageIdentifier>* mb =
				new GadgetContainerMessage<GadgetMessageIdentifier>();

		mb->getObjectPtr()->id = GADGET_MESSAGE_CLOSE;

		if (controller_->output_ready(mb) < 0) {
			return GADGET_FAIL;
		}

		GDEBUG("Calling close in base class  with flags %d\n", flags);
		return Gadget::close(flags);
	}

protected:
	virtual int process(ACE_Message_Block *m)
	{
		m->release();
		return 0;
	}

	virtual int next_step(ACE_Message_Block *m)
	{
		m->release();
		return 0;
	}

	virtual int process_config(ACE_Message_Block * m) {
		m->release();
		return 0;
	}

};
}

#endif /* ENDGADGET_H_ */
