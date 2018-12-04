#include "Message.h"

#include "GadgetContainerMessage.h"

namespace Gadgetron::Core {


    template<class T> GadgetContainerMessageBase* TypedMessage<T>::to_container_message() {
                return new GadgetContainerMessage<T>(this->take_data());
    }
}