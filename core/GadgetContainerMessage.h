#pragma once

#include <string>
#include "Message.h"
#include "LegacyACE.h"
#include <typeinfo>
#include "log.h"

namespace Gadgetron {


    class GadgetContainerMessageBase : public ACE_Message_Block {
    public:
        virtual std::unique_ptr<Core::MessageChunk> take_message() = 0;
        virtual ~GadgetContainerMessageBase() = default;
    };



    template<class T>
    class GadgetContainerMessage : public GadgetContainerMessageBase {

    public:

        template<class... ARGS>
        explicit GadgetContainerMessage(ARGS&&... xs){
            message = std::make_unique<Core::TypedMessageChunk<T>>(std::forward<ARGS>(xs)...);
            data = &message->data;
        }

         ~GadgetContainerMessage() override = default;

        std::unique_ptr<Core::MessageChunk> take_message() override {
            data = nullptr;
            return std::move(message);
        }

        T *getObjectPtr() {
            return data;
        }

        GadgetContainerMessage<T>* duplicate() {
            if (data) {
                return new GadgetContainerMessage<T>(*this->data);
            } else {
                throw std::runtime_error("Tried to duplicated message after data was taken");
            }
        }

    private:
        std::unique_ptr<Core::TypedMessageChunk<T>> message;
        T* data;
    };


    template<class T>
    GadgetContainerMessage<T> *AsContainerMessage(ACE_Message_Block *mb) {

        if (mb && typeid(*mb) == typeid(GadgetContainerMessage<T>)){
            return reinterpret_cast<GadgetContainerMessage<T>*>(mb);
        }
        return nullptr;
    }
}
