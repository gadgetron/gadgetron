#pragma once

#include <string>
#include "Message.h"
#include "LegacyACE.h"

namespace Gadgetron {


    class GadgetContainerMessageBase : public ACE_Message_Block {
    public:


        virtual std::unique_ptr<Core::Message> take_message() = 0;

        virtual ~GadgetContainerMessageBase(){};


        std::type_index type(){
            return message->type();
        }

    protected:


        std::unique_ptr<Core::Message> message;
    };



    template<class T>
    class GadgetContainerMessage : public GadgetContainerMessageBase {

    public:
        /**
         *  Constructor, passing on input arguments to the contained class.
         * @param xs Variadic arguments to the contained class
         */

        template<typename... X>
        GadgetContainerMessage(X &&... xs)  {
            data= new T(std::move(xs)...);
            message = std::make_unique<Core::TypedMessage<T>>(std::unique_ptr<T>(data));
        }

        GadgetContainerMessage(std::unique_ptr<Core::TypedMessage<T>>&& msg ){
            auto msg_data = msg->take_data();
            data = msg_data.get();
            message = std::make_unique<Core::TypedMessage<T>>(std::move(msg_data));
        }

        GadgetContainerMessage(std::unique_ptr<T>&& msg_data ){
            data = msg_data.get();
            message = std::make_unique<Core::TypedMessage<T>>(std::move(msg_data));
        }

        virtual ~GadgetContainerMessage() {

        }

        virtual std::unique_ptr<Core::Message> take_message() override {
            data = nullptr;
            return std::unique_ptr<Core::TypedMessage<T>>(std::move(message));
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

        T* data;
    };


    template<class T>
    GadgetContainerMessage<T> *AsContainerMessage(ACE_Message_Block *mb) {

        if (typeid(mb) == typeid(GadgetContainerMessage<T>*)){
            return reinterpret_cast<GadgetContainerMessage<T>*>(mb);
        }
        return nullptr;
    }
}
