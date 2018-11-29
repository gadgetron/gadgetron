#pragma once

#include <string>
#include "Message.h"

namespace Gadgetron {


    class ACE_Message_Block {

    public:
        ACE_Message_Block(){};
        ACE_Message_Block(const std::string &s) : buffer(s) {

        }

        const char *rd_ptr() { return buffer.c_str(); };

        virtual ~ACE_Message_Block() {};
    private:
        std::string buffer;
    };


    class GadgetContainerMessageBase : public ACE_Message_Block {
    public:
        GadgetContainerMessageBase(){};

        virtual ~GadgetContainerMessageBase() {};

        virtual void *release() {
            delete (this); // Seppuku
            return nullptr;
        }


    };



    template<class T>
    class GadgetContainerMessage : public GadgetContainerMessageBase, private Core::TypedMessage<T> {

    public:
        /**
         *  Constructor, passing on input arguments to the contained class.
         * @param xs Variadic arguments to the contained class
         */
        template<typename... X>
        GadgetContainerMessage(X... xs) : Core::TypedMessage<T>(T(xs...)) {

        }

        virtual ~GadgetContainerMessage() {
        }

        T *getObjectPtr() {
            return this->data.get();
        }

        GadgetContainerMessageBase* cont(){return cont_element;}
        void  cont(GadgetContainerMessageBase* ptr){cont_element = ptr;}


        GadgetContainerMessage<T>* duplicate() {
            return new GadgetContainerMessage<T>(*this->data);
        }


    private:

        GadgetContainerMessageBase* cont_element = nullptr;



    };


    template<class T>
    GadgetContainerMessage<T> *AsContainerMessage(ACE_Message_Block *mb) {
        if (typeid(mb) == typeid(GadgetContainerMessage<T> *))
            return reinterpret_cast<GadgetContainerMessage<T> * >(mb);
        return nullptr;
    }
}
