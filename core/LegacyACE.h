#pragma once

#include <list>
#include <memory>
#include <string>

namespace Gadgetron {

    class ACE_Message_Block {

    public:
        ACE_Message_Block() = default;

        explicit ACE_Message_Block(std::string s) : buffer(std::move(s)) {

        }

        const char *rd_ptr() { return buffer.c_str(); };

        virtual ~ACE_Message_Block() {
            if (cont_element) {
                cont_element->release();
            }
        };

        virtual void *release() {
            delete (this); // Seppuku
            return nullptr;
        }

        ACE_Message_Block *cont() { return cont_element; }

        void cont(ACE_Message_Block *ptr) { cont_element = ptr; }


    private:
        ACE_Message_Block *cont_element = nullptr;

    private:
        std::string buffer;
    };




    struct ACE_MT_SYNCH {
    };

    namespace ACE_Message_Queue_Base {
        constexpr int DEFAULT_HWM = 0;
        constexpr int DEFAULT_LWM = 0;
    }

    template<class T>
    class ACE_Message_Queue {
        virtual void not_implemented_ever() = 0;
    };


    template<>
    class ACE_Message_Queue<ACE_MT_SYNCH> {

    public:

        [[deprecated]]
        ACE_Message_Queue(int, int);

        [[deprecated]]
        ACE_Message_Queue();

        void high_water_mark(size_t bsize) {};

        void low_water_mark(size_t bsize) {};

        size_t message_count();

        int dequeue_head(ACE_Message_Block *&);

        int dequeue_tail(ACE_Message_Block *&);

        int enqueue_tail(ACE_Message_Block *);

        int enqueue_head(ACE_Message_Block *);

        void flush();

    private:

        std::list<std::unique_ptr<ACE_Message_Block>> queue;


    public:
        class ITERATOR {
        public:
            explicit ITERATOR(ACE_Message_Queue &);

            int advance();

            int next(ACE_Message_Block *&);

        private:

            decltype(queue.begin()) list_iterator;
            decltype(queue.begin()) list_end;
        };


        friend ITERATOR;

    };


}