
#include "LegacyACE.h"

namespace Gadgetron {

    ACE_Message_Queue<ACE_MT_SYNCH>::ACE_Message_Queue(int, int) {

    }

    ACE_Message_Queue<ACE_MT_SYNCH>::ACE_Message_Queue() {

    }

    size_t ACE_Message_Queue<ACE_MT_SYNCH>::message_count() {
        return queue.size();
    }

    int ACE_Message_Queue<ACE_MT_SYNCH>::dequeue_head(ACE_Message_Block *&msg) {
        msg = queue.front().release();
        queue.pop_front();
        return 0;
    }

    int ACE_Message_Queue<ACE_MT_SYNCH>::dequeue_tail(ACE_Message_Block *&msg) {
        msg = queue.back().release();
        queue.pop_back();
        return 0;
    }

    int ACE_Message_Queue<ACE_MT_SYNCH>::enqueue_tail(ACE_Message_Block *msg) {
        queue.emplace_back(msg);
        return 0;
    }

    int ACE_Message_Queue<ACE_MT_SYNCH>::enqueue_head(ACE_Message_Block *msg) {
        queue.emplace_front(msg);
        return 0;
    }

    void ACE_Message_Queue<ACE_MT_SYNCH>::flush() {
        queue = std::list<std::unique_ptr<ACE_Message_Block>>();
    }


    ACE_Message_Queue<ACE_MT_SYNCH>::ITERATOR::ITERATOR(ACE_Message_Queue &q) : list_iterator(q.queue.begin()),
                                                                                list_end(q.queue.end()) {

    }

    int ACE_Message_Queue<ACE_MT_SYNCH>::ITERATOR::advance() {
        list_iterator++;
        if (list_iterator == list_end) { return 0; }
        return 1;
    }

    int ACE_Message_Queue<ACE_MT_SYNCH>::ITERATOR::next(ACE_Message_Block *& msg) {
        if (list_iterator == list_end) { return 0; }
        msg = list_iterator->get();
        return 1;
    }
}
