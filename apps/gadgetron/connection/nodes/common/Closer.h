#pragma once

#include <memory>

#include "log.h"

namespace Gadgetron::Server::Connection::Nodes {

    template<class CLOSED>
    class Closer{
    public:
        ~Closer() { closed->close(); };
        explicit Closer(CLOSED& closed) : closed(closed) {}
    private:
        CLOSED& closed;
    };

    template<class CLOSED>
    static Closer<CLOSED> make_closer(CLOSED& cls) { return Closer<CLOSED>(cls); }

}
