#pragma once

#include <memory>
#include <typeindex>
#include <boost/dll.hpp>

#include "Message.h"
#include "Channel.h"
#include "Node.h"

namespace Gadgetron::Core {

    class Writer {
    public:
        virtual ~Writer() = default;

        virtual bool accepts(const Message &) = 0;

        virtual void write(std::ostream &stream, Message message) = 0;
    };


    template<class ...ARGS>
    class TypedWriter : public Writer {
    public:
        ~TypedWriter() override = default;

        bool accepts(const Message &) override;

        void write(std::ostream &stream, Message message) override;

    protected:
        virtual void serialize(std::ostream &stream, const ARGS& ...) = 0;
    };

}

namespace Gadgetron::Core {

    template<class ...ARGS>
    bool TypedWriter<ARGS...>::accepts(const Message &message) {
        return convertible_to<ARGS...>(message);
    }


    namespace {
        namespace gadgetron_writer_detail {
            template<class F, size_t... Is>
            constexpr auto index_apply_impl(F f,
                                            std::index_sequence<Is...>) {
                return f(std::integral_constant<size_t, Is>{}...);
            }

            template<size_t N, class F>
            constexpr auto index_apply(F f) {
                return index_apply_impl(f, std::make_index_sequence<N>{});
            }
        }
    }
}



    template<class ...ARGS>
    void Gadgetron::Core::TypedWriter<ARGS...>::write(std::ostream &stream,  Message message) {

        std::tuple<ARGS...> arg_tuple = force_unpack<ARGS...>(std::move(message));

        gadgetron_writer_detail::index_apply<sizeof...(ARGS)>(
                [&](auto... Is) { this->serialize(stream, std::move(std::get<Is>(arg_tuple))...); });
    }


#include "Writer.hpp"

#include <boost/dll/alias.hpp>
#define GADGETRON_WRITER_EXPORT(WriterClass)                        \
std::unique_ptr<Gadgetron::Core::Writer> writer_factory_##WriterClass() {            \
    return std::make_unique<WriterClass>();                         \
}                                                                   \
                                                                    \
BOOST_DLL_ALIAS(                                                    \
        writer_factory_##WriterClass,                               \
        writer_factory_export_##WriterClass                         \
)
