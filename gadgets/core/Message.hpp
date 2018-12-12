#include "Message.h"

#include "GadgetContainerMessage.h"
#include <boost/optional.hpp>
#include <boost/hana.hpp>

namespace Gadgetron::Core {


    template<class T>
    GadgetContainerMessageBase *TypedMessage<T>::to_container_message() {
        return new GadgetContainerMessage<T>(this->take_data());
    }


    namespace {

        namespace gadgetron_detail {

            namespace hana = boost::hana;

              template<class T>
              constexpr bool is_optional(hana::basic_type<T> const&){
                  return true;
              }

            template<class T>
            constexpr auto is_optional(hana::basic_type<boost::optional<T>> const&){
                return hana::bool_c<true>;
            }





            template<class T, class Iterator>
            auto convertible(const hana::basic_type<T> &, Iterator it) {
                if (typeid(TypedMessage < T > ) == typeid(*it))
                    return hana::make_tuple(true, ++it);
                return hana::make_tuple(false, ++it);
            }

            template<class T, class Iterator>
            auto convertible(const hana::basic_type<boost::optional<T>> &, Iterator it) {
                if (typeid(TypedMessage < T > ) == typeid(*it))
                    return hana::make_tuple(true, ++it);
                return hana::make_tuple(true, it);
            }

            template<class TYPES>
            bool convertible_to_impl(const TYPES &types, const MessageTuple &messageTuple) {
                using namespace hana::literals;
                if (hana::count_if(types, [](auto a){ return hana::not_(is_optional(a));}) <= messageTuple.messages().size()) {
                    auto val = hana::fold_left(types,hana::make_tuple(true,messageTuple.messages().begin()),[&](auto tuple, const auto& specific_type ){
                        if (tuple[1_c] == messageTuple.messages().end())
                            return tuple;
                        auto result =  convertible(specific_type,tuple[1_c]);
                        return hana::make_tuple(tuple[0_c] && result[0_c],result[1_c]);
                    });
                    return val[0_c];
                }

                return false;
            }

            template<class T>
            T* reinterpret_message(std::unique_ptr<Message> &message) {
                return static_cast<TypedMessage<T> *>(message.get())->take_data().release();
            }


            template<class T, class Iterator>
            auto convert(const hana::basic_type<T> &, Iterator it) {
                return hana::make_tuple(reinterpret_message<T>(*it), ++it);
            }

            template<class T, class Iterator>
            auto convert(const hana::basic_type<boost::optional<T>> &, Iterator it) {
                if (typeid(TypedMessage < T > ) == typeid(*it)) {
                    auto content = std::unique_ptr<T>(reinterpret_message<T>(*it));
                    return hana::make_tuple(new boost::optional<T>(std::move(*content)),it);
                }
                return hana::make_tuple(new boost::optional<T>(boost::none), it);
            }


            template<class TYPES>
            auto messageTuple_to_tuple(const TYPES &types, MessageTuple &messageTuple) {

                using namespace hana::literals;
                auto messages = messageTuple.take_messages();
                auto result = hana::fold_left(types, hana::make_tuple(hana::make_tuple(),messages.begin()),
                [](auto tuple,auto& specific_type ){
                    auto result = convert(specific_type,tuple[1_c]);
                    return hana::make_tuple(hana::append(tuple[0_c],result[0_c]), result[1_c]);
                });

                return hana::unpack(result[0_c],[](auto ...xs){
                    return std::make_tuple(std::unique_ptr<std::remove_pointer_t<decltype(xs)>>(xs)...);
                });


            }

            template<class ...ARGS>
            std::tuple<std::unique_ptr<ARGS>...> message_to_tuple(Message &message) {
                return gadgetron_detail::messageTuple_to_tuple(hana::tuple_t<ARGS...>,static_cast<MessageTuple&>(message));
            }
        }
    }

    template<class ...ARGS>
    bool convertible_to(const Message &message) {
        if (typeid(message) == typeid(MessageTuple)) {
            return gadgetron_detail::convertible_to_impl(boost::hana::tuple_t<ARGS...>,
                                                         static_cast<const MessageTuple &>(message));
        }

        return false;
    }

    template<class T>
    bool convertible_to(const Message &message) {
        return typeid(message) == typeid(TypedMessage<T>);
    }

    template<class ...ARGS>
    std::tuple<std::unique_ptr<ARGS>...> force_unpack(std::unique_ptr<Message> &message) {
        return gadgetron_detail::message_to_tuple<ARGS...>(*message);
    }


    template<class T>
    std::unique_ptr<T> force_unpack(std::unique_ptr<Message> &message) {
        return gadgetron_detail::reinterpret_message<T>(message);
    }

    template<class ...ARGS>
    std::tuple<std::unique_ptr<ARGS>...> unpack(std::unique_ptr<Message> &message) {
        if (convertible_to<ARGS...>(*message)) {
            return force_unpack<ARGS...>(message);
        }
    }

    template<class T>
    std::unique_ptr<T> unpack(std::unique_ptr<Message> &message) {
        if (convertible_to<T>(*message)) {
            return gadgetron_detail::reinterpret_message<T>(message);
        }
        return std::unique_ptr<T>();

    }

}