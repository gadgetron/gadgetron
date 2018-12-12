#include "Message.h"

#include "GadgetContainerMessage.h"
#include <boost/optional.hpp>

namespace Gadgetron::Core {


    template<class T>
    GadgetContainerMessageBase *TypedMessage<T>::to_container_message() {
        return new GadgetContainerMessage<T>(this->take_data());
    }


    namespace {

        namespace gadgetron_detail {

            template<class ... ARGS>
            struct count_non_optional;

            template<>
            struct count_non_optional<>{
                static constexpr size_t value = 0;
            };

            template<class T, class ...ARGS>
            struct count_non_optional<T,ARGS...>{
                static constexpr size_t value = count_non_optional<ARGS...>::value + 1;
            };

            template<class T, class ...ARGS>
            struct count_non_optional<boost::optional<T>,ARGS...>{
                static constexpr size_t value = count_non_optional<ARGS...>::value + 1;
            };


            template<unsigned int I, class ...ARGS>
            struct tuple_converter;

            template<unsigned int I>
            struct tuple_converter<I> {
                static bool accept(const MessageTuple&){
                    return true;
                }
            };



            template<unsigned int I, class T, class ...REST>
            struct tuple_converter<I,T,REST...> {
                static bool accept(const MessageTuple &messagetuple) {
                    auto &ptr = messagetuple.messages()[I];
                    if (typeid(*ptr) == typeid(T)) {
                        return tuple_converter< I + 1, REST...>::accept(messagetuple);
                    }
                    return false;
                }
            };

            template<unsigned int I, class T, class ...REST>
            struct tuple_converter<I,boost::optional<T>,REST...> {
                static bool accept(const MessageTuple &messagetuple) {
                    auto &ptr = messagetuple.messages()[I];
                    if (typeid(*ptr) == typeid(T)) {
                        return tuple_converter< I + 1, REST...>::accept(messagetuple);
                    }
                    return tuple_converter<I,REST...>::accept(messagetuple);
                }
            };



            template<class ...ARGS>
            bool convertible_to_impl(const MessageTuple &messageTuple) {
                if (count_non_optional<ARGS...>::value <= messageTuple.messages().size()) {
                    return tuple_converter<0, ARGS...>::accept(messageTuple);
                }

                return false;
            }

            template<class T>
            std::unique_ptr<T> reinterpret_message(std::unique_ptr<Message> &message) {
                return static_cast<TypedMessage<T> *>(message.get())->take_data();
            }

            template<class ...ARGS>
            struct tuple_maker {
                template<size_t ...S>
                static auto from_messages(std::vector<std::unique_ptr<Message>> &messages, std::index_sequence<S...>) {
                    return std::make_tuple(reinterpret_message<ARGS>(messages[S])...);
                }
            };

            template<class ...ARGS, typename Indices = std::make_index_sequence<sizeof...(ARGS)>>
            std::tuple<std::unique_ptr<ARGS>...> messageTuple_to_tuple(MessageTuple &messageTuple) {
                auto messages = messageTuple.take_messages();
                return tuple_maker<ARGS...>::from_messages(messages, Indices{});
            }

            template<class ...ARGS>
            std::tuple<std::unique_ptr<ARGS>...> message_to_tuple(Message &message) {
                return gadgetron_detail::messageTuple_to_tuple<ARGS...>(static_cast<MessageTuple&>(message));
            }
        }
    }

    template<class ...REST>
    std::enable_if_t<(sizeof...(REST) > 1), bool> convertible_to(const Message &message) {
        if (typeid(message) == typeid(MessageTuple)) {
            return gadgetron_detail::convertible_to_impl<REST...>(static_cast<const MessageTuple&>(message));
        }

        return false;
    }

    template<class T>
    bool convertible_to(const Message &message) {
        return typeid(message) == typeid(TypedMessage<T>);
    }

    template<class ...ARGS>
    std::enable_if_t<(sizeof...(ARGS) > 1),std::tuple<std::unique_ptr<ARGS>...>> force_unpack(std::unique_ptr<Message> &message) {
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