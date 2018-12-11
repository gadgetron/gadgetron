#include "Message.h"

#include "GadgetContainerMessage.h"

namespace Gadgetron::Core {


    template<class T>
    GadgetContainerMessageBase *TypedMessage<T>::to_container_message() {
        return new GadgetContainerMessage<T>(this->take_data());
    }


    namespace {

        namespace gadgetron_detail {




            template<unsigned int I>
            bool convertible_to_impl(const MessageTuple &) {
                return true;
            }
            template<unsigned int I, class T, class ...REST>
            bool convertible_to_impl(const MessageTuple &messagetuple) {
                auto &ptr = messagetuple.messages()[I];
                if (typeid(*ptr) == typeid(T)) {
                    return convertible_to_impl<I + 1, REST...>(messagetuple);
                }

                return false;

            }


            template<class ...ARGS>
            bool convertible_to_impl(const MessageTuple &messageTuple) {
                if (sizeof...(ARGS) <= messageTuple.messages().size()) {
                    return convertible_to_impl<0, ARGS...>(messageTuple);
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
    bool convertible_to(const Message &message) {
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