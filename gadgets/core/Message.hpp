#include "Message.h"

#include "GadgetContainerMessage.h"
#include <boost/optional.hpp>
#include <boost/hana.hpp>

#include <iostream>
#include <boost/core/demangle.hpp>
namespace Gadgetron::Core {


    template<class T>
    GadgetContainerMessageBase *TypedMessage<T>::to_container_message() {
        return new GadgetContainerMessage<T>(this->take_data());
    }


    namespace {

        namespace gadgetron_detail {

            namespace hana = boost::hana;

            template<class T>
            constexpr bool is_optional(hana::basic_type<T> const &) {
                return true;
            }

            template<class T>
            constexpr auto is_optional(hana::basic_type<boost::optional<T>> const &) {
                return hana::bool_c<true>;
            }


            template<class T, class Iterator>
            auto convertible(const hana::basic_type<T> &, Iterator it, Iterator it_end) {
                if (it == it_end) return hana::make_tuple(false, it);
                if (typeid(TypedMessage < T > ) == typeid(**it))
                    return hana::make_tuple(true, ++it);
                return hana::make_tuple(false, ++it);
            }

            template<class T, class Iterator>
            auto convertible(const hana::basic_type<boost::optional<T>> &, Iterator it, Iterator it_end) {
                if (it != it_end && typeid(TypedMessage < T > ) == typeid(**it))
                    return hana::make_tuple(true, ++it);
                return hana::make_tuple(true, it);
            }

            template<class TYPES, class Iterator>
            bool convertible_to_impl(const TYPES &types,const  Iterator& it, const Iterator& it_end) {
                using namespace hana::literals;
                if (hana::count_if(types,
                                   [](auto a) { return hana::not_(is_optional(a)); }) <=
                    std::distance(it,it_end)) {

                    auto val = hana::fold_left(types, hana::make_tuple(true, it),
                                               [&](auto tuple, const auto &specific_type) {
                                                   auto result = convertible(specific_type, tuple[1_c],
                                                                             it_end);
                                                   return hana::make_tuple(tuple[0_c] && result[0_c], result[1_c]);
                                               });
                    return val[0_c];
                }

                return false;
            }

            template<class T>
            T *reinterpret_message(Message &message) {
                return static_cast<TypedMessage<T>&>(message).take_data().release();
            }


            template<class T, class Iterator>
            auto convert(const hana::basic_type<T> &, Iterator it, Iterator end_it) {
                auto val = reinterpret_message<T>(**it);
                return hana::make_tuple(val, ++it);
            }

            template<class T, class Iterator>
            auto convert(const hana::basic_type<boost::optional<T>> &, Iterator it, Iterator end_it) {
                if (it != end_it && typeid(TypedMessage < T > ) == typeid(*it)) {
                    auto content = std::unique_ptr<T>(reinterpret_message<T>(**it));
                    return hana::make_tuple(new boost::optional<T>(std::move(*content)), it);
                }
                return hana::make_tuple(new boost::optional<T>(boost::none), it);
            }


            template<class TYPES, class Iterator>
            auto messageTuple_to_tuple(const TYPES &types, const Iterator& it, const Iterator& it_end) {

                using namespace hana::literals;
                auto result = hana::fold_left(types, hana::make_tuple(hana::make_tuple(), it),
                                              [&](auto tuple, auto &specific_type) {
                                                  auto result = convert(specific_type, tuple[1_c],it_end);
                                                  return hana::make_tuple(hana::append(tuple[0_c], result[0_c]),
                                                                          result[1_c]);
                                              });

                return hana::unpack(result[0_c], [](auto ...xs) {
                    return std::make_tuple(std::unique_ptr<std::remove_pointer_t<decltype(xs)>>(xs)...);
                });


            }

            template<class ...ARGS>
            std::tuple<std::unique_ptr<ARGS>...> message_to_tuple(Message &message) {
                 if (typeid(message) == typeid(MessageTuple)) {
                    auto messages = static_cast<MessageTuple&>(message).take_messages();
                    return gadgetron_detail::messageTuple_to_tuple(hana::tuple_t<ARGS...>,
                                                           messages.begin(),messages.end());
                 }
                 auto* m_ptr = &message;
                 return gadgetron_detail::messageTuple_to_tuple(hana::tuple_t<ARGS...>,&m_ptr,&m_ptr+1);

            }
        }
    }

    template<class ...ARGS>
    std::enable_if_t<(sizeof...(ARGS) > 1), bool> convertible_to(const Message &message) {
        if (typeid(message) == typeid(MessageTuple)) {
            auto& messages  = static_cast<const MessageTuple &>(message).messages();
            return gadgetron_detail::convertible_to_impl(boost::hana::tuple_t<ARGS...>,
                                                         messages.begin(),messages.end());
        }

        auto* m_ptr = &message;
        return gadgetron_detail::convertible_to_impl(boost::hana::tuple_t<ARGS...>,&m_ptr,&m_ptr+1);
    }

    template<class T>
    bool convertible_to(const Message &message) {
        return typeid(message) == typeid(TypedMessage<T>);
    }

    template<class ...ARGS>
    std::enable_if_t<(sizeof...(ARGS) > 1), std::tuple<std::unique_ptr<ARGS>...>>
    force_unpack(std::unique_ptr<Message> &message) {
        return gadgetron_detail::message_to_tuple<ARGS...>(*message);
    }


    template<class T>
    std::unique_ptr<T> force_unpack(std::unique_ptr<Message> message) {
        return std::unique_ptr<T>(gadgetron_detail::reinterpret_message<T>(*message));
    }

    template<class ...ARGS>
    std::enable_if_t<(sizeof...(ARGS) > 1),std::tuple<std::unique_ptr<ARGS>...>> unpack(std::unique_ptr<Message> message) {
        if (convertible_to<ARGS...>(*message)) {
            return force_unpack<ARGS...>(message);
        }
        return std::tuple<std::unique_ptr<ARGS>...>();
    }

    template<class T>
    std::unique_ptr<T> unpack(std::unique_ptr<Message>&& message) {
        if (convertible_to<T>(*message)) {
            return force_unpack<T>(std::move(message));
        }
        return std::unique_ptr<T>();

    }

}
