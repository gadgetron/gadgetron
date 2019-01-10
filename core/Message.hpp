#include "Message.h"

#include "GadgetContainerMessage.h"
#include <boost/optional.hpp>
#include <boost/hana.hpp>

#include <iostream>
#include <boost/core/demangle.hpp>
#include "Types.h"

namespace Gadgetron::Core {


    template<class T>
    GadgetContainerMessageBase *TypedMessage<T>::to_container_message() {
        return new GadgetContainerMessage<T>(this->take_data());
    }


    namespace {

        namespace gadgetron_detail {


            namespace hana = boost::hana;

            struct detail {
                template<class Iterator>
                static bool convertible(Iterator it, const Iterator &it_end) {
                    return true;
                }

                template<class Iterator, class T, class ...TYPES>
                static bool convertible(Iterator it, const Iterator &it_end, const hana::basic_type<T> &,
                                        const hana::basic_type<TYPES> &... xs) {
                    if (it == it_end) return false;
                    if (typeid(TypedMessage<T>) == typeid(**it)) {
                        return convertible(++it, it_end, xs...);
                    }
                    return false;
                }

                template<class Iterator, class T, class ...TYPES>
                static bool
                convertible(Iterator it, const Iterator &it_end, const hana::basic_type<optional<T>> &,
                            const hana::basic_type<TYPES> &... xs) {
                    if (convertible(it, it_end, hana::type_c<T>, xs...)) return true;
                    return convertible(it, it_end, xs...);
                }

                template<class Iterator, class... TTYPES>
                static bool
                convertible(Iterator it, const Iterator &it_end, const hana::basic_type<tuple<TTYPES...>> &) {
                    return convertible(it, it_end, hana::type_c<TTYPES>...);
                }


                template<class Iterator, class ...VTYPES, class ...TYPES>
                static bool
                convertible(Iterator it, const Iterator &it_end, const hana::basic_type<variant<VTYPES...>> &) {
                    constexpr auto vtypes = hana::tuple_t<VTYPES...>;
                    return hana::fold(vtypes, false, [&](bool result, auto type) {
                        return result || convertible(it, it_end, type);
                    });
                }


                template<class T>
                static std::unique_ptr<T> reinterpret_message(Message &message) {
                    return static_cast<TypedMessage<T> &>(message).take_data();
                }

                template<class Iterator>
                static auto convert(Iterator& it, const Iterator &it_end) {
                    return hana::make_tuple();
                }

                template<class Iterator, class T, class... TYPES>
                static auto convert(Iterator& it, const Iterator &it_end, const hana::basic_type<T>,
                                    const hana::basic_type<TYPES> &...xs) {
                    auto val = reinterpret_message<T>(**it);
                    return hana::prepend(convert(++it, it_end, xs...), std::move(val));
                }


                template<class Iterator, class T, class... TYPES>
                static auto convert(Iterator& it, const Iterator &it_end, const hana::basic_type<boost::optional<T>> &,
                                    const hana::basic_type<TYPES> &... xs) {

                    if (it != it_end && typeid(TypedMessage<optional<T>>) == typeid(**it)) {
                        auto val = reinterpret_message<optional<T>>(**it);
                        return hana::prepend(convert(++it, it_end, xs...), std::move(val));
                    }

                    if (convertible(it, it_end, hana::basic_type<T>(), xs...)) {
                        auto val = reinterpret_message<T>(**it);
                        return hana::prepend(convert(++it, it_end, xs...), std::make_unique<optional<T>>(std::move(*val)));
                    }

                    return hana::prepend(convert(++it, it_end, xs...), std::make_unique<optional<T>>());
                }

                template<class Iterator, class... TTYPES>
                static auto convert(Iterator it, const Iterator &it_end, const hana::basic_type<tuple<TTYPES...>>&) {

                    if (it != it_end && typeid(TypedMessage<tuple<TTYPES...>>) == typeid(**it)) {
                        auto val = reinterpret_message<tuple<TTYPES...>>(**it);
                        return hana::tuple(std::move(val));
                    }

                    return convert(it,it_end, hana::type_c<TTYPES>... );
                }

                template<class Iterator, class... VTYPES>
                static auto convert(Iterator it, const Iterator &it_end, const hana::basic_type<variant<VTYPES...>> &) {

                    constexpr auto vtypes = hana::tuple_t<VTYPES...>;
                    variant<std::unique_ptr<VTYPES>...> variation;
                    bool variant_found = hana::fold(vtypes, false, [&](bool result, auto type) {
                        if (result) return true;
                        if(convertible(it, it_end, type)) {
                            variation = convert(it, it_end, type);
                            return true;
                        }
                        return false;
                    });
                    if (!variant_found)
                        throw std::runtime_error("Tried to convert message to variant, but no legal type was found");



                }


                template<class Iterator, class... TYPES>
                static auto
                messageTuple_to_tuple(Iterator it, const Iterator &it_end, const hana::basic_type<TYPES> &... xs) {

                    auto result = convert(it, it_end, xs...);
                    return hana::unpack(result, [](auto ...xs) {
                        return std::make_tuple(std::unique_ptr<std::remove_pointer_t<decltype(xs)>>(xs)...);
                    });


                }

                template<class ...ARGS>
                static std::tuple<std::unique_ptr<ARGS>...> message_to_tuple(Message &message) {
                    if (typeid(message) == typeid(MessageTuple)) {
                        auto messages = static_cast<MessageTuple &>(message).take_messages();
                        return messageTuple_to_tuple(messages.begin(), messages.end(), hana::basic_type<ARGS>()...);
                    }
                    auto *m_ptr = &message;
                    return messageTuple_to_tuple(&m_ptr, &m_ptr + 1, hana::basic_type<ARGS>()...);

                }
            };
        }
    }

    template<class ...ARGS>
    bool convertible_to(const Message &message) {
        if (typeid(message) == typeid(MessageTuple)) {
            auto &messages = static_cast<const MessageTuple &>(message).messages();
            return gadgetron_detail::detail::convertible(messages.begin(), messages.end(),
                                                         boost::hana::basic_type<ARGS>()...);
        }

        auto *m_ptr = &message;
        return gadgetron_detail::detail::convertible(&m_ptr, &m_ptr + 1, boost::hana::type<ARGS>()...);
    }

    template<class ...ARGS>
    std::enable_if_t<(sizeof...(ARGS) > 1), std::tuple<std::unique_ptr<ARGS>...>>
    force_unpack(std::unique_ptr<Message> message) {
        return gadgetron_detail::detail::message_to_tuple<ARGS...>(*message);
    }


    template<class T>
    std::unique_ptr<T> force_unpack(std::unique_ptr<Message> message) {
        auto tup = gadgetron_detail::detail::message_to_tuple<T>(*message);
        return std::move(std::get<0>(tup));
    }

    template<class ...ARGS>
    std::enable_if_t<(sizeof...(ARGS) > 1), std::tuple<std::unique_ptr<ARGS>...>>
    unpack(std::unique_ptr<Message> message) {
        if (convertible_to<ARGS...>(*message)) {
            return force_unpack<ARGS...>(std::move(message));
        }
        return std::tuple<std::unique_ptr<ARGS>...>();
    }

    template<class T>
    std::unique_ptr<T> unpack(std::unique_ptr<Message> &&message) {
        if (convertible_to<T>(*message)) {
            return force_unpack<T>(std::move(message));
        }
        return std::unique_ptr<T>();

    }

}
