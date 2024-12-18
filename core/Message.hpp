#include "Message.h"

#include <boost/hana.hpp>
#include <boost/core/demangle.hpp>
#include <iostream>


namespace Gadgetron::Core {


    template<class T>
    GadgetContainerMessageBase* TypedMessageChunk<T>::to_container_message() {
        return new GadgetContainerMessage<T>(std::move(data));
    }


    template<class T>
    std::unique_ptr<MessageChunk> TypedMessageChunk<T>::clone() const {
        return std::make_unique<TypedMessageChunk<T>>(data);
    }

    namespace {
        namespace gadgetron_message_detail {

            template<class T>
            std::unique_ptr<MessageChunk> make_message(T &&input) {
                return std::make_unique<TypedMessageChunk<std::remove_reference_t<T>>>(std::forward<T>(input));
            }

            struct MessageMaker {
                template<class... VARGS, class ...REST>
                static void
                add_messages(std::vector<std::unique_ptr<MessageChunk>> &messages, std::variant<VARGS...> var,
                            REST &&... args) {
                    std::visit([&](auto val) { add_messages(messages, val, std::forward<REST>(args)...); },
                                         var);
                }

                template<class... TARGS, class ...REST>
                static void
                add_messages(std::vector<std::unique_ptr<MessageChunk>> &messages, std::tuple<TARGS...> opt,
                             REST &&... args) {
                    std::apply([&](auto... targs) {
                            (..., add_messages(messages, std::move(targs), std::forward<REST...>(args)...));
                        }, opt);
                }

                template<class T, class ...REST>
                static void
                add_messages(std::vector<std::unique_ptr<MessageChunk>> &messages, std::optional <T> opt, REST &&... args) {
                    if (opt) messages.emplace_back(make_message(*opt));
                    add_messages(messages, std::forward<REST>(args)...);
                }

                template<class T, class ...REST>
                static void
                add_messages(std::vector<std::unique_ptr<MessageChunk>> &messages, T &&arg, REST &&... args) {
                    messages.emplace_back(make_message(std::forward<T>(arg)));
                    add_messages(messages, std::forward<REST>(args)...);
                }

                static void add_messages(std::vector<std::unique_ptr<MessageChunk>> &messages) {
                }

            };

            template<class... ARGS>
            std::vector<std::unique_ptr<MessageChunk>> make_messages(ARGS... args) {
                std::vector<std::unique_ptr<MessageChunk>> messages;
                MessageMaker::add_messages(messages, std::forward<ARGS>(args)...);
                return messages;
            }
        }
    }


    namespace {

        namespace gadgetron_message_detail {


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
                    if (typeid(TypedMessageChunk<T>) == typeid(**it)) {
                        return convertible(++it, it_end, xs...);
                    }
                    return false;
                }


                template<class T, class S>
                static hana::tuple<T, S> combine(T &&val1, S &&val2) {
                    return hana::tuple<T, S>(std::forward<T>(val1), std::forward<S>(val2));
                }

                template<class T, class... SARGS>
                static hana::tuple<T, SARGS...> combine(T &&val1, hana::tuple<SARGS...> &&val2) {
                    return hana::prepend(val2, std::forward<T>(val1));
                }


                template<class Iterator, class T, class ...TYPES>
                static bool
                convertible(Iterator it, const Iterator &it_end, const hana::basic_type<std::optional < T>>&,
                const hana::basic_type<TYPES> &... xs
                ) {
                    if (convertible(it, it_end, hana::type_c<T>, xs...)) return true;
                    return convertible(it, it_end, xs...);
                }

                template<class Iterator, class... TTYPES>
                static bool
                convertible(Iterator it, const Iterator &it_end, const hana::basic_type<std::tuple < TTYPES...>>&) {
                    return convertible(it, it_end, hana::type_c<TTYPES>...);
                }


                template<class Iterator, class ...VTYPES >
                static bool
                convertible(Iterator it, const Iterator &it_end, const hana::basic_type<std::variant < VTYPES...>>&) {
                    constexpr auto vtypes = hana::tuple_t<VTYPES...>;
                    return hana::fold(vtypes, false, [&](bool result, auto type) {
                        return result || convertible(it, it_end, type);
                    });
                }


                template<class T>
                static TypedMessageChunk <T> &reinterpret_message(MessageChunk &message) {
                    return static_cast<TypedMessageChunk<T> &>(message);
                }

                template<class Iterator, class T>
                static T convert(Iterator &it, const Iterator &it_end, const hana::basic_type<T>&) {
                    return std::move(reinterpret_message<T>(**it).data);
                }

                template<class Iterator, class T>
                static std::optional <T> convert(Iterator &it, const Iterator &it_end, const hana::basic_type<std::optional < T>>

                ) {
                    if (convertible(it, it_end, hana::type_c<T>)) return reinterpret_message<T>(**it).data;
                    return std::optional<T>();
                }

                template<class Iterator, class T, class... TYPES>
                static hana::tuple<T, TYPES...> convert(Iterator &it, const Iterator &it_end, const hana::basic_type<T>&,
                                                        const hana::basic_type<TYPES> &...xs) {
                    auto &value = reinterpret_message<T>(**it).data;
                    return combine(std::move(value), convert(++it, it_end, xs...));
                }


                template<class Iterator, class T, class... TYPES>
                static hana::tuple<std::optional < T>, TYPES...>
                convert(Iterator & it,  const Iterator &it_end, const hana::basic_type<std::optional<T>> &,
                    const hana::basic_type<TYPES> &... xs
                ) {

                    if (convertible(it, it_end, hana::basic_type<T>(), xs...)) {
                        auto &val = reinterpret_message<T>(**it).data;
                        return combine(std::optional<T>(std::move(val)), convert(++it, it_end, xs...));
                    }
                    return combine(std::optional<T>(), convert(it, it_end, xs...));
                }

                template<class Iterator, class... TTYPES>
                static std::tuple<TTYPES...>
                convert(Iterator it, const Iterator &it_end, const hana::basic_type<std::tuple < TTYPES...>>&) {
                    return hana::unpack(convert(it, it_end, hana::type_c<TTYPES>...), [](auto ...xs) {
                        return std::make_tuple(std::move(xs)...);
                    });
                }

                template<class Iterator, class... VTYPES>
                static std::variant<VTYPES...>
                convert(Iterator it, const Iterator &it_end, const hana::basic_type<std::variant < VTYPES...>>&) {

                    constexpr auto vtypes = hana::tuple_t<VTYPES...>;
                    std::variant < VTYPES...> variation;
                    bool variant_found = hana::fold(vtypes, false, [&](bool result, auto type) {
                        if (result) return true;
                        if (convertible(it, it_end, type)) {
                            variation = convert(it, it_end, type);
                            return true;
                        }
                        return false;
                    });
                    if (!variant_found)
                        throw std::runtime_error("Tried to convert message to variant, but no legal type was found");


                    return variation;
                }


                template<class Iterator, class... TYPES>
                static auto
                messageTuple_to_tuple(std::enable_if_t<(sizeof...(TYPES) > 1), Iterator> it, const Iterator &it_end,
                                      const hana::basic_type<TYPES> &... xs) {

                    auto result = convert(it, it_end, xs...);
                    return hana::unpack(result, [](auto ...xs) {
                        return std::make_tuple(std::move(xs)...);
                    });


                }

                template<class Iterator, class TYPE>
                static auto
                messageTuple_to_tuple(Iterator it, const Iterator &it_end, const hana::basic_type<TYPE> &x) {
                    return convert(it, it_end, x);
                }

                template<class ...ARGS>
                static auto message_to_tuple(Message &message) {
                    auto messages = message.take_messages();
                    return messageTuple_to_tuple(messages.begin(), messages.end(), hana::basic_type<ARGS>()...);

                }
            };
        }
    }

    template<class ...ARGS>
    bool convertible_to(const Message &message) {
        auto &messages = message.messages();
        return gadgetron_message_detail::detail::convertible(messages.begin(), messages.end(),
                                                     boost::hana::basic_type<ARGS>()...);

    }

    template<class ...ARGS>
    std::enable_if_t<(sizeof...(ARGS) > 1), std::tuple<ARGS...>>
    force_unpack(Message message) {
        return gadgetron_message_detail::detail::message_to_tuple<ARGS...>(message);
    }


    template<class T>
    T force_unpack(Message message) {
        return gadgetron_message_detail::detail::message_to_tuple<T>(message);
    }

    template<class ...ARGS>
    std::enable_if_t<(sizeof...(ARGS) > 1), std::optional < std::tuple<ARGS...>>>
    unpack(Message &&message) {
    if (convertible_to<ARGS...>(message)) {
        return force_unpack<ARGS...>(std::move(message));
    }
    return std::nullopt;
}

template<class T>
std::optional <T> unpack(Message &&message) {
    if (convertible_to<T>(message)) {
        return force_unpack<T>(std::move(message));
    }
    return std::nullopt;

}


}

template<class... ARGS>
Gadgetron::Core::Message::Message(ARGS &&... args) : messages_(
        gadgetron_message_detail::make_messages<ARGS...>(std::forward<ARGS>(args)...)) {


}

