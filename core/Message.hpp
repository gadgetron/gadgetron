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
        return new GadgetContainerMessage<T>(std::move(data));
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


                template<class T, class S>
                static hana::tuple<T,S> combine(T&& val1, S&& val2){
                    return hana::tuple<T,S>(std::forward<T>(val1),std::forward<S>(val2));
                }

                template<class T, class... SARGS>
                static hana::tuple<T,SARGS...> combine(T&& val1, hana::tuple<SARGS...>&& val2){
                    return hana::prepend(val2, std::forward<T>(val1));
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
                static TypedMessage<T>& reinterpret_message(Message &message) {
                    return static_cast<TypedMessage<T> &>(message);
                }

                template<class Iterator,class T>
                static T convert(Iterator& it, const Iterator &it_end, const hana::basic_type<T>) {
                    return std::move(reinterpret_message<T>(**it).data);
                }

                template<class Iterator,class T>
                static optional<T> convert(Iterator& it, const Iterator &it_end, const hana::basic_type<optional<T>>) {
                    if (convertible(it,it_end,hana::type_c<T>)) return reinterpret_message<T>(**it).data;
                    return optional<T>();
                }

                template<class Iterator, class T, class... TYPES>
                static hana::tuple<T,TYPES...> convert(Iterator& it, const Iterator &it_end, const hana::basic_type<T>,
                                    const hana::basic_type<TYPES> &...xs) {
                    auto& value = reinterpret_message<T>(**it).data;
                    return combine(std::move(value),convert(++it, it_end, xs...));
                }


                template<class Iterator, class T, class... TYPES>
                static hana::tuple<optional<T>,TYPES...> convert(Iterator& it, const Iterator &it_end, const hana::basic_type<boost::optional<T>> &,
                                    const hana::basic_type<TYPES> &... xs) {

                   if (convertible(it, it_end, hana::basic_type<T>(), xs...)) {
                        auto& val = reinterpret_message<T>(**it).data;
                        return combine(optional<T>(std::move(val)), convert(++it, it_end, xs...));
                    }
                    return combine(optional<T>(),convert(it,it_end,xs...));
                }

                template<class Iterator, class... TTYPES>
                static tuple<TTYPES...> convert(Iterator it, const Iterator &it_end, const hana::basic_type<tuple<TTYPES...>>&) {
                    return hana::unpack(convert(it,it_end, hana::type_c<TTYPES>... ), [](auto ...xs){
                        return std::make_tuple(std::move(xs)...);
                    });
                }

                template<class Iterator, class... VTYPES>
                static variant<VTYPES...> convert(Iterator it, const Iterator &it_end, const hana::basic_type<variant<VTYPES...>> &) {

                    constexpr auto vtypes = hana::tuple_t<VTYPES...>;
                    variant<VTYPES...> variation;
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


                    return variation;
                }


                template<class Iterator, class... TYPES>
                static auto
                messageTuple_to_tuple(std::enable_if_t<(sizeof...(TYPES) > 1),Iterator> it, const Iterator &it_end, const hana::basic_type<TYPES> &... xs) {

                    auto result = convert(it, it_end, xs...);
                    return hana::unpack(result, [](auto ...xs) {
                        return std::make_tuple(std::move(xs)...);
                    });


                }

                template<class Iterator, class TYPE>
                static auto
                messageTuple_to_tuple(Iterator it, const Iterator &it_end, const hana::basic_type<TYPE> & x) {
                    return convert(it, it_end, x);
                }

                template<class ...ARGS>
                static auto message_to_tuple(Message &message) {
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
    std::enable_if_t<(sizeof...(ARGS) > 1), std::tuple<ARGS...>>
    force_unpack(Message& message) {
        return gadgetron_detail::detail::message_to_tuple<ARGS...>(message);
    }


    template<class T>
    T force_unpack(Message& message) {
        return gadgetron_detail::detail::message_to_tuple<T>(message);
    }

    template<class ...ARGS>
    std::enable_if_t<(sizeof...(ARGS) > 1), std::tuple<ARGS...>>
    unpack(Message& message) {
        if (convertible_to<ARGS...>(message)) {
            return force_unpack<ARGS...>(message);
        }
        throw std::runtime_error("Unable to unpack messages");
    }

    template<class T>
    T unpack(Message& message) {
        if (convertible_to<T>(message)) {
            return force_unpack<T>(message);
        }
        throw std::runtime_error("Unable to unpack message");

    }

}
