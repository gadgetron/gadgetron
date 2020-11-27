//
// Created by dchansen on 4/9/19.
//

#pragma once

namespace Gadgetron { namespace Core {

    template <class T> constexpr bool is_trivially_copyable_v = std::is_trivially_copyable<T>::value;

    template <class T, class V> constexpr bool is_same_v = std::is_same<T, V>::value;

    template <class T, class V> constexpr bool is_convertible_v = std::is_convertible<T, V>::value;

    template <bool... ARGS> struct all_of;

    template <> struct all_of<> : std::true_type {};

    template <bool... ARGS> struct all_of<false, ARGS...> : std::false_type {};

    template <bool... ARGS> struct all_of<true, ARGS...> : all_of<ARGS...> {};


    template <bool... ARGS> constexpr bool all_of_v = all_of<ARGS...>::value;


    template <bool... ARGS> struct any_of;
    template <> struct any_of<> : std::false_type{};

    template <bool... ARGS> struct any_of<true, ARGS...> : std::true_type {};

    template <bool... ARGS> struct any_of<false, ARGS...> : any_of<ARGS...> {};

    template <bool... ARGS> constexpr bool any_of_v = any_of<ARGS...>::value;

    template<class T>
    constexpr bool is_floating_point_v = std::is_floating_point<T>::value;
}}
