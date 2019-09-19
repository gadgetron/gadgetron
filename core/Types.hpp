#pragma once

namespace Gadgetron { namespace Core {

    // Adapted from cppreference.com

    namespace { namespace gadgetron_detail {
        template <class F, class Tuple, std::size_t... I>
        inline constexpr decltype(auto) apply_impl(F&& f, Tuple&& t, std::index_sequence<I...>) {
            return f(get<I>(std::forward<Tuple>(t))...);
        }
    }}

    template <class F, class TupleLike> constexpr decltype(auto) apply(F&& f, TupleLike&& t) {
        return gadgetron_detail::apply_impl(std::forward<F>(f), std::forward<TupleLike>(t),
            std::make_index_sequence<std::tuple_size<std::remove_reference_t<TupleLike>>::value>{});
    }

}}
