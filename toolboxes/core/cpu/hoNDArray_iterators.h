#pragma once

#include <numeric>
#include <boost/iterator/iterator_facade.hpp>

namespace Gadgetron {

    template<class T>
    class hoNDArraySpanIterator : public boost::iterator_facade<hoNDArraySpanIterator<T>, hoNDArray<T>, boost::random_access_traversal_tag, hoNDArray<T>> {
    public:

        hoNDArray<T> dereference() const;
        void increment();
        void decrement();
        void advance(ptrdiff_t n);
        bool equal(const hoNDArraySpanIterator<T> &other) const;
        ptrdiff_t distance_to(const hoNDArraySpanIterator<T> &other) const;

        hoNDArray<T> &source;
        size_t index = 0;

        const std::vector<size_t> span_dimensions;
        const size_t span_size, number_of_spans;

        hoNDArraySpanIterator(hoNDArray<T> &source, size_t number_of_span_dimensions);

    private:
        std::vector<size_t> make_span_dimensions(size_t number_of_span_dimensions);
        size_t make_span_size();
        size_t make_number_of_spans();
    };

    template<class T>
    struct span_range_ {
        hoNDArray<T> &source;
        const size_t number_of_span_dimensions;

        hoNDArraySpanIterator<T> begin() {
            return hoNDArraySpanIterator<T>{source, number_of_span_dimensions};
        }

        hoNDArraySpanIterator<T> end() {
            hoNDArraySpanIterator<T> iterator{source, number_of_span_dimensions};
            iterator.index = iterator.number_of_spans;
            return iterator;
        }
    };

    template<class T>
    auto begin(span_range_<T> range) { return range.begin(); }

    template<class T>
    auto end(span_range_<T> range) { return range.end(); }

    template<class T>
    span_range_<T> spans(hoNDArray<T> &array, size_t number_of_span_dimensions) {
        return span_range_<T>{array, number_of_span_dimensions};
    }

    template<class T>
    hoNDArraySpanIterator<T>::hoNDArraySpanIterator(
            hoNDArray<T> &source,
            size_t number_of_span_dimensions
    ) : source(source),
        span_dimensions(make_span_dimensions(number_of_span_dimensions)),
        span_size(make_span_size()),
        number_of_spans(make_number_of_spans()) {}

    template<class T>
    std::vector<size_t> hoNDArraySpanIterator<T>::make_span_dimensions(size_t number_of_span_dimensions) {
        std::vector<size_t> dimensions(number_of_span_dimensions);
        std::copy_n(source.dimensions().begin(), number_of_span_dimensions, dimensions.begin());
        return dimensions;
    }

    template<class T>
    size_t hoNDArraySpanIterator<T>::make_span_size() {
        return std::accumulate(span_dimensions.begin(), span_dimensions.end(), size_t(1), std::multiplies<>());
    }

    template<class T>
    size_t hoNDArraySpanIterator<T>::make_number_of_spans() {
        return source.size() / span_size;
    }

template<class T>
    hoNDArray<T> hoNDArraySpanIterator<T>::dereference() const {
        return hoNDArray<T>(span_dimensions, source.data() + index * span_size);
    }

    template<class T>
    void hoNDArraySpanIterator<T>::increment() { index++; }

    template<class T>
    void hoNDArraySpanIterator<T>::decrement() { index--; }

    template<class T>
    void hoNDArraySpanIterator<T>::advance(ptrdiff_t n) { index += n; }

    template<class T>
    bool hoNDArraySpanIterator<T>::equal(const hoNDArraySpanIterator<T> &other) const {
        return (source == other.source) && (index == other.index);
    }

    template<class T>
    ptrdiff_t hoNDArraySpanIterator<T>::distance_to(const hoNDArraySpanIterator<T> &other) const {
        return ptrdiff_t(other.index) - ptrdiff_t(index);
    }
}