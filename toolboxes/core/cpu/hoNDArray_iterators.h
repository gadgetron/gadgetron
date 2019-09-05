#pragma once

#include <numeric>
#include <boost/iterator/iterator_facade.hpp>

namespace Gadgetron {

    template<class T, class Value>
    class hoNDArraySpanIterator : public boost::iterator_facade<hoNDArraySpanIterator<T, Value>, Value, boost::random_access_traversal_tag, Value> {
    public:
        Value dereference() const;
        void increment();
        void decrement();
        void advance(ptrdiff_t n);
        bool equal(const hoNDArraySpanIterator<T, Value> &other) const;
        ptrdiff_t distance_to(const hoNDArraySpanIterator<T, Value> &other) const;

        Value &source;
        size_t index = 0;

        const std::vector<size_t> span_dimensions;
        const size_t span_size, number_of_spans;

        hoNDArraySpanIterator(Value &source, size_t number_of_span_dimensions);

    private:
        std::vector<size_t> make_span_dimensions(size_t number_of_span_dimensions);
        size_t make_span_size();
        size_t make_number_of_spans();
    };


    template<class T, class Value>
    struct span_range_ {
        Value &source;
        const size_t number_of_span_dimensions;

        hoNDArraySpanIterator<T, Value> begin() {
            return hoNDArraySpanIterator<T, Value>{source, number_of_span_dimensions};
        }

        hoNDArraySpanIterator<T, Value> end() {
            hoNDArraySpanIterator<T, Value> iterator{source, number_of_span_dimensions};
            iterator.index = iterator.number_of_spans;
            return iterator;
        }
    };

    template<class T, class Value>
    auto begin(span_range_<T, Value> range) { return range.begin(); }

    template<class T, class Value>
    auto end(span_range_<T, Value> range) { return range.end(); }

    template<class T>
    span_range_<T, hoNDArray<T>> spans(hoNDArray<T> &array, size_t number_of_span_dimensions) {
        return span_range_<T, hoNDArray<T>>{array, number_of_span_dimensions};
    }

    template<class T>
    span_range_<T, const hoNDArray<T>> spans(const hoNDArray<T> &array, size_t number_of_span_dimensions) {
        return span_range_<T, const hoNDArray<T>>{array, number_of_span_dimensions};
    }

    template<class T, class Value>
    hoNDArraySpanIterator<T, Value>::hoNDArraySpanIterator(
            Value &source,
            size_t number_of_span_dimensions
    ) : source(source),
        span_dimensions(make_span_dimensions(number_of_span_dimensions)),
        span_size(make_span_size()),
        number_of_spans(make_number_of_spans()) {}

    template<class T, class Value>
    std::vector<size_t> hoNDArraySpanIterator<T, Value>::make_span_dimensions(size_t number_of_span_dimensions) {
        std::vector<size_t> dimensions(number_of_span_dimensions);
        std::copy_n(source.dimensions().begin(), number_of_span_dimensions, dimensions.begin());
        return dimensions;
    }

    template<class T, class Value>
    size_t hoNDArraySpanIterator<T, Value>::make_span_size() {
        return std::accumulate(span_dimensions.begin(), span_dimensions.end(), size_t(1), std::multiplies<>());
    }

    template<class T, class Value>
    size_t hoNDArraySpanIterator<T, Value>::make_number_of_spans() {
        return source.size() / span_size;
    }

    template<class T, class Value>
    Value hoNDArraySpanIterator<T, Value>::dereference() const {
        return hoNDArray<T>(span_dimensions, const_cast<T *>(source.data()) + index * span_size);
    }

    template<class T, class Value>
    void hoNDArraySpanIterator<T, Value>::increment() { index++; }

    template<class T, class Value>
    void hoNDArraySpanIterator<T, Value>::decrement() { index--; }

    template<class T, class Value>
    void hoNDArraySpanIterator<T, Value>::advance(ptrdiff_t n) { index += n; }

    template<class T, class Value>
    bool hoNDArraySpanIterator<T, Value>::equal(const hoNDArraySpanIterator<T, Value> &other) const {
        return (&source == &other.source) && (index == other.index);
    }

    template<class T, class Value>
    ptrdiff_t hoNDArraySpanIterator<T, Value>::distance_to(const hoNDArraySpanIterator<T, Value> &other) const {
        return ptrdiff_t(other.index) - ptrdiff_t(index);
    }
}