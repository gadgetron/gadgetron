#pragma once

template<class Iter>
class CyclicIterator {
public:
    CyclicIterator(const Iter &start, const Iter &end) : start(start), current(start), end(end) {};

    auto operator*() { return *current; }

    bool operator==(const CyclicIterator &other) { return current == other.current; }

    CyclicIterator &operator++() {
        current++;
        if (current == end) current = start;

        return *this;
    }

private:
    Iter current;
    const Iter start, end;
};

template<class Iter>
CyclicIterator<Iter> make_cyclic(const Iter &begin, const Iter &end) {
    if (begin == end) throw std::runtime_error("Empty range provided");
    return CyclicIterator<Iter>(begin, end);
};
