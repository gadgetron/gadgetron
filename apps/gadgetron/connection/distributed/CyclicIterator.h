#pragma once

template<class Iter>
class CyclicIterator {
public:


    CyclicIterator(const Iter &start, const Iter &end) : start(start), current(start), end(end) {};

    auto operator*() {
        return *current;
    }


    bool operator==(const CyclicIterator &other) {
        return current == other.current;
    }

    CyclicIterator &operator++() {
        current++;
        if (current == end) current = start;

        return *this;
    }

private:
    const Iter start;
    Iter current;
    const Iter end;

};


template<class Iter>
CyclicIterator<Iter> make_cyclic(const Iter &begin, const Iter &end) {
    return CyclicIterator<Iter>(begin, end);
};

