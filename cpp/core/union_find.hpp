#pragma once

#include <cstddef>
#include <vector>

namespace smart_core {

template <class Index>
struct UnionFind {
    std::vector<Index> parent;

    void reset(std::size_t n) {
        parent.resize(n);
        for (std::size_t i = 0; i < n; ++i) {
            parent[i] = static_cast<Index>(i);
        }
    }

    bool empty() const {
        return parent.empty();
    }

    Index find(Index x) {
        while (parent[static_cast<std::size_t>(x)] != x) {
            parent[static_cast<std::size_t>(x)] =
                parent[static_cast<std::size_t>(
                    parent[static_cast<std::size_t>(x)])];
            x = parent[static_cast<std::size_t>(x)];
        }
        return x;
    }

    Index& operator[](std::size_t i) {
        return parent[i];
    }

    const Index& operator[](std::size_t i) const {
        return parent[i];
    }

    std::size_t bytes() const {
        return parent.size() * sizeof(Index);
    }
};

} // namespace smart_core
