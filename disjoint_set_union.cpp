template <std::equality_comparable T = int>
struct disjoint_set_union : std::vector<int> {
    struct hasher {
        static std::uint64_t mix(std::uint64_t x) {
            static const std::uint64_t seed =
                std::chrono::steady_clock::now().time_since_epoch().count();
            x += seed + 0x9e3779b97f4a7c15;
            x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9;
            x = (x ^ (x >> 27)) * 0x94d049bb133111eb;
            return x ^ (x >> 31);
        }

        template <class U>
        static std::size_t operator()(const U& x) {
            if constexpr (std::integral<U>) {
                return mix(x);
            } else if constexpr (requires { std::tuple_size<U>::value; }) {
                std::uint64_t h = 0;
                std::apply([&](const auto&... part) {
                    ((h = mix(h ^ hasher{}(part))), ...);
                }, x);
                return h;
            } else {
                return mix(std::hash<U>{}(x));
            }
        }
    };

    int dense = 0, comps = 0;
    std::vector<T> ids;
    std::unordered_map<T, int, hasher> pos;
    using std::vector<int>::size;

    disjoint_set_union() = default;

    disjoint_set_union(int n) requires std::signed_integral<T>
        : std::vector<int>(n, -1), dense(n), comps(n), ids(n) {
        std::ranges::iota(ids, 0);
    }

    int index(const T& x) {
        if constexpr (std::signed_integral<T>) {
            if (x >= 0 && x < dense) return x;
        }
        auto [it, added] = pos.try_emplace(x, size());
        if (added) {
            push_back(-1);
            ids.push_back(x);
            comps++;
        }
        return it->second;
    }

    int root(int v) {
        if ((*this)[v] < 0) return v;
        return (*this)[v] = root((*this)[v]);
    }

    T find(const T& x) { return ids[root(index(x))]; }
    bool same(const T& a, const T& b) { return root(index(a)) == root(index(b)); }
    int size(const T& x) { return -(*this)[root(index(x))]; }

    bool join(const T& a, const T& b) {
        int u = root(index(a));
        int v = root(index(b));
        if (u == v) return false;
        if ((*this)[u] > (*this)[v]) std::swap(u, v);
        (*this)[u] += (*this)[v];
        (*this)[v] = u;
        comps--;
        return true;
    }

    std::vector<std::vector<T>> groups() {
        std::vector<std::vector<T>> g(size());
        for (auto [v, x] : std::views::enumerate(ids)) {
            g[root(v)].push_back(x);
        }
        std::erase_if(g, std::ranges::empty);
        return g;
    }
};
