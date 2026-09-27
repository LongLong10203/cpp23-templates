template <std::signed_integral T = int>
struct disjoint_set_union : std::vector<int> {
    struct hasher {
        static std::size_t operator()(std::uint64_t x) {
            static const std::uint64_t seed =
                std::chrono::steady_clock::now().time_since_epoch().count();
            x += seed + 0x9e3779b97f4a7c15;
            x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9;
            x = (x ^ (x >> 27)) * 0x94d049bb133111eb;
            return x ^ (x >> 31);
        }
    };

    int dense, comps;
    std::vector<T> ids;
    std::unordered_map<T, int, hasher> pos;
    using std::vector<int>::size;

    disjoint_set_union(int n = 0) : std::vector<int>(n, -1), dense(n), comps(n), ids(n) {
        std::ranges::iota(ids, 0);
    }

    int index(T x) {
        if (x >= 0 && x < dense) return x;
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

    T find(T x) { return ids[root(index(x))]; }
    bool same(T a, T b) { return root(index(a)) == root(index(b)); }
    int size(T x) { return -(*this)[root(index(x))]; }

    bool join(T a, T b) {
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
