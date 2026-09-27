template <class T>
struct range_fenwick_tree : std::vector<std::array<T, 2>> {
    using base = std::vector<std::array<T, 2>>;

    class proxy {
        range_fenwick_tree& tree;
        std::size_t l, r;

    public:
        proxy(range_fenwick_tree& t, std::size_t lo, std::size_t hi) : tree(t), l(lo), r(hi) {
            tree.check(l, r);
        }
        proxy(const proxy&) = default;

        operator T() const { return tree.sum(l, r); }
        proxy& operator+=(const T& k) { tree.add(l, r, k); return *this; }
        proxy& operator-=(const T& k) { tree.add(l, r, -k); return *this; }
        proxy& operator++() { return *this += T(1); }
        proxy& operator--() { return *this -= T(1); }

        proxy& operator=(const proxy& o) { return *this = T(o); }
        proxy& operator=(const T& v) {
            if (l != r) throw std::invalid_argument("range_fenwick_tree: assignment needs l == r");
            return *this += v - T(*this);
        }

        T operator++(int) {
            T old = *this;
            ++*this;
            return old;
        }
        T operator--(int) {
            T old = *this;
            --*this;
            return old;
        }
    };

    range_fenwick_tree() = default;
    explicit range_fenwick_tree(std::size_t n) : base(n) {}
    range_fenwick_tree(std::initializer_list<T> a) : range_fenwick_tree(std::span(a)) {}

    template <std::ranges::forward_range R>
        requires std::convertible_to<std::ranges::range_reference_t<R>, T>
    explicit range_fenwick_tree(R&& a) : base(static_cast<std::size_t>(std::ranges::distance(a))) {
        auto& t = static_cast<base&>(*this);
        T prev{};
        for (std::size_t i = 0; const T x : a) {
            const T d = x - prev;
            prev = x;
            t[i][0] += d;
            t[i][1] += d * static_cast<T>(i);
            if (std::size_t j = i | (i + 1); j < t.size()) {
                t[j][0] += t[i][0];
                t[j][1] += t[i][1];
            }
            ++i;
        }
    }

    proxy operator[](std::size_t l, std::size_t r) { return {*this, l, r}; }
    proxy operator[](std::size_t i) { return {*this, i, i}; }
    T operator[](std::size_t l, std::size_t r) const { check(l, r); return sum(l, r); }
    T operator[](std::size_t i) const { return (*this)[i, i]; }

private:
    void check(std::size_t l, std::size_t r) const {
        if (l > r) throw std::invalid_argument("range_fenwick_tree: l > r");
        if (r >= this->size()) throw std::out_of_range("range_fenwick_tree: r >= size()");
    }

    T prefix(std::size_t p) const {
        const auto& t = static_cast<const base&>(*this);
        T s0{}, s1{};
        for (std::size_t i = p; i > 0; i &= i - 1) {
            s0 += t[i - 1][0];
            s1 += t[i - 1][1];
        }
        return static_cast<T>(p) * s0 - s1;
    }

    void update(std::size_t p, const T& v) {
        auto& t = static_cast<base&>(*this);
        for (std::size_t i = p; i < t.size(); i |= i + 1) {
            t[i][0] += v;
            t[i][1] += v * static_cast<T>(p);
        }
    }

    T sum(std::size_t l, std::size_t r) const { return prefix(r + 1) - prefix(l); }
    void add(std::size_t l, std::size_t r, const T& k) { update(l, k); update(r + 1, -k); }
};
