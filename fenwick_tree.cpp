template <class T>
struct fenwick_tree : std::vector<T> {
    using base = std::vector<T>;

    class proxy {
        fenwick_tree& tree;
        std::size_t i;

    public:
        proxy(fenwick_tree& t, std::size_t idx) : tree(t), i(idx) { tree.check(i, i); }
        proxy(const proxy&) = default;

        operator T() const { return tree.sum(i, i); }
        proxy& operator+=(const T& k) { tree.add(i, k); return *this; }
        proxy& operator-=(const T& k) { tree.add(i, -k); return *this; }
        proxy& operator++() { return *this += T(1); }
        proxy& operator--() { return *this -= T(1); }
        proxy& operator=(const proxy& o) { return *this = T(o); }
        proxy& operator=(const T& v) { return *this += v - T(*this); }

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

    fenwick_tree() = default;
    explicit fenwick_tree(std::size_t n) : base(n) {}
    fenwick_tree(std::initializer_list<T> a) : fenwick_tree(std::span(a)) {}

    template <std::ranges::forward_range R>
        requires(!requires { typename std::remove_cvref_t<R>::proxy; } &&
                 std::convertible_to<std::ranges::range_reference_t<R>, T>)
    explicit fenwick_tree(R&& a) : base(static_cast<std::size_t>(std::ranges::distance(a))) {
        auto& t = static_cast<base&>(*this);
        for (std::size_t i = 0; const T x : a) {
            t[i] += x;
            if (std::size_t j = i | (i + 1); j < t.size()) t[j] += t[i];
            ++i;
        }
    }

    proxy operator[](std::size_t i) { return {*this, i}; }
    T operator[](std::size_t i) const { return (*this)[i, i]; }
    T operator[](std::size_t l, std::size_t r) const { check(l, r); return sum(l, r); }

private:
    void check(std::size_t l, std::size_t r) const {
        if (l > r) throw std::invalid_argument("fenwick_tree: l > r");
        if (r >= this->size()) throw std::out_of_range("fenwick_tree: r >= size()");
    }

    T prefix(std::size_t p) const {
        const auto& t = static_cast<const base&>(*this);
        T s{};
        for (; p > 0; p &= p - 1) s += t[p - 1];
        return s;
    }

    void add(std::size_t p, const T& v) {
        auto& t = static_cast<base&>(*this);
        for (; p < t.size(); p |= p + 1) t[p] += v;
    }

    T sum(std::size_t l, std::size_t r) const { return prefix(r + 1) - prefix(l); }
};
