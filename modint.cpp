template <std::uint32_t M> requires (M <= (1u << 31))
struct modint {
    using value_type = std::uint32_t;

    static constexpr value_type mod() noexcept {
        if constexpr (M != 0) {
            return M;
        } else {
            assert(runtime_mod_ != 0);
            return runtime_mod_;
        }
    }
    static void set_mod(std::int64_t m) requires (M == 0) {
        assert(1 <= m && m <= (1LL << 31));
        runtime_mod_ = static_cast<value_type>(m);
    }

    constexpr modint() noexcept = default;
    template <std::integral T>
    constexpr modint(T x) noexcept : v_(reduce(x)) {}

    static constexpr modint raw(value_type x) noexcept {
        assert(x < mod());
        modint r;
        r.v_ = x;
        return r;
    }

    constexpr value_type val() const noexcept { return v_; }
    template <std::integral T>
    constexpr explicit operator T() const noexcept { return static_cast<T>(v_); }

    constexpr modint& operator+=(modint rhs) noexcept {
        v_ += rhs.v_;
        if (v_ >= mod()) v_ -= mod();
        return *this;
    }
    constexpr modint& operator-=(modint rhs) noexcept {
        v_ -= rhs.v_;
        if (v_ >= mod()) v_ += mod();
        return *this;
    }
    constexpr modint& operator*=(modint rhs) noexcept {
        v_ = static_cast<value_type>(std::uint64_t{v_} * rhs.v_ % mod());
        return *this;
    }
    constexpr modint& operator/=(modint rhs) noexcept { return *this *= rhs.inv(); }

    constexpr modint& operator++() noexcept { if (++v_ == mod()) v_ = 0; return *this; }
    constexpr modint& operator--() noexcept { if (v_ == 0) v_ = mod(); --v_; return *this; }
    constexpr modint operator++(int) noexcept { modint old = *this; ++*this; return old; }
    constexpr modint operator--(int) noexcept { modint old = *this; --*this; return old; }

    constexpr modint operator+() const noexcept { return *this; }
    constexpr modint operator-() const noexcept { return raw(v_ == 0 ? 0u : mod() - v_); }

    friend constexpr modint operator+(modint a, modint b) noexcept { return a += b; }
    friend constexpr modint operator-(modint a, modint b) noexcept { return a -= b; }
    friend constexpr modint operator*(modint a, modint b) noexcept { return a *= b; }
    friend constexpr modint operator/(modint a, modint b) noexcept { return a /= b; }
    friend constexpr bool operator==(modint, modint) noexcept = default;

    template <std::integral T>
    [[nodiscard]] constexpr modint pow(T e) const noexcept {
        modint base = *this, result = 1;
        auto n = static_cast<std::uint64_t>(e);
        if constexpr (std::signed_integral<T>) {
            if (e < 0) {
                base = inv();
                n = 0 - n;
            }
        }
        for (; n != 0; n >>= 1, base *= base)
            if (n & 1) result *= base;
        return result;
    }

    [[nodiscard]] constexpr modint inv() const noexcept {
        std::int64_t a = v_, b = mod(), x = 1, y = 0;
        while (b != 0) {
            const std::int64_t q = a / b;
            a = std::exchange(b, a - q * b);
            x = std::exchange(y, x - q * y);
        }
        assert(a == 1);
        return modint(x);
    }

    friend std::ostream& operator<<(std::ostream& os, modint x) { return os << x.v_; }
    friend std::istream& operator>>(std::istream& is, modint& x) {
        std::int64_t t;
        if (is >> t) x = modint(t);
        return is;
    }

private:
    value_type v_ = 0;
    static inline value_type runtime_mod_ = 0;

    template <std::integral T>
    static constexpr value_type reduce(T x) noexcept {
        if constexpr (std::signed_integral<T>) {
            using W = std::common_type_t<T, std::int64_t>;
            const W r = static_cast<W>(x) % static_cast<W>(mod());
            return static_cast<value_type>(r < 0 ? r + static_cast<W>(mod()) : r);
        } else {
            return static_cast<value_type>(x % mod());
        }
    }
};

using mint = modint<998244353>;
