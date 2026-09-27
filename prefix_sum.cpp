template <class T>
struct prefix_sum : std::vector<T> {
    using std::vector<T>::operator[];

    prefix_sum() = default;

    template <std::ranges::input_range R>
        requires (!std::same_as<std::remove_cvref_t<R>, prefix_sum>)
    explicit prefix_sum(R&& a) {
        if constexpr (std::ranges::sized_range<R>)
            this->reserve(std::ranges::size(a));
        for (T s{}; auto&& x : a)
            this->push_back(s += x);
    }

    prefix_sum(std::initializer_list<T> a) : prefix_sum(std::views::all(a)) {}

    T operator[](std::size_t l, std::size_t r) const {
        if (l > r) throw std::invalid_argument("prefix_sum: l > r");
        if (r >= this->size()) throw std::out_of_range("prefix_sum: r >= size()");
        return (*this)[r] - (l ? (*this)[l - 1] : T{});
    }
};

template <std::ranges::input_range R>
prefix_sum(R&&) -> prefix_sum<std::ranges::range_value_t<R>>;
