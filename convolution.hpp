#ifndef CONVOLUTION_HPP
#define CONVOLUTION_HPP

#include <algorithm>
#include <bit>
#include <cmath>
#include <complex>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numbers>
#include <ranges>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

namespace convolution_detail {

template <typename T>
concept Coefficient = (std::integral<T> || std::floating_point<T>)
                   && !std::same_as<T, bool>;

using Real = long double;
using Complex = std::complex<Real>;

template <typename T>
void bit_reverse(std::vector<T>& values) {
    const std::size_t count = values.size();
    std::size_t reversed_index = 0;
    for (std::size_t index = 1; index < count; ++index) {
        std::size_t bit = count >> 1;
        while ((reversed_index & bit) != 0) {
            reversed_index ^= bit;
            bit >>= 1;
        }
        reversed_index ^= bit;
        if (index < reversed_index) {
            std::swap(values[index], values[reversed_index]);
        }
    }
}

inline std::vector<Complex> make_roots(std::size_t count) {
    std::vector<Complex> roots(count / 2);
    const Real angle = 2 * std::numbers::pi_v<Real> / static_cast<Real>(count);
    for (std::size_t index = 0; index < roots.size(); ++index) {
        roots[index] = std::polar(Real{1}, angle * static_cast<Real>(index));
    }
    if (count >= 4) {
        roots[count / 4] = Complex{0, 1};
    }
    return roots;
}

inline void fft(std::vector<Complex>& values, const std::vector<Complex>& roots, bool inverse) {
    const std::size_t count = values.size();
    bit_reverse(values);
    for (std::size_t block_size = 2; block_size <= count; block_size <<= 1) {
        const std::size_t half_block = block_size / 2;
        const std::size_t root_stride = count / block_size;
        for (std::size_t block_start = 0; block_start < count; block_start += block_size) {
            for (std::size_t offset = 0; offset < half_block; ++offset) {
                const std::size_t left_index = block_start + offset;
                const std::size_t right_index = left_index + half_block;
                const Complex root = inverse ? std::conj(roots[offset * root_stride])
                                             : roots[offset * root_stride];
                const Complex left = values[left_index];
                const Complex right = values[right_index] * root;

                values[left_index] = left + right;
                values[right_index] = left - right;
            }
        }
    }
    if (inverse) {
        const Real scale = static_cast<Real>(count);
        for (Complex& value : values) {
            value /= scale;
        }
    }
}

consteval bool is_prime(std::uint32_t value) {
    if (value < 2 || value % 2 == 0) {
        return value == 2;
    }
    for (std::uint32_t divisor = 3; divisor <= value / divisor; divisor += 2) {
        if (value % divisor == 0) {
            return false;
        }
    }
    return true;
}

constexpr std::uint32_t power(std::uint32_t base, std::uint32_t exponent,
                              std::uint32_t modulus) {
    std::uint64_t result = 1;
    std::uint64_t factor = base;
    while (exponent != 0) {
        if ((exponent & 1) != 0) {
            result = result * factor % modulus;
        }
        factor = factor * factor % modulus;
        exponent >>= 1;
    }
    return static_cast<std::uint32_t>(result);
}

template <std::uint32_t MOD>
struct NTT {
    static_assert(is_prime(MOD), "MOD must be prime");

    static constexpr std::uint32_t order = MOD > 1
        ? std::uint32_t{1} << std::countr_zero(MOD - 1) : 1;
    static constexpr std::uint32_t root = [] consteval {
        if constexpr (MOD <= 2) {
            return std::uint32_t{1};
        }
        for (std::uint32_t candidate = 2; candidate < MOD; ++candidate) {
            const std::uint32_t value = power(candidate, (MOD - 1) / order, MOD);
            if (power(value, order / 2, MOD) != 1) {
                return value;
            }
        }
        return std::uint32_t{0};
    }();
    static_assert(MOD < 2 || (power(root, order, MOD) == 1
                  && (order == 1 || power(root, order / 2, MOD) != 1)),
                  "Invalid NTT root");

    static void transform(std::vector<std::uint32_t>& values, bool inverse) {
        const std::size_t count = values.size();
        bit_reverse(values);

        const std::uint32_t transform_root = inverse ? power(root, MOD - 2, MOD) : root;
        for (std::size_t block_size = 2; block_size <= count; block_size <<= 1) {
            const std::size_t half_block = block_size / 2;
            const std::uint32_t block_root = power(
                transform_root, order / static_cast<std::uint32_t>(block_size), MOD);

            for (std::size_t block_start = 0; block_start < count; block_start += block_size) {
                std::uint64_t root_power = 1;
                for (std::size_t offset = 0; offset < half_block; ++offset) {
                    const std::size_t left_index = block_start + offset;
                    const std::size_t right_index = left_index + half_block;
                    const std::uint64_t left = values[left_index];
                    const std::uint64_t right = values[right_index] * root_power % MOD;
                    const std::uint64_t sum = left + right;

                    values[left_index] = static_cast<std::uint32_t>(sum >= MOD ? sum - MOD : sum);
                    values[right_index] = static_cast<std::uint32_t>(
                        left >= right ? left - right : left + MOD - right);
                    root_power = root_power * block_root % MOD;
                }
            }
        }

        if (inverse) {
            const std::uint64_t scale = power(static_cast<std::uint32_t>(count), MOD - 2, MOD);
            for (std::uint32_t& value : values) {
                value = static_cast<std::uint32_t>(value * scale % MOD);
            }
        }
    }
};

}

template <convolution_detail::Coefficient T>
[[nodiscard]] std::vector<T> convolution(const std::vector<T>& first, const std::vector<T>& second) {
    using convolution_detail::Real;
    using convolution_detail::Complex;

    if (first.empty() || second.empty()) {
        return {};
    }

    const std::size_t max_transform_size = std::bit_floor(std::vector<Complex>{}.max_size());
    if (first.size() > max_transform_size
        || second.size() - 1 > max_transform_size - first.size()) {
        throw std::length_error("Convolution exceeds the supported size");
    }
    const std::size_t result_size = first.size() + second.size() - 1;
    const std::size_t transform_size = std::bit_ceil(result_size);
    std::vector<Complex> first_transform(transform_size);
    std::vector<Complex> second_transform(transform_size);
    const auto to_complex = [](T value) {
        const Real converted = static_cast<Real>(value);
        if (!std::isfinite(converted)) {
            throw std::invalid_argument("Convolution inputs must be finite");
        }
        return Complex{converted};
    };
    std::ranges::transform(first, first_transform.begin(), to_complex);
    std::ranges::transform(second, second_transform.begin(), to_complex);
    const std::vector<Complex> roots = convolution_detail::make_roots(transform_size);

    convolution_detail::fft(first_transform, roots, false);
    convolution_detail::fft(second_transform, roots, false);
    for (auto&& [left, right] : std::views::zip(first_transform, second_transform)) {
        left *= right;
    }
    convolution_detail::fft(first_transform, roots, true);

    std::vector<T> result(result_size);
    for (auto&& [transformed, coefficient] : std::views::zip(first_transform, result)) {
        Real value = transformed.real();
        if (!std::isfinite(value) || !std::isfinite(transformed.imag())) {
            throw std::overflow_error("FFT produced a non-finite coefficient");
        }
        if constexpr (std::integral<T>) {
            value = std::round(value);
            const Real upper_bound = std::ldexp(Real{1}, std::numeric_limits<T>::digits);
            const Real lower_bound = std::is_signed_v<T> ? -upper_bound : Real{0};
            if (value < lower_bound || value >= upper_bound) {
                throw std::overflow_error("Convolution coefficient exceeds the result type");
            }
        } else {
            const Real maximum = static_cast<Real>(std::numeric_limits<T>::max());
            if (value < -maximum || value > maximum) {
                throw std::overflow_error("Convolution coefficient exceeds the result type");
            }
        }
        coefficient = static_cast<T>(value);
    }
    return result;
}

template <std::uint32_t MOD, convolution_detail::Coefficient T>
    requires std::integral<T>
[[nodiscard]] auto convolution(const std::vector<T>& first, const std::vector<T>& second) {
    using Transform = convolution_detail::NTT<MOD>;
    using Result = std::conditional_t<
        (static_cast<std::uintmax_t>(std::numeric_limits<T>::max()) >= MOD - 1),
        T, std::uint32_t>;

    if (first.empty() || second.empty()) {
        return std::vector<Result>{};
    }

    const std::size_t max_transform_size = std::min(
        static_cast<std::size_t>(Transform::order),
        std::bit_floor(std::vector<std::uint32_t>{}.max_size()));
    if (first.size() > max_transform_size
        || second.size() - 1 > max_transform_size - first.size()) {
        throw std::length_error("Convolution length is not supported by MOD");
    }
    const std::size_t result_size = first.size() + second.size() - 1;
    const std::size_t transform_size = std::bit_ceil(result_size);

    std::vector<std::uint32_t> first_transform(transform_size);
    std::vector<std::uint32_t> second_transform(transform_size);
    const auto normalize = [](T value) {
        if constexpr (std::signed_integral<T>) {
            std::intmax_t reduced = static_cast<std::intmax_t>(value) % static_cast<std::intmax_t>(MOD);
            if (reduced < 0) {
                reduced += MOD;
            }
            return static_cast<std::uint32_t>(reduced);
        } else {
            return static_cast<std::uint32_t>(static_cast<std::uintmax_t>(value) % MOD);
        }
    };
    std::ranges::transform(first, first_transform.begin(), normalize);
    std::ranges::transform(second, second_transform.begin(), normalize);

    Transform::transform(first_transform, false);
    Transform::transform(second_transform, false);
    for (auto&& [left, right] : std::views::zip(first_transform, second_transform)) {
        left = static_cast<std::uint32_t>(static_cast<std::uint64_t>(left) * right % MOD);
    }
    Transform::transform(first_transform, true);

    first_transform.resize(result_size);
    return std::vector<Result>(std::from_range, first_transform);
}

#endif
