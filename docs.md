# Usage guide

## Getting started

Include the headers you need directly; there is no separate library to build or link. Add the repository root to your compiler's include path:

```sh
clang++ -std=c++23 -O2 -I /path/to/cpp23-templates main.cpp -o main
```

Verified with Apple Clang 21.0.0 and the bundled libc++ (`_LIBCPP_VERSION=220106`) on Apple Silicon, using `-std=c++23`. No extra compiler installation is needed for that toolchain.

The disjoint-set union uses an ordinary indexed loop for grouping, so it does not require `std::views::enumerate`. Other C++23 features used by the templates, including multidimensional subscripting, `std::ranges::iota`, `std::views::zip`, and `std::from_range`, are supported by the verified toolchain. Older compiler or standard-library versions may still need an upgrade.

All templates are in the global namespace. The examples below are standalone programs.

| Header | Purpose |
| --- | --- |
| [disjoint_set_union.hpp](disjoint_set_union.hpp) | Maintain connected components with union and find. |
| [prefix_sum.hpp](prefix_sum.hpp) | Answer range sums over fixed data. |
| [fenwick_tree.hpp](fenwick_tree.hpp) | Point updates and range sums. |
| [range_fenwick_tree.hpp](range_fenwick_tree.hpp) | Range additions and range sums. |
| [modint.hpp](modint.hpp) | Arithmetic under a fixed or runtime modulus. |
| [convolution.hpp](convolution.hpp) | Polynomial multiplication using FFT or modular NTT. |

## Indexing and storage

Prefix sums and both Fenwick trees use zero-based indices and **inclusive** ranges: `a[l, r]` includes both endpoints. Range accesses throw `std::invalid_argument` when `l > r` and `std::out_of_range` when `r >= size()`.

These containers inherit from `std::vector`, but their backing storage represents computed data. Avoid modifying it through inherited operations such as `push_back`, `resize`, or iterators; rebuild the container or use its documented updates instead.

Fenwick updates return a proxy object. Read into a variable of the element type, such as `long long sum = tree[0, 2]`, to take a value snapshot; `auto` can retain a proxy instead.

## Disjoint-set union

```cpp
#include "disjoint_set_union.hpp"
#include <cassert>

int main() {
    disjoint_set_union<> dsu(5); // Vertices 0 through 4.
    bool merged = dsu.join(0, 1);
    assert(merged);
    assert(dsu.same(0, 1));
    assert(dsu.size(0) == 2);
    assert(dsu.comps == 4);
    auto groups = dsu.groups();
    assert(groups.size() == 4);
}
```

- `join(a, b)` merges components and returns whether a merge occurred.
- `find(x)` returns a representative label; `same(a, b)` checks connectivity.
- `size(x)` returns the component size; `size()` returns the total number of registered labels.
- `comps` stores the component count; `groups()` returns the labels in each component.
- An empty `disjoint_set_union<T>` registers labels on demand, including during queries. Labels must support equality and hashing; integral labels, hashable types such as `std::string`, and tuples/pairs of supported labels work. Include the headers for your label types.

The size constructor is available for signed integral labels. Union by size and path compression give near-constant amortized union/find operations, assuming constant-time hashing.

## Prefix sums

```cpp
#include "prefix_sum.hpp"
#include <cassert>

int main() {
    prefix_sum<long long> sums{2, 4, 6, 8};
    assert((sums[1, 3] == 18));
    assert(sums[2] == 12); // Prefix total, not the original element.
}
```

Construct from an initializer list or input range. Construction takes O(n), and range queries take O(1), with O(n) storage. Rebuild after changing the original data. Single-index access reads the stored prefix total and follows `std::vector::operator[]` bounds rules.

## Fenwick tree

```cpp
#include "fenwick_tree.hpp"
#include <cassert>

int main() {
    fenwick_tree<long long> tree{1, 2, 3, 4};
    tree[1] += 5;
    tree[3] = 10;
    long long sum = tree[1, 3];
    assert(sum == 20);
}
```

Construct from an initializer list, a forward range, or a size: `fenwick_tree<long long> tree(10)` creates ten zeros. `tree[i]` reads or updates an element; `tree[l, r]` reads a range sum. Point updates support assignment, `+=`, `-=`, `++`, and `--`.

Construction takes O(n), reads and updates take O(log n), and storage is O(n).

## Range Fenwick tree

```cpp
#include "range_fenwick_tree.hpp"
#include <cassert>

int main() {
    range_fenwick_tree<long long> tree{1, 2, 3, 4};
    tree[1, 3] += 5; // Add 5 to each element in the range.
    tree[0] = 10;
    long long sum = tree[0, 3];
    assert(sum == 34);
}
```

Construction options match `fenwick_tree`. Use `tree[l, r] += value` or `-= value` to update every element in a range. Reading `tree[l, r]` returns the sum. Assignment requires a single element (`tree[i]` or `tree[i, i]`); assigning to a larger range throws `std::invalid_argument`.

Construction takes O(n), reads and updates take O(log n), and storage is O(n).

## Modular arithmetic

```cpp
#include "modint.hpp"
#include <cassert>

int main() {
    mint a = 2; // Alias for modint<998244353>.
    assert(a.pow(10).val() == 1024);
    assert((a / a).val() == 1);

    using dynamic_mint = modint<0>;
    dynamic_mint::set_mod(17);
    dynamic_mint b = -3;
    assert(b.val() == 14);
}
```

- `modint<M>` uses a fixed modulus for `1 <= M <= 2^31`; `modint<0>` uses a shared runtime modulus in the same range.
- Integral inputs are normalized. Arithmetic, comparisons for equality, increment/decrement, and stream input/output are supported.
- `val()` returns the stored residue; `mod()` returns the modulus.
- `pow(e)` supports integral exponents, including negative exponents when an inverse exists.
- `inv()` and division require the value being inverted to be coprime to the modulus.
- Set the runtime modulus before use, and do not change it while continuing to use existing values.
- `raw(x)` skips normalization and requires `x < mod()`.

Addition, subtraction, and multiplication take O(1). Inversion takes O(log M), and nonnegative exponentiation takes O(log e).

## Convolution

```cpp
#include "convolution.hpp"
#include <cassert>
#include <vector>

int main() {
    std::vector<int> a{1, 2, 3}, b{4, 5};
    const std::vector<int> expected{4, 13, 22, 15};
    auto ordinary = convolution(a, b);
    auto modular = convolution<998244353>(a, b);
    assert(ordinary == expected);
    assert(modular == expected);
}
```

The input vectors must have the same arithmetic element type; `bool` is excluded. Coefficients are ordered from the constant term upward. Nonempty inputs of lengths n and m produce n + m - 1 coefficients; an empty input produces an empty result.

- `convolution(a, b)` uses a floating-point FFT and returns the input element type. Integral outputs are rounded, so large coefficients can lose precision. Inputs must be finite; detected non-finite or out-of-range results throw exceptions.
- `convolution<MOD>(a, b)` requires integral inputs and a compile-time prime modulus. It computes exact coefficients modulo `MOD`. The transform length (the smallest power of two at least n + m - 1) must divide `MOD - 1`; unsupported lengths throw `std::length_error`. For `998244353`, the maximum transform length is `2^23`.
- The modular overload returns the input element type if it can represent every residue, otherwise `std::uint32_t`.

Both overloads use O(k log k) arithmetic operations and O(k) storage for transform length k.
