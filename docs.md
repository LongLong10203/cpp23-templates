# Quick reference

Include the header you need, for example `#include "fenwick_tree.hpp"`.
Compile from the repository root:

```sh
clang++ -std=c++23 -O2 -I. main.cpp -o main
```

All indices are zero-based. Ranges `[l, r]` include both endpoints.
The snippets below go inside your program.

## disjoint_set_union.hpp

```cpp
disjoint_set_union<> dsu(5); // Vertices 0 to 4.
dsu.join(0, 1);
bool connected = dsu.same(0, 1); // true
int count = dsu.size(0);        // 2
auto groups = dsu.groups();
```

## prefix_sum.hpp

```cpp
prefix_sum<long long> sums{2, 4, 6};
long long total = sums[1, 2]; // 10
long long prefix = sums[2];  // 12
```

For fixed data; rebuild when values change.

## fenwick_tree.hpp

```cpp
fenwick_tree<long long> tree{1, 2, 3};
tree[1] += 5;
tree[0] = 4;
long long total = tree[0, 2]; // 14
```

## range_fenwick_tree.hpp

```cpp
range_fenwick_tree<long long> tree{1, 2, 3};
tree[0, 2] += 5;             // Add to every element in the range.
tree[0] = 4;                // Assignment is for single elements.
long long total = tree[0, 2]; // 19
```

For both Fenwick trees, read sums into `T` (such as `long long`), not `auto`.

## modint.hpp

```cpp
mint a = 2;                 // Modulus 998244353.
mint b = a.pow(10);          // 1024
mint c = b / a;              // 512
auto value = c.val();
```

Use `modint<M>` for another modulus. For a runtime modulus, call
`modint<0>::set_mod(M)` before use. Division requires an invertible divisor.

## convolution.hpp

```cpp
std::vector<int> a{1, 2}, b{3, 4};
auto result = convolution(a, b);            // {3, 10, 8}
auto modular = convolution<998244353>(a, b); // Same result modulo 998244353.
```

Ordinary FFT can lose precision for large coefficients. The modular overload
requires a prime modulus; the padded power-of-two result length must divide `MOD - 1`.
