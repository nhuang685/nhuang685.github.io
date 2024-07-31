---
layout: post
title: "NOI 2013 - 书法家"
date: 2024-07-30 16:35:01 -0700
categories: noi dp
math: true
# classes: wide
---

Problem Link: [LOJ2668](https://loj.ac/p/2668)

# Problem
This problem requires us to write "NOI" in an $N \times M$ grid, with
restrictions on what coverings have "NOI." in them.
We want to maximize the sum of the values of each cell covered by "NOI."

# Solution
Note that we can separate rectangles with width greater than $1$ into multiple
rectangles of width $1$ without largely changing the problem. This leads us to
use Dynamic Programming to solve the problem.

Let $dp[i][st]$ represent the maximum value when we have wrote on the
$1..i$th column and are in writing state $st$. The states are the following:

0. Left of `N`
1. 1st rectangle of `N`
2. intermediate rectangles of `N`
3. Last rectangle of `N`
4. Between `N` and `O`
5. First column of `O`
6. Intermediate columns of `O`
7. Last column of `O`
8. Between `O` and `I`
9. Left columns of `I`
10. Middle columns of `I`
11. Right columns of `I`
12. Right of `I`

For states $1, 2, 3, 5, 7, 10$, we can add two dimensions $l$ and $r$
that represent the lowest and highest row covered on column $i$.
Likewise, for states $6, 9$, and $11$, we can add two dimensions $l$ and $r$
that represent the two rows covered on column $i$. Using these
dimensions and the restrictions described in the problem, we can create
the dp transitions below. Here, $p$ refers to the value array.

$$
\begin{align*}
dp[i][0] &= 0, \\
dp[i][1][l][r] &= \max(dp[i - 1][0], dp[i - 1][1][l][r]) +
                  \sum_{j = l}^r p[i][j], \\
dp[i][2][l][r] &= \max(\max_{1 \leq l' < l} dp[i - 1][1][l'][r],
                       \max_{l \leq l' \leq r + 1, \max(l', r) \leq r' \leq n}
                          dp[i - 1][2][l'][r']) +
                  \sum_{j = l}^r p[i][j], \\
dp[i][3][l][r] &= \max(\max_{l \leq r' < r} dp[i - 1][2][l][r'],
                       dp[i - 1][3][l][r]) +
                  \sum_{j = l}^r p[i][j], \\
dp[i][4] &= \max(\max_{1 \leq l \leq r \leq n} dp[i - 1][3][l][r],
                       dp[i - 1][4]), \\
dp[i][5][l][r] &= dp[i - 1][4] + \sum_{j = l}^r p[i][j], \\
dp[i][6][l][r] &= \max(dp[i - 1][5][l][r], dp[i - 1][6][l][r]) +
                  (p[i][l] + p[i][r]), \\
dp[i][7][l][r] &= dp[i - 1][6][l][r] + \sum_{j = l}^r p[i][j], \\
dp[i][8] &= \max(\max_{1 \leq l \leq r \leq n} dp[i - 1][7][l][r],
                       dp[i - 1][8]), \\
dp[i][9][l][r] &= \max(dp[i - 1][8],
                       dp[i - 1][9][l][r]) + (p[i][l] + p[i][r]), \\
dp[i][10][l][r] &= \max(dp[i - 1][9][l][r], dp[i - 1][10][l][r]) +
                   \sum_{j = l}^r p[i][j], \\
dp[i][11][l][r] &= \max(dp[i - 1][10][l][r], dp[i - 1][11][l][r]) +
                   (p[i][l] + p[i][r]), \\
dp[i][12] &= \max(\max_{1 \leq l \leq r \leq n} dp[i - 1][11][l][r],
                        dp[i - 1][12]),
\end{align*}
$$

This takes $O(N^2M)$ with the appropiate precomputation. To avoid memory limit,
you only need to store the previous and current column's dp, although I'm not
sure if this is needed.

# Implementation
```cpp
/**
 * @author n685
 * @brief
 * @date 2024-07-30 10:53:38
 *
 *
 */
#include <bits/stdc++.h>

#ifdef LOCAL
#include "dd/debug.h"
#else
#define dbg(...) 42
#define dbgP(...) 420
#define dbgRP(...) 420420
void nline() {}
void bar() {}
#endif

template <class T> constexpr T INF = T{};
template <std::floating_point T>
constexpr T INF<T> = std::numeric_limits<T>::infinity();
template <> constexpr int INF<int> = 0x3f3f3f3f; // 1061109567
template <>
constexpr int64_t INF<int64_t> = 0x3f3f3f3f3f3f3f3f; // 4557430888798830399

using Grid = std::vector<std::vector<int>>;
using Arr = std::array<Grid, 13>;

int n, m;
Grid p;
Grid ps;
Arr pre, cur;
int chmax(int &a, int b) {
  if (a < b) {
    a = b;
  }
  return a;
}
int psum(int xl, int xr, int yl, int yr) {
  return ps[xr][yr] - ps[xl - 1][yr] - ps[xr][yl - 1] + ps[xl - 1][yl - 1];
}
void tr0(int /*i*/) { cur[0][0][0] = pre[0][0][0]; }
void tr1(int i) {
  for (int l = 1; l <= n; ++l) {
    for (int r = l; r <= n; ++r) {
      chmax(cur[1][l][r],
            std::max(pre[0][0][0], pre[1][l][r]) + psum(i, i, l, r));
    }
  }
}
void tr2(int i) {
  // case 1: from 1
  {
    std::vector mx(n + 1, -INF<int>);
    for (int l = 2; l <= n; ++l) {
      for (int r = l; r <= n; ++r) {
        chmax(mx[r], pre[1][l - 1][r]);
        chmax(cur[2][l][r], mx[r] + psum(i, i, l, r));
      }
    }
  }
  // case 2: from 2
  std::vector<int> mxl(n + 1,
                       -INF<int>); // mx for each l' that satisfies r' >= r
  for (int r = n; r >= 1; --r) {
    for (int l = 1; l <= r; ++l) {
      chmax(mxl[l], pre[2][l][r]);
    }
    int mx = (r == n ? -INF<int> : mxl[r + 1]);
    for (int l = r; l >= 1; --l) {
      chmax(mx, mxl[l]);
      chmax(cur[2][l][r], mx + psum(i, i, l, r));
    }
  }
}
void tr3(int i) {
  for (int l = 1; l <= n; ++l) {
    int mx = -INF<int>;
    for (int r = l; r <= n; ++r) {
      chmax(cur[3][l][r], std::max(pre[3][l][r], mx) + psum(i, i, l, r));
      chmax(mx, pre[2][l][r]);
    }
  }
}
void tr4(int /*i*/) {
  chmax(cur[4][0][0], pre[4][0][0]);
  for (int l = 1; l <= n; ++l) {
    for (int r = l; r <= n; ++r) {
      chmax(cur[4][0][0], pre[3][l][r]);
    }
  }
}
void tr5(int i) {
  for (int l = 1; l <= n; ++l) {
    for (int r = l + 2; r <= n; ++r) {
      chmax(cur[5][l][r], pre[4][0][0] + psum(i, i, l, r));
    }
  }
}
void tr6(int i) {
  for (int l = 1; l <= n; ++l) {
    for (int r = l + 2; r <= n; ++r) {
      chmax(cur[6][l][r],
            std::max(pre[6][l][r], pre[5][l][r]) + p[i][l] + p[i][r]);
    }
  }
}
void tr7(int i) {
  for (int l = 1; l <= n; ++l) {
    for (int r = l + 2; r <= n; ++r) {
      chmax(cur[7][l][r], pre[6][l][r] + psum(i, i, l, r));
    }
  }
}
void tr8(int /*i*/) {
  chmax(cur[8][0][0], pre[8][0][0]);
  for (int l = 1; l <= n; ++l) {
    for (int r = l; r <= n; ++r) {
      chmax(cur[8][0][0], pre[7][l][r]);
    }
  }
}
void tr9(int i) {
  for (int l = 1; l <= n; ++l) {
    for (int r = l + 2; r <= n; ++r) {
      chmax(cur[9][l][r],
            std::max(pre[8][0][0], pre[9][l][r]) + p[i][l] + p[i][r]);
    }
  }
}
void tr10(int i) {
  for (int l = 1; l <= n; ++l) {
    for (int r = l + 2; r <= n; ++r) {
      chmax(cur[10][l][r],
            std::max(pre[9][l][r], pre[10][l][r]) + psum(i, i, l, r));
    }
  }
}
void tr11(int i) {
  for (int l = 1; l <= n; ++l) {
    for (int r = l + 2; r <= n; ++r) {
      chmax(cur[11][l][r],
            std::max(pre[10][l][r], pre[11][l][r]) + p[i][l] + p[i][r]);
    }
  }
}
void tr12(int /*i*/) {
  chmax(cur[12][0][0], pre[12][0][0]);
  for (int l = 1; l <= n; ++l) {
    for (int r = l; r <= n; ++r) {
      chmax(cur[12][0][0], pre[11][l][r]);
    }
  }
}

int main() {
#ifndef LOCAL
  std::ios::sync_with_stdio(false);
  std::cin.tie(nullptr);
#endif

  std::cin >> n >> m;
  p = std::vector(m + 1, std::vector<int>(n + 1));
  for (int j = n; j >= 1; --j) {
    for (int i = 1; i <= m; ++i) {
      std::cin >> p[i][j];
    }
  }
  ps = p;
  for (int i = 1; i <= m; ++i) {
    for (int j = 1; j <= n; ++j) {
      ps[i][j] += ps[i - 1][j] + ps[i][j - 1] - ps[i - 1][j - 1];
    }
  }

  pre.fill(std::vector(n + 1, std::vector(n + 1, -INF<int>)));
  cur.fill(std::vector(n + 1, std::vector(n + 1, -INF<int>)));
  pre[0][0][0] = 0;
  for (int i = 1; i <= m; ++i) {
    tr0(i);
    tr1(i);
    tr2(i);
    tr3(i);
    tr4(i);
    tr5(i);
    tr6(i);
    tr7(i);
    tr8(i);
    tr9(i);
    tr10(i);
    tr11(i);
    tr12(i);
    pre.swap(cur);
    cur.fill(std::vector(n + 1, std::vector(n + 1, -INF<int>)));
  }
  int ans = pre[12][0][0];
  for (int l = 1; l <= n; ++l) {
    for (int r = l; r <= n; ++r) {
      chmax(ans, pre[11][l][r]);
    }
  }
  std::cout << ans << '\n';
}
```
