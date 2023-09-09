//
// Random Number Generation Library
//
// Changmook Chun. (c) 2003--2023.
//
// Parts of this library is written by D. E. Knuth.
//

#ifndef __RANDOM_H_
#define __RANDOM_H_

#include <memory>

namespace random_number {
// uniform dist.(int)
int dist_uniform_int (int*, int n);

// uniform dist.(double)
int dist_uniform_double (double*, int n);

// normal dist.
int dist_normal_polar_rejection (double*, int n);

// gamma dist.
int dist_gamma (double*, int n, double a, double r);

// Brownian bridge
int brownian_bridge_bisection (
    double*, double*, double z_0, double z_N, int N
);
} // namespace RandomNumber

#endif // __RANDOM_H_
