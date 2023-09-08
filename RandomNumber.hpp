//
// Random Number Generation Library
//
// Changmook Chun. (c) 2003. All rights reserved.
//
// >> IMPORTANT <<
// Parts of this library is written by D.E.Knuth.
//

#ifndef __RANDOM_H_
#define __RANDOM_H_

#include <memory>

namespace random_number {
// uniform dist.(int)
int dist_uniform_int (std::unique_ptr<int[]>&, int n);

// uniform dist.(double)
int dist_uniform_double (std::unique_ptr<double[]>&, int n);

// normal dist.
int dist_normal_polar_rejection (std::unique_ptr<double[]>&, int n);

// gamma dist.
int dist_gamma (std::unique_ptr<double[]>&, int n, double a, double r);

// Brownian bridge
int brownian_bridge_bisection (
    std::unique_ptr<double[]>&, std::unique_ptr<double[]>&, double z_0, double z_N, int N
);
} // namespace RandomNumber

#endif // __RANDOM_H_
