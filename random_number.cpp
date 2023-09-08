//
// Random Number Generation Library
//
// Changmook Chun. (c) 2003--2023.
//
// Parts of this library is written by D. E. Knuth.
//

#include "random.h"
#include "random_number.hpp"
#include <math.h>
#include <stdlib.h>
#include <time.h>

#define RANF_ARR_NEXT()                                                        \
    (*_ranf_arr_ptr >= 0 ? *_ranf_arr_ptr++ : ranf_arr_cycle())

namespace random_number {

constexpr int KK = 100;        // the long lag
constexpr int LL = 37;         // the short lag
constexpr int MM = (1L << 30); // the modulus
constexpr int QUALITY = 1009;  // recommended quality level for high-res use
constexpr int TT = 70;         // guaranteed separation between streams

auto _ran_u = std::make_unique<double[]> (KK); // the generator state
auto _ranf_arr_buf = std::make_unique<double[]> (QUALITY);
double _ranf_arr_sentinel = -1.0;
double* _ranf_arr_ptr = &_ranf_arr_sentinel; // the next random fraction, or -1

auto _ran_x = std::make_unique<int[]> (KK); // the generator state
auto _ran_arr_buf = std::make_unique<int[]> (QUALITY);
int _ran_arr_sentinel = -1;
int* _ran_arr_ptr = &_ran_arr_sentinel; // the next random number, or -1

// *****************************************************************************
//                                                              Utility routines
// *****************************************************************************

// Units bit of x
bool is_odd (const int x) { return x & 1; }

// Subtraction mod MM
int mod_diff (const int x, const int y) { return (x - y) & (MM - 1); }

// (x+y) mod 1.0
double mod_sum (const double x, const double y) {
    return (x + y) - static_cast<int> (x + y);
}

// Forward Declarations of private functions
void _brownian_bridge_bisect (int j, int J, std::unique_ptr<double[]>&, std::unique_ptr<double[]>&, std::unique_ptr<double[]>&);

void _ran_array (std::unique_ptr<int[]>&, int n);
void _ran_start (int32_t seed);
int _ran_arr_cycle ();
void _ranf_array (std::unique_ptr<double[]>&, int n);
void _ranf_start (int32_t seed);
double _ranf_arr_cycle ();

/* It is explained in Seminumerical Algorithms, 3rd edition, Section 3.6
 * (or in the errata to the 2nd edition --- see
 *  http://www-cs-faculty.stanford.edu/~knuth/taocp.html
 * in the changes to Volume 2 on pages 171 and following).
 */

/* N.B. The MODIFICATIONS introduced in the 9th printing (2002) are
 * included here; there's no backwards compatibility with the original.
 */

/* This version also adopts Brendan McKay's suggestion to
 * accommodate naive users who forget to call _ran_start(seed).
 */

void _ranf_array (std::unique_ptr<double[]>& aa, int n) {
    int i, j;
    for (j = 0; j < KK; j++)
        aa[j] = _ran_u[j];
    for (; j < n; j++)
        aa[j] = mod_sum (aa[j - KK], aa[j - LL]);
    for (i = 0; i < LL; i++, j++)
        _ran_u[i] = mod_sum (aa[j - KK], aa[j - LL]);
    for (; i < KK; i++, j++)
        _ran_u[i] = mod_sum (aa[j - KK], _ran_u[i - LL]);
}

/* the following routines are adapted from exercise 3.6--15
 * after calling _ranf_start, get new randoms by, e.g., "x=RANF_ARR_NEXT()"
 */

void _ranf_start (int seed) {
    int t, s, j;
    auto u = std::make_unique<double[]> (KK + KK - 1);
    double ulp = (1.0 / (1L << 30)) / (1L << 22); // 2 to the -52
    double ss = 2.0 * ulp * ((seed & 0x3fffffff) + 2);

    for (j = 0; j < KK; j++) {
        u[j] = ss; // bootstrap the buffer
        ss += ss;
        if (ss >= 1.0)
            ss -= 1.0 - 2 * ulp; // cyclic shift of 51 bits
    }
    u[1] += ulp; // make u[1] (and only u[1]) "odd"
    for (s = seed & 0x3fffffff, t = TT - 1; t;) {
        for (j = KK - 1; j > 0; j--)
            u[j + j] = u[j], u[j + j - 1] = 0.0; // "square"
        for (j = KK + KK - 2; j >= KK; j--) {
            u[j - (KK - LL)] = mod_sum (u[j - (KK - LL)], u[j]);
            u[j - KK] = mod_sum (u[j - KK], u[j]);
        }
        if (is_odd (s)) { // "multiply by z"
            for (j = KK; j > 0; j--)
                u[j] = u[j - 1];
            u[0] = u[KK]; // shift the buffer cyclically
            u[LL] = mod_sum (u[LL], u[KK]);
        }
        if (s)
            s >>= 1;
        else
            t--;
    }
    for (j = 0; j < LL; j++)
        _ran_u[j + KK - LL] = u[j];
    for (; j < KK; j++)
        _ran_u[j - LL] = u[j];
    for (j = 0; j < 10; j++)
        _ranf_array (u, KK + KK - 1); // warm things up
    _ranf_arr_ptr = &_ranf_arr_sentinel;
}

double _ranf_arr_cycle () {
    if (_ranf_arr_ptr != &_ranf_arr_sentinel) {
        _ranf_array (_ranf_arr_buf, QUALITY);
    } else {
        _ranf_start (314159LL); // the user forgot to initialize
    }
    _ranf_arr_buf[100] = -1;
    _ranf_arr_ptr = _ranf_arr_buf.get() + 1;
    return _ranf_arr_buf[0];
}

void _ran_array (std::unique_ptr<int[]>& aa, int n) {
    int i, j;
    for (j = 0; j < KK; j++)
        aa[j] = _ran_x[j];
    for (; j < n; j++)
        aa[j] = mod_diff (aa[j - KK], aa[j - LL]);
    for (i = 0; i < LL; i++, j++)
        _ran_x[i] = mod_diff (aa[j - KK], aa[j - LL]);
    for (; i < KK; i++, j++)
        _ran_x[i] = mod_diff (aa[j - KK], _ran_x[i - LL]);
}

void _ran_start (int seed) {
    int t, j;
    auto x = std::make_unique<int[]> (KK + KK - 1); // the preparation buffer
    int ss = (seed + 2) & (MM - 2);
    for (j = 0; j < KK; j++) {
        x[j] = ss; // bootstrap the buffer
        ss <<= 1;
        if (ss >= MM)
            ss -= MM - 2; // cyclic shift 29 bits
    }
    x[1]++; // make x[1] (and only x[1]) odd
    for (ss = seed & (MM - 1), t = TT - 1; t;) {
        for (j = KK - 1; j > 0; j--)
            x[j + j] = x[j], x[j + j - 1] = 0; // "square"
        for (j = KK + KK - 2; j >= KK; j--)
            x[j - (KK - LL)] = mod_diff (x[j - (KK - LL)], x[j]),
                        x[j - KK] = mod_diff (x[j - KK], x[j]);
        if (is_odd (ss)) { // "multiply by z"
            for (j = KK; j > 0; j--)
                x[j] = x[j - 1];
            x[0] = x[KK]; // shift the buffer cyclically
            x[LL] = mod_diff (x[LL], x[KK]);
        }
        if (ss)
            ss >>= 1;
        else
            t--;
    }
    for (j = 0; j < LL; j++)
        _ran_x[j + KK - LL] = x[j];
    for (; j < KK; j++)
        _ran_x[j - LL] = x[j];
    for (j = 0; j < 10; j++)
        _ran_array (x, KK + KK - 1); // warm things up
    _ran_arr_ptr = &_ran_arr_sentinel;
}

int _ran_arr_cycle () {
    if (_ran_arr_ptr != &_ran_arr_sentinel) {
        _ran_array (_ran_arr_buf, QUALITY);
    } else {
        _ran_start (314159LL); // the user forgot to initialize
    }
    _ran_arr_buf[100] = -1;
    _ran_arr_ptr = _ran_arr_buf.get() + 1;
    return _ran_arr_buf[0];
}

// uniform dist.(int)
// This function generate n uniformly distributed random integers.
int dist_uniform_int (std::unique_ptr<int[]>& aa, int n) {
    _ran_start (static_cast<int32_t> (time (NULL)));
    _ran_array (aa, n);
    return 0;
}

// uniform dist.(double)
// This function generates n uniformly distributed random real numbers.
int dist_uniform_double (std::unique_ptr<double[]>& aa, int n) {
    _ranf_start (static_cast<int32_t> (time (NULL)));
    _ranf_array (aa, n);
    return 0;
}

int dist_normal_polar_rejection (std::unique_ptr<double[]>& aa, int n) {
    /*
     This function generates n normally distributed random numbers.
     The mean and standard deviation of the distribution is
     0 and 1, respectively.

     To see the detailed description of this algorithm,
     consult D.E.Knuth's ``The Art of Computer Programming Volume 2''
     section 3.4.1 on pp.119.

     aa[]: array that will contain generated random sequence
     n: int that tells the size of aa[]
     */
    int m, i, idx = 0;
    double v1, v2, s, tmp;
    int done = 0;

    if (is_odd (n)) {
        m = n + 1;
    } else {
        m = n;
    }
    auto _a = std::make_unique<double[]> (2 * m); // uniform deviates
    // NOTE: _a is twice as large as _aa.
    // see D.E.Knuth(pp.122) for the rationale.
    auto _aa = std::make_unique<double[]> (m); // normal deviates

    while (!done) {
        // build an array of uniformly random numbers
        // _ranf_start((int)time(NULL));
        _ranf_start (static_cast<int32_t> (time (NULL)));
        _ranf_array (_a, 2 * m);

        // calculate normally distributed numbers
        // from the uniformly random numbers
        for (i = 0; i < m; i++) {
            v1 = _a[i * 2] * 2. - 1.;
            v2 = _a[i * 2 + 1] * 2. - 1.;
            s = v1 * v1 + v2 * v2;
            if (s >= 1)
                continue;

            tmp = sqrt (-2 * log (s) / s);
            _aa[idx++] = v1 * tmp;
            _aa[idx++] = v2 * tmp;

            if (idx > m) {
                done = 1;
                break;
            }
        }
    }
    for (i = 0; i < n; i++) {
        aa[i] = _aa[i];
    }

    return 0;
}

int dist_gamma (std::unique_ptr<double[]>& aa, int n, double a, double r) {
    /*
     int dist_gamma(double aa[], int n, double a, double r)
     generates random deviates from gamma distribution

     Function
     Generates random deviates from the gamma distribution whose
     density is

     (A**R)/Gamma(R) * X**(R-1) * Exp(-A*X)

     Arguments
     a --> Location parameter of Gamma distribution
     r --> Shape parameter of Gamma distribution

     Method
     Renames SGAMMA from TOMS as slightly modified by BWB to use RANF
     instead of SUNIF.

     For details see:
     (Case R >= 1.0)
     Ahrens, J.H. and Dieter, U.
     Generating Gamma Variates by a
     Modified Rejection Technique.
     Comm. ACM, 25,1 (Jan. 1982), 47 - 54.

     Algorithm GD
     (Case 0.0 <= R <= 1.0)
     Ahrens, J.H. and Dieter, U.
     Computer Methods for Sampling from Gamma,
     Beta, Poisson and Binomial Distributions.
     Computing, 12 (1974), 223-246/
     Adapted algorithm GS.
     */

    for (int i = 0; i < n; i++) {
        aa[i] = gengam (a, r);
    }

    return 0;
}

int brownian_bridge_bisection (
    std::unique_ptr<double[]>& tt,
    std::unique_ptr<double[]>& zz,
    double z_0,
    double z_N,
    int N
) {
    /*
     Constructs a Brownian bridge from z_0 to z_N
     this function uses bisection construction method that gives
     ``better control'' over the generated sample paths.

     tt[]: array of time history. from 0 to N, i.e., size = N+1.
     zz[]: array of z_t (Wiener process). size = N+1.
     z_0: starting point of z_t at tt[0].
     z_N: destination of z_t at tt[N].
     N: the number of steps in tt and zz.
     */

    // prepare data
    // 1. build an array containing epsilon
    auto eps = std::make_unique<double[]> (N + 1); // array of epsilon
    _ranf_start ((int)time (NULL));
    _ranf_array (eps, N + 1); // normal, mean = 0, std = 1

    // 2. put z_0 and z_N at both ends of array, zz
    zz[0] = z_0;
    zz[N] = z_N;

    // Now, all the other elements of zz will be filled by
    // recursive bisection routine.

    // Call recursive bisection routine
    _brownian_bridge_bisect (0, N, tt, zz, eps);

    return 0;
}

void _brownian_bridge_bisect (
    int j,
    int J,
    std::unique_ptr<double[]>& tt,
    std::unique_ptr<double[]>& zz,
    std::unique_ptr<double[]>& eps
) {
    /*
     Recursive bisection routine called by brownian_bridge_bisection

     int j;         left index, i.e., j- in the text
     int J;         right index, i.e., j+ in the text
     double* tt;     array of time
     double* zz;     array of z
     double* eps;     array of epsilon, standard normal random sequence
     */

    int jm;              // mid-way between j and J
    double td, tld, trd; // time differences
    // td: t_(j+) - t_(j-)
    // tld: t_(jm) - t_(j-)
    // trd: t_(j+) - t_(jm)

    // Evaluate the condition to end the recursion
    if (J == j + 1)
        return;

    // calculate index, jm, the point where bisection occurs
    jm = (J - j) / 2 + j;

    // calculate time differences
    // this enables more efficient calculation of the bisection
    td = tt[J] - tt[j];
    tld = tt[jm] - tt[j];
    trd = tt[J] - tt[jm];

    // the result of the bisection
    zz[jm] = (trd * zz[j] + tld * zz[J]) / td + sqrt (trd * tld / td) * eps[jm];

    // now, call this bisection method (recursion) for the left part of jm
    _brownian_bridge_bisect (j, jm, tt, zz, eps);

    // now, call this bisection method (recursion) for the right part of jm
    _brownian_bridge_bisect (jm, J, tt, zz, eps);
}

} // namespace random_number
