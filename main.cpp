#include <algorithm>
#include <iomanip>
#include <iostream>
#include <vector>
#include "random_number.hpp"

int main (int argc, const char * argv[]) {
    constexpr  size_t dim = 2009;
    std::vector<int32_t> aa_int(dim);
    
    // random_number::dist_normal_polar_rejection(aa.data(), dim);
    random_number::random_int_uniform rng_int{310952L};
    for (int m = 0; m < 2009; ++m) rng_int.ran_array(aa_int.data(), 1009);
    std::cout << rng_int.ran_x(0) << '\n'; // 995235265
    
    rng_int.ran_start(310952L);
    for (int m = 0; m < 1009; ++m) rng_int.ran_array(aa_int.data(), 2009);
    std::cout << rng_int.ran_x(0) << '\n'; // 995235265
    
//    std::for_each (aa_int.cbegin(), aa_int.cbegin() + 20, [](const auto& a){
//        std::cout << a << '\n';
//    });
    
    std::vector<double> aa_dbl(dim);
    random_number::random_double_uniform rng_dbl{310952L};
    for (int m = 0; m < 2009; ++m) rng_dbl.ranf_array(aa_dbl.data(), 1009);
    std::cout << std::setprecision(20);
    std::cout << rng_dbl.ran_u(0) << '\n'; // 0.36410514377569680455
    
    rng_dbl.ranf_start(310952L);
    for (int m = 0; m < 1009; ++m) rng_dbl.ranf_array(aa_dbl.data(), 2009);
    std::cout << rng_dbl.ran_u(0) << '\n'; // 0.36410514377569680455
    
    return 0;
}
