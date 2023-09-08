#include <iostream>
#include <memory>
#include "random_number.hpp"

int main (int argc, const char * argv[]) {
    constexpr  size_t dim = 8000;
    std::unique_ptr<double[]> aa = std::make_unique<double[]>(dim);
    
    random_number::dist_normal_polar_rejection(aa, dim);
    
    for (size_t i = 0; i != dim; i++) {
        std::cout << aa[i] << std::endl;
    }
    
    return 0;
}
