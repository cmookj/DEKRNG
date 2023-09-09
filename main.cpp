#include <algorithm>
#include <iostream>
#include <vector>
#include "random_number.hpp"

int main (int argc, const char * argv[]) {
    constexpr  size_t dim = 8000;
    std::vector<double> aa(dim);
    
    random_number::dist_normal_polar_rejection(aa.data(), dim);
    
    std::for_each (aa.cbegin(), aa.cend(), [](const auto& a){
        std::cout << a << '\n';
    });
    
    return 0;
}
