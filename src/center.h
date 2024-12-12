#ifndef CENTER_H
#define CENTER_H

#include <vector>
#include <unordered_map>
#include <string>
#include <cmath>

void centerConversion(
    std::vector<double>& beta, 
    std::vector<double>& mb_v,
    std::vector<double>& rb_v,
    bool robust,
    std::unordered_map<std::string, double> meanValues, 
    std::vector<std::string> interaction_names, 
    size_t dim,
    size_t nExp,
    int centerIn,
    int centerOut);

#endif // CENTER_H
