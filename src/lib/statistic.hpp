#pragma once

#include <algorithm>
#include <vector>

// Function for sum
double sum(std::vector<double> const &v) { return accumulate(v.begin(), v.end(), 0.0); }

// Function for average
double avg(std::vector<double> const &v) { return sum(v) / v.size(); }

// Function for variance
double variance(std::vector<double> const &v) {
    double s2 = accumulate(v.begin(), v.end(), 0.0, [](double a, double b) { return a + b * b; });
    double mean = avg(v);
    return (s2 / v.size()) - mean * mean;
}


// Function for sum
double sum(std::vector<unsigned> const &v) {
    return accumulate(v.begin(), v.end(), static_cast<unsigned>(0));
}

// Function for mean
double mean(std::vector<unsigned> const &v) {
    double return_value = 0.0;
    for (unsigned i = 0; i < v.size(); i++) { return_value += i * v[i]; }
    return return_value / sum(v);
};