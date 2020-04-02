#pragma once

#include <numeric>
#include <vector>

class SummaryStatistic {
    double long abs_s{0}, s{0}, s2{0};
    u_long size{0};

  public:
    SummaryStatistic() = default;

    void add(double const &v) {
        if (v > 1e12) { std::cerr << s << std::endl; }
        s += v;
        abs_s += std::abs(v);
        s2 += v * v;
        size++;
    }

    double mean() const { return static_cast<double>(s / size); }
    double variance() const { return static_cast<double>(s2 / size) - mean() * mean(); }
    double std() const { return sqrt(variance()); }

    double abs_mean() const { return static_cast<double>(abs_s / size); }
    double abs_variance() const { return static_cast<double>(s2 / size) - abs_mean() * abs_mean(); }
    double abs_std() const { return sqrt(abs_variance()); }
};

// Function for sum
double sum(std::vector<double> const &v) { return std::accumulate(v.begin(), v.end(), 0.0); }

// Function for average
double avg(std::vector<double> const &v) { return sum(v) / v.size(); }

// Function for variance
double variance(std::vector<double> const &v) {
    double s2 =
        std::accumulate(v.begin(), v.end(), 0.0, [](double a, double b) { return a + b * b; });
    double mean = avg(v);
    return (s2 / v.size()) - mean * mean;
}

// Variance of a vector
double variance(std::vector<double> const &v, double mean) {
    double s2 =
        std::accumulate(v.begin(), v.end(), 0.0, [](double a, double b) { return a + b * b; });
    return (s2 / v.size()) - mean * mean;
}

// Function for sum
double sum(std::vector<u_long> const &v) {
    return accumulate(v.begin(), v.end(), static_cast<u_long>(0));
}

// Function for mean
double mean(std::vector<u_long> const &v) {
    double return_value = 0.0;
    for (u_long i = 0; i < v.size(); i++) { return_value += i * v[i]; }
    return return_value / sum(v);
}