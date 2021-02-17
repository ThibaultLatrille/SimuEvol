#pragma once

#include <random>
// Random generator engine with seed 0.
double seed{0};
std::default_random_engine generator(seed);
std::default_random_engine generator_brownian(seed);

std::normal_distribution<double> normal_distrib(0.0, 1.0);