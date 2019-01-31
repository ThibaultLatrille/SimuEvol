#pragma once

#include <random>
// Random generator engine with seed 0.
double seed{0};
std::default_random_engine generator(seed);

std::normal_distribution<double> normal_distrib(0.0, 1.0);