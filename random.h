//
// Created by raphr on 16/11/2018.
//

#ifndef EBARC_RANDOM_H
#define EBARC_RANDOM_H

#include <random>

namespace rnd {
    extern std::mt19937_64 rng;
    unsigned int set_seed();
    unsigned int set_seed(const unsigned int);
    bool bernoulli(const double = 0.5);
    size_t binomial(const size_t, const double = 0.5);
    size_t random_int(const size_t);
    size_t poisson(const double);
    double normal(const double, const double);
}

#endif //EBARC_RANDOM_H
