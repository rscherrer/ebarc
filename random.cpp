//
// Created by raphr on 16/11/2018.
// This file contains useful random number generating functions
//

#include "random.h"
#include <random>
#include <chrono>

namespace rnd {

    // Random number generator
    std::mt19937_64 rng;

    // Function to set a seed
    // without input argument
    unsigned int set_seed() {
        unsigned int seed = static_cast<unsigned int>(std::chrono::high_resolution_clock::now().time_since_epoch().count());
        return seed;
    }

    // with input argument
    unsigned int set_seed(const unsigned int seed) {
        rng.seed(seed);
        return seed;
    }

    // Bernoulli event
    bool bernoulli(const double p) {
        return std::bernoulli_distribution(p)(rng);
    }

    // Binomial distribution
    size_t binomial(size_t n, const double p) {
        return std::binomial_distribution<size_t>(n,p)(rng);
    }

    // Sample a random integer
    size_t random_int(const size_t n) {
        return std::uniform_int_distribution<size_t>(0u, n - 1u)(rng);
    }

    // Sample in a Poisson distribution
    size_t poisson(const double k) {
        return std::poisson_distribution<size_t>(k)(rng);
    }

    // Sample from a normal distribution
    double normal(const double mu, const double sigma) {
        return sigma == 0.0 ? mu : std::normal_distribution<double>(mu, sigma)(rng);
    }

}


