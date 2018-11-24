//
// Created by raphr on 16/11/2018.
//

#ifndef EBARC_INDIVIDUAL_H
#define EBARC_INDIVIDUAL_H

#include <array>

// Class individual

// Attributes
// ecological trait value
// attack rates
// habitat
// ecotype
// Functions
// disperse

// Inline square function
inline double sqr(double x) { return x * x;}

class individual {

public:

    // Constructors
    individual(individual const * parent);
    individual();

    // Action functions
    void disperse() const; // dispersion function

    // Functions to access private variables
    double get_trait_value() const;
    bool get_habitat() const;
    bool get_ecotype() const;
    double get_payoff() const;
    std::pair<double, double> get_attack_rates() const;

    // Functions to set mutable private variables
    void set_ecotype(const std::pair<double, double> &ecotypeCutoff) const;
    void set_payoff(const std::pair<double, double> &resources) const;

private:
    double ecologicalTraitValue;
    std::pair<double, double> attackRates;
    mutable bool habitat;
    mutable bool ecotype;
    mutable double payoff;
};

// Define an alias for a constant pointer to an individual
typedef individual const * pointer_ind;

// Function to sort individuals along the trade-off curve
bool sort_along_tradeoff (const std::pair<double, double> &attackRates1, const std::pair<double, double> &attackRates2);

#endif //EBARC_INDIVIDUAL_H
