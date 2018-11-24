//
// Created by raphr on 16/11/2018.
//

#include "individual.h"
#include "random.h"
#include <cmath>

// Call to external variables
extern double ecoSelCoeff;
extern double mutationRate;
extern double mutationalEffect;

// Constructor functions
individual::individual(pointer_ind parent) {

    // Initialize trait value
    ecologicalTraitValue = parent->get_trait_value();

    // Mutate
    bool isMutation = rnd::bernoulli(mutationRate);
    if(isMutation) {
        ecologicalTraitValue = rnd::normal(ecologicalTraitValue, mutationalEffect);
    }

    // Initialize attack rates
    attackRates.first  = exp(-ecoSelCoeff * sqr(ecologicalTraitValue + 1.0));
    attackRates.second = exp(-ecoSelCoeff * sqr(ecologicalTraitValue - 1.0));

    // Initialize habitat
    habitat = parent->get_habitat();

}

// Constructor used for the initial population
individual::individual() {

    // Initialize trait value
    ecologicalTraitValue = 1.0;

    // Initialize attack rates
    attackRates.first  = exp(-ecoSelCoeff * sqr(ecologicalTraitValue + 1.0));
    attackRates.second = exp(-ecoSelCoeff * sqr(ecologicalTraitValue - 1.0));

    // Initialize habitat
    habitat = true;
}

// Dispersion function
void individual::disperse() const {

    habitat = !habitat; // flip to the other habitat

}

// Functions to access private variables
double individual::get_trait_value() const {
    return(ecologicalTraitValue);
}

bool individual::get_habitat() const {
    return(habitat);
}

bool individual::get_ecotype() const {
    return(ecotype);
}

double individual::get_payoff() const {
    return(payoff);
}

std::pair<double, double> individual::get_attack_rates() const {
    return(attackRates);
}

// Functions to set mutable private variables
void individual::set_ecotype(const std::pair<double, double> &ecotypeCutoff) const {

    // If the individual is further than the cutoff along the tradeoff curve, it is a resource 1 specialist (ecotype 0)
    // otherwise it is a resource 2 specialist (ecotype 1)
    ecotype = !sort_along_tradeoff(attackRates, ecotypeCutoff);

}

void individual::set_payoff(const std::pair<double, double> &resources) const {

    payoff = ecotype ? attackRates.second * resources.second : attackRates.first * resources.first;

}

// Function to sort individuals along the trade-off curve
bool sort_along_tradeoff (const std::pair<double, double> &attackRates1, const std::pair<double, double> &attackRates2) {

    // Locate individuals along the trade-off curve
    bool isInd1OnLeftSide = attackRates1.second > attackRates1.first;
    bool isInd2OnLeftSide = attackRates2.second > attackRates2.first;

    // Determine whether individual 1 is before individual 2 along the trade-off curve
    if(isInd1OnLeftSide & isInd2OnLeftSide)
        return (attackRates1.first > attackRates2.first);
    else if(!isInd1OnLeftSide & !isInd2OnLeftSide)
        return (attackRates1.second < attackRates2.second);
    else
        return (isInd1OnLeftSide & !isInd2OnLeftSide);

}