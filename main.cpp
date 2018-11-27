#include <iostream>
#include <cstdlib>
#include <set>
#include <list>
#include <utility>
#include <cassert>
#include <fstream>
#include <sstream>
#include <string>
#include <exception>
#include "random.h"
#include "individual.h"

/*=======================================================================================================
                                           Global parameters
========================================================================================================*/

// Global parameters
double dispersalRate = 0.001;
size_t initialPopSize = 100u;
double outflow = 0.001;
double habitatAsymmetry = 0.5;
double survivalRate = 0.5;
double basalGrowth = 4.0;
size_t Tmax = 10000u;
size_t timeToRecord = 100u;
size_t timeForScreenshot = 1000u;

// Global constants that are externally called
double mutationRate = 0.001;
double mutationalEffect = 0.01;
double ecoSelCoeff = 1.0;

std::ofstream logFile, datFile;

// Set a seed
unsigned int seed;

/*=======================================================================================================
                                 Population-level life-cycle functions
========================================================================================================*/

// Dispersal function
void dispersal(std::list<pointer_ind> &population) {

    // If dispersal rate is high, sample binary events for each individual in the population
    if(dispersalRate > 0.5) {
        for(pointer_ind ind : population)
            if(rnd::bernoulli(dispersalRate)) ind -> disperse();
    }

    // If not, use a binomial to sample the number of dispersers
    else {
        const size_t popSize = population.size();
        size_t nDispersers = rnd::binomial(popSize, dispersalRate);
        if(nDispersers == 0) return;
        std::set<size_t> dispersers;

        // Sample each disperser
        while(dispersers.size() < nDispersers)
            dispersers.insert(rnd::random_int(popSize));

        // Implement dispersion throughout the set of dispersers
        auto it = population.cbegin();
        size_t j = 0u;
        for(size_t i : dispersers) {
            std::advance(it, i - j);
            (*it)->disperse();
            j = i;
        }
    }
}

// Resource utilization, where individual receive their payoff
void resource_utilization(std::list<pointer_ind> &population) {

    // Sort individuals by habitat
    auto iti = population.begin();
    for(auto itj = population.end(); iti != itj;) {
        // If individual belongs to habitat 0, keep it here
        if(!(*iti)->get_habitat()) {
            ++iti;
        }

        // Otherwise, if it belongs to habitat 1, swap it towards the end of the list
        else {
            --itj;
            std::swap(*iti, *itj);
        }
    }

    auto habitatLimit = iti; // index of the first individual in habitat 1

    // Prepare to store attack rates within habitat
    std::list<std::pair<double, double>> allAttackRates;

    // Repeat within each habitat
    for(size_t hab = 0u; hab < 2u; ++hab) {

        // Prepare to store consumed resources
        std::pair<double, double> consumedResources {0.0, 0.0};
        std::pair<double, double> resources;

        // Set iterator limits matching the current habitat
        auto from = hab == 0u ? population.begin() : habitatLimit;
        auto to = hab == 0u ? habitatLimit : population.end();

        // Store the attack rates of all individuals
        for(auto it = from; it != to; ++it) {

            auto currentAttackRates = (*it) -> get_attack_rates();
            allAttackRates.push_back(currentAttackRates);

            // Everybody feeds on resource 1 for now, will be reassigned later
            consumedResources.first += currentAttackRates.first;

        }

        // Initialize equilibrium resource concentrations
        resources.first = (hab == 0u ? 1.0 : 1.0 - habitatAsymmetry) / (1 + outflow * consumedResources.first);
        resources.second = (hab == 1u ? 1.0 : 1.0 - habitatAsymmetry) / (1 + outflow * consumedResources.second);

        // Sort the list of attack rates
        allAttackRates.sort(sort_along_tradeoff);

        // Prepare to record the ecotype cutoff
        std::pair<double, double> ecotypeCutoff;

        // Loop through individuals in the local habitat and flip them to resource 2 until it is not advantageous to do so anymore
        for(std::pair<double, double> currentAttackRates : allAttackRates) {
            bool isFlippingAdvantageous = currentAttackRates.first * resources.first < currentAttackRates.second * resources.second;

            // If you get a better payoff by switching to the resource 2
            if(isFlippingAdvantageous) {

                // Update consumed resources
                consumedResources.first -= currentAttackRates.first;
                consumedResources.second += currentAttackRates.second;

                // Update equilibrium concentrations
                resources.first = (hab == 0u ? 1.0 : 1.0 - habitatAsymmetry) / (1 + outflow * consumedResources.first);
                resources.second = (hab == 1u ? 1.0 : 1.0 - habitatAsymmetry) / (1 + outflow * consumedResources.second);

                // Update the ecotype cut-off
                ecotypeCutoff.first = currentAttackRates.first;
                ecotypeCutoff.second = currentAttackRates.second;

            }
            // Break the loop when it is not advantageous anymore to switch resource (we found the ecotype cutoff)
            else
                break;

        } // end of loop through individuals in the focal habitat

        // Assign payoffs and ecotypes to all individuals
        for(auto it = from; it != to; ++it) {
            (*it)->set_ecotype(ecotypeCutoff);
            (*it)->set_payoff(resources);
        }

        // Reset habitat-specific variables
        consumedResources.first = consumedResources.second = resources.first = resources.second = 0.0;

    } // end of focal habitat

}

// Reproduction function
std::list<pointer_ind> reproduction(std::list<pointer_ind> &population) {

    // Each individual can reproduce clonally
    // The number of offspring is Poisson distributed with mean = payoff
    // Append each newborn to the population

    // Prepare to store the newborns
    std::list<pointer_ind> newpopulation;

    for(auto ind : population) {

        const double expectedOffspring = ind->get_payoff() * basalGrowth;
        double nOffspring = rnd::poisson(expectedOffspring);

        // Append the offspring to the population (with mutation)
        while(nOffspring) {
            newpopulation.push_back(new individual(ind));
            --nOffspring;
        }

    }

    return newpopulation;

}

// Survival function
void survival(std::list<pointer_ind> &population, std::list<pointer_ind> &newpopulation) {

    // Each individual has a certain probability of surviving to the next generation
    // Append the survivors to the newborns
    for(auto ind : population) {
        bool isSurvivor = rnd::bernoulli(survivalRate);
        if(isSurvivor)  newpopulation.push_back(ind);
    }

}

/*=======================================================================================================
                                 Input-output functions
========================================================================================================*/

// Function to read and set parameters from an input file
void read_parameters(const std::string &inputFilename) {

    std::clog << "Reading parameters from file " << inputFilename << std::endl;

    // Open parameter file
    std::ifstream ifs(inputFilename);
    if(!ifs.is_open())
        throw std::runtime_error("unable to open input file in read_parameters()");

    // Prepare to receive input
    std::string nextInput;

    // Read seed
    ifs >> nextInput;
    if(nextInput == "rng_seed_clock") seed = rnd::set_seed();
    else if(nextInput == "rng_seed_user") {
        ifs >> seed;
        rnd::set_seed(seed);
    }
    else throw std::logic_error("\'rng_seed_clock\' or \'rng_seed_user <arg>\' expected at first line of parameter file\n");

    // Read input until there is nothing to read anymore
    while(ifs >> nextInput) {
        if(nextInput == "dispersalRate") {
            ifs >> dispersalRate;
        } else if(nextInput == "initialPopSize") {
            ifs >> initialPopSize;
        } else if(nextInput == "outflow") {
            ifs >> outflow;
        } else if(nextInput == "habitatAsymmetry") {
            ifs >> habitatAsymmetry;
        } else if(nextInput == "survivalRate") {
            ifs >> survivalRate;
        } else if(nextInput == "basalGrowth") {
            ifs >> basalGrowth;
        } else if(nextInput == "Tmax") {
            ifs >> Tmax;
        } else if(nextInput == "timeToRecord") {
            ifs >> timeToRecord;
        } else if(nextInput == "timeForScreenshot") {
            ifs >> timeForScreenshot;
        } else if(nextInput == "mutationRate") {
            ifs >> mutationRate;
        } else if(nextInput == "mutationalEffect") {
            ifs >> mutationalEffect;
        } else if(nextInput == "ecoSelCoeff") {
            ifs >> ecoSelCoeff;
        }
    }

    std::clog << "Parameters succesfully read" << std::endl;
}

// Function to record population summary data
void record_data(const size_t &timestep, const std::list<pointer_ind> &population) {

    // Send population summary data to the data file
    datFile << timestep << ';';

    // Prepare to record pop sizes and mean trait values
    size_t popSize = population.size();
    size_t popSizeH0E0 = 0u;
    size_t popSizeH0E1 = 0u;
    size_t popSizeH1E0 = 0u;
    size_t popSizeH1E1 = 0u;
    double sumTraitValuesH0E0 = 0.0;
    double sumTraitValuesH0E1 = 0.0;
    double sumTraitValuesH1E0 = 0.0;
    double sumTraitValuesH1E1 = 0.0;
    double sumOfSquaresE0 = 0.0;
    double sumOfSquaresE1 = 0.0;
    double resourceH0E0 = 0.0;
    double resourceH0E1 = 0.0;
    double resourceH1E0 = 0.0;
    double resourceH1E1 = 0.0;

    // Loop through the population and count individuals and accumulate trait values
    // Also record one individual from each category, to later calculate resource concentrations from payoffs
    // Also accumulate sum of square trait values, to later calculate variance
    bool foundH0E0, foundH0E1, foundH1E0, foundH1E1;
    foundH0E0 = foundH0E1 = foundH1E0 = foundH1E1 = false;
    for(auto ind : population) {
        const bool isHabitat1 = ind->get_habitat();
        const bool isEcotype1 = ind->get_ecotype();
        if(isHabitat1)
            if(isEcotype1) {
                ++popSizeH1E1;
                sumTraitValuesH1E1 += ind->get_trait_value();
                sumOfSquaresE1 += ind->get_trait_value() * ind->get_trait_value();
                if(!foundH1E1) {
                    foundH1E1 = true;
                    resourceH1E1 = ind->get_payoff() / ind->get_attack_rates().second;
                }
            }
            else {
                ++popSizeH1E0;
                sumTraitValuesH1E0 += ind->get_trait_value();
                sumOfSquaresE0 += ind->get_trait_value() * ind->get_trait_value();
                if(!foundH1E0) {
                    foundH1E0 = true;
                    resourceH1E0 = ind->get_payoff() / ind->get_attack_rates().first;
                }
            }
        else
            if(isEcotype1) {
                ++popSizeH0E1;
                sumTraitValuesH0E1 += ind->get_trait_value();
                sumOfSquaresE1 += ind->get_trait_value() * ind->get_trait_value();
                if(!foundH0E1) {
                    foundH0E1 = true;
                    resourceH0E1 = ind->get_payoff() / ind->get_attack_rates().second;
                }
            }
            else {
                ++popSizeH0E0;
                sumTraitValuesH0E0 += ind->get_trait_value();
                sumOfSquaresE0 += ind->get_trait_value() * ind->get_trait_value();
                if(!foundH0E0) {
                    foundH0E0 = true;
                    resourceH0E0 = ind->get_payoff() / ind->get_attack_rates().first;
                }
            }
    }

    // Calculate the mean trait values with the accumulated values
    const double meanTraitValue = (sumTraitValuesH0E0 + sumTraitValuesH0E1 + sumTraitValuesH1E0 + sumTraitValuesH1E1) / popSize;
    const double meanTraitValueH0E0 = sumTraitValuesH0E0 / popSizeH0E0;
    const double meanTraitValueH0E1 = sumTraitValuesH0E1 / popSizeH0E1;
    const double meanTraitValueH1E0 = sumTraitValuesH1E0 / popSizeH1E0;
    const double meanTraitValueH1E1 = sumTraitValuesH1E1 / popSizeH1E1;

    // Calculate sample variance within each ecotype and total variance
    double varianceE0 = sumOfSquaresE0 - sqr(sumTraitValuesH0E0 + sumTraitValuesH1E0) / (popSizeH0E0 + popSizeH1E0);
    double varianceE1 = sumOfSquaresE1 - sqr(sumTraitValuesH0E1 + sumTraitValuesH1E1) / (popSizeH0E1 + popSizeH1E1);
    double totalVariance = (sumOfSquaresE0 + sumOfSquaresE1) - sqr(sumTraitValuesH0E0 + sumTraitValuesH1E0 + sumTraitValuesH0E1 + sumTraitValuesH1E1) / popSize;
    varianceE0 /= (popSizeH0E0 + popSizeH1E0 - 1.0);
    varianceE1 /= (popSizeH0E1 + popSizeH1E1 - 1.0);
    totalVariance /= (popSize - 1.0);

    // Calculate isolation statistics
    const double spatialIsolation = (popSizeH0E0 * popSizeH1E1 - popSizeH0E1 * popSizeH1E0) / sqrt((popSizeH0E0 + popSizeH0E1) * (popSizeH1E0 + popSizeH1E1) * (popSizeH0E0 + popSizeH1E0) * (popSizeH0E1 + popSizeH1E1));
    const double ecologicalIsolation = 1.0 - ((popSizeH0E0 + popSizeH1E0) * varianceE0 + (popSizeH0E1 + popSizeH1E1) * varianceE1) / (popSize * totalVariance);

    // Write to data file
    datFile << popSize << ';' << popSizeH0E0 << ';' << popSizeH0E1 << ';' << popSizeH1E0 << ';' << popSizeH1E1 << ';';
    datFile << meanTraitValue << ';' << meanTraitValueH0E0 << ';' << meanTraitValueH0E1 << ';' << meanTraitValueH1E0 << ';' << meanTraitValueH1E1 << ';';
    datFile << resourceH0E0 << ';' << resourceH0E1 << ';' << resourceH1E0 << ';' << resourceH1E1 << ';';
    datFile << spatialIsolation << ';' << ecologicalIsolation << ';';

    datFile << std::endl;

}

// Function to take a screenshot of the population
void take_screenshot(const size_t &timestep, const std::list<pointer_ind> &population) {

    // Open a file for the screenshot
    std::ofstream screenshotFile;
    std::ostringstream osscreenshot;
    osscreenshot << "screenshot_" << seed << '_' << timestep << ".csv";
    screenshotFile.open(osscreenshot.str());
    if(!screenshotFile.is_open())
        throw std::runtime_error("unable to open file");

    // Write the header of the screenshot file
    std::vector<std::string> headers = {
            "ecotype",
            "habitat",
            "attackRate1",
            "attackRate2",
            "traitValue"
    };
    for(auto currHeader : headers) {
        screenshotFile << currHeader << ';';
    }
    screenshotFile << std::endl;

    // Write down information for every individual
    for(auto ind : population) {
        screenshotFile << ind->get_ecotype() << ';';
        screenshotFile << ind->get_habitat() << ';';
        screenshotFile << ind->get_attack_rates().first << ';';
        screenshotFile << ind->get_attack_rates().second << ';';
        screenshotFile << ind->get_trait_value() << std::endl;
    }

    // Close screenshot file
    screenshotFile.close();

}

/*=======================================================================================================
                                            Main function
========================================================================================================*/

int main(int argc, char * argv[]) {

    try {

        /*=======================================================================================================
                                                Initialization
        ========================================================================================================*/

        // Read in parameters
        if(argc == 1) {
            seed = rnd::set_seed(); // use default parameters values and use clock to set random seed
        }
        else if(argc == 2)
            read_parameters(argv[1]);
        else throw std::runtime_error("invalid number of program arguments in main()");

        // Open log and data files
        std::ostringstream oss;
        oss << "simulation_ebarc_" << seed;
        logFile.open(oss.str() + ".log");
        datFile.open(oss.str() + ".dat");
        if(!(logFile.is_open() && datFile.is_open()))
            throw std::runtime_error("unable to open output file in main()");

        // Write the header of the data file
        std::vector<std::string> headerDataFile = {
                "time",
                "popSize",
                "popSizeH0E0",
                "popSizeH0E1",
                "popSizeH1E0",
                "popSizeH1E1",
                "meanTrait",
                "meanTraitH0E0",
                "meanTraitH0E1",
                "meanTraitH1E0",
                "meanTraitH1E1",
                "resource0H0",
                "resource1H0",
                "resource0H1",
                "resource1H1",
                "spatialIsolation",
                "ecologicalIsolation"
        };
        for(size_t i = 0u; i < headerDataFile.size(); ++i)
            datFile << headerDataFile[i] << ';';
        datFile << std::endl;

        logFile << "Data file created\n";

        // Create a monomorphic population
        std::list<pointer_ind> population; // list of pointers to individuals
        for(size_t i = 0u; i < initialPopSize; ++i) {
            population.push_back(new individual());
        }

        // Initialize time
        size_t timestep = 0u;

        /*=======================================================================================================
                                                    Simulation loop
        ========================================================================================================*/

        logFile << "Simulation begins...\n";

        // Simulation loop
        while(timestep <= Tmax) {

            // Dispersal
            dispersal(population);

            // Resource utilization
            resource_utilization(population);

            // Is it time to get population summary data?
            if(timestep % timeToRecord == 0u)
                record_data(timestep, population);

            // Is it time to take a screenshot of the population?
            if(timestep % timeForScreenshot == 0u)
                take_screenshot(timestep, population);

            // Reproduction
            std::list<pointer_ind> newpopulation = reproduction(population);

            // Death (append survivors to the new population)
            survival(population, newpopulation);

            // Is population extinct?
            if(newpopulation.empty()) {
                logFile << "Population went extinct at t = " << timestep << std::endl;
                break;
            }

            // Erase old population
            population = newpopulation;

            // Advance time
            ++timestep;

        }

        logFile << "Simulation done!" << std::endl;

        /*=======================================================================================================
                                                    Finalization
        ========================================================================================================*/

        // Free allocated memory
        while(!population.empty()) {
            delete population.back();
            population.pop_back();
        }

        // Close output files
        logFile.close();
        datFile.close();

    }

    // Catch my errors
    catch(const std::exception &error) {
        std::cerr << "exception: " << error.what() << '\n';
        logFile << "exception: " << error.what() << '\n';
        exit(EXIT_FAILURE);
    }

    return EXIT_SUCCESS;
}