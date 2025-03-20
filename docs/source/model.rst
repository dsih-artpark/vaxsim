Model Description
=================

Introduction
------------
The SIRSV model is an epidemiological framework designed to simulate the spread of infectious diseases, with a particular focus on Foot and Mouth Disease (FMD). This model expands upon the classical SIR (Susceptible-Infectious-Recovered) framework by introducing a fourth compartment for vaccinated individuals (V), thereby considering the effects of vaccination and immunity waning on disease transmission dynamics. Unlike traditional compartmental models, this approach allows for the re-vaccination of previously vaccinated individuals during each round, while systematically tracking the immunity decay time of each individual in vaccinated state.

The model consists of four compartments:

* S: Susceptible population
* I: Infected population
* R: Recovered population
* V: Vaccinated population

Immunity Waning
---------------
Weibull distribution is used to generate decay times for immunity in the vaccinated population. This distribution is characterised by two parameters: shape (k) and scale (lambda). The shape parameter determines the rate of decay, while the scale parameter determines the overall duration of immunity. By incorporating this distribution, the model can simulate the gradual loss of immunity over time, reflecting the real-world scenario of immunity waning post-vaccination.

Re-vaccination Strategies
-------------------
1. Random Vaccination
   * Random selection of animals for vaccination
   * Uniform coverage across population

2. Targeted Vaccination
   * Priority based on immunity decay time
   * Targets animals with lowest immunity