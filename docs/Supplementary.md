### Model Fit

Parameters are optimised using the 'SLSQP' method (Han-Powell quasi-Newton method) from the `minimize` function in `SciPy`, incorporating bounds and constraints to maintain the total sum of the fraction of the population within each compartment—susceptible (S), infected (I), recovered (R), and vaccinated (V)—summing to 1. The solver settings include a relative tolerance (`rtol`) and absolute tolerance (`atol`), both set to \(10^{-8}\), ensuring precision, and a maximum step size (`max_step`) of 0.1 to control integration granularity.

#### Loss Function

$$\sum_{t \in T_{P}^{\star}} \biggl(\frac{R(t)+V(t)}{1-I(t)} - F_{p}^{\star}(t)\biggr)^2 + \sum_{t \in T_{D}^{\star}} \biggl(\frac{R(t)}{1-I(t)} - F_{D}^{\star}(t)\biggr)^2$$
where $T_{P}$ is the set of time points where percent protected data is available, $T_{D}$ is the set of time points where DIVA data is available, $F_{P}^{\star}(t)$ is the fraction of SP(structural protein) protected animals and $F_{D}^{\star}(t)$ is the fraction of NSP(non-structural protein) DIVA positive animals.


## Model Assumptions and Limitations

- The model assumes that individuals within the population mix homogeneously. The total population size remains constant throughout the simulation, and individuals only transition between compartments without any changes in the overall population size.
- Vaccination is assumed to have an immediate effect, with individuals in the Susceptible compartment moving directly to the Vaccinated compartment without any delay. The model does not account for the possibility of partial immunity or delayed onset of immunity following vaccination.
- Immunity is assumed to wane following a sigmoid decay, where immunity duration decreases over time based on decay times drawn from a Weibull distribution. This models the variability in immunity loss, with individuals transitioning from protected to susceptible compartments as their decay time reaches zero.
- The estimation of \(R_0\) may be affected by the spatial difference in the onset of the outbreak wave. Such differences can influence local transmission rates and the timing of interventions.
- The seromonitoring data is only available for cattle and buffalo; other species are ignored in the analysis. Future iterations of the model may incorporate serological data from additional species as it becomes available.
- We assume that animals showing antibody titres above the threshold will not be susceptible to infection during the wave. However, this assumption may or may not be valid in the field. The effectiveness of vaccination against FMD is influenced by many factors, including vaccine quality (potency and antigenic relevance), the way vaccination is implemented (e.g., regime, cold chain and coverage), the weight of infection that must be blocked (e.g., livestock densities and contact structures), and how well the vaccination is supported by other control measures (e.g., movement controls, biosecurity).
