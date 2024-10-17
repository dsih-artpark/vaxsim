## *VaxSim*: Scenario Analysis tool for Vaccination Strategies to control Foot and Mouth Disease (FMD) transmission

---------------------------------------------------------------

[![MIT license](https://img.shields.io/badge/License-MIT-blue.svg)](https://github.com/epiforecasts/EpiNow2/blob/main/LICENSE.md/)

*VaxSim* is a python package for simulating the spread of Foot and Mouth Disease (FMD) using a discrete compartmental model (SIRSV framework), along with control interventions such as periodic pulsed vaccination campaigns.

Installation
-------------
You can install vaxsim using [Poetry](https://python-poetry.org/) by running the following command:

```bash
poetry add git+https://github.com/dsih-artpark/vaxsim.git
```
-------------

### Model Parameters

| **Parameter** | **Description**                                                      | **Baseline** | **Min**  | **Max**  | **Unit**         | **Source**                        |
|---------------|----------------------------------------------------------------------|--------------|----------|----------|------------------|-----------------------------------|
| $\beta$       | Transmission rate                                                   | 0.12         | 1e-06    | -        | days<sup>-1</sup>| model fit                        |
| $\gamma$      | Recovery rate                                                       | 0.07         | 0.065    | 0.122    | days<sup>-1</sup>| [Yadav et al. (2019)](https://doi.org/10.3389/fvets.2019.00263) |
| $\mu$         | Effective vaccination rate                                          | 0.016        | 1e-06    | 0.3333   | days<sup>-1</sup>| model fit                        |
| $\lambda_{\text{vax}}$ | Scale parameter for immunity decay in vaccinated individuals | 160          | 120      | 180      | days             | [Singh et al. (2019)](https://doi.org/10.3390/vaccines7030090)  |
| $k_{\text{vax}}$ | Shape parameter for immunity decay in vaccinated individuals      | 3            | 1        | 6        | -                | model fit                        |
| $\lambda_{\text{rec}}$ | Scale parameter for immunity decay in recovered individuals  | 1380         | 365      | 1800     | days             | [Pomeroy et al. (2015)](https://doi.org/10.1371/journal.pone.0136642) |
| $k_{\text{rec}}$ | Shape parameter for immunity decay in recovered individuals       | 3            | 1        | 6        | -                | model fit                        |
| $N$           | Total population                                                    | 1000000      | -        | -        | -                | -                                 |
| $S_0$         | Initial susceptible population                                       | 699998       | 400000   | 800000   | -                | model fit                        |
| $I_0$         | Initial infected population                                         | 2            | 1        | 10       | -                | model fit                        |
| $R_0$         | Initial recovered population                                         | 150000       | 100000   | 300000   | -                | model fit                        |
| $V_0$         | Initial vaccinated population                                        | 150000       | 100000   | 300000   | -                | model fit                        |
| days          | Simulation duration                                                 | 1100         | -        | -        | days             | -                                 |
| $\omega$      | Duration between successive pulses of vaccination                    | 335          | 1        | 1100     | days             | vaccination schedule datasets (GoK)     |
| $\tau$        | Duration of each vaccination pulse                                   | 30           | 1        | 1100     | days             | vaccination progress datasets (GoK)     |


### Running model scenarios

Run the script from the command line as follows:

```python3
python run.py --scenario <scenario_name> --model_type <model_type>
```

To run the baseline scenario with default model

```python3
python run.py
```

List of scenarios are given below:
| Scenario Name | Remarks |
|----------------|---------|
| baseline       | 2021 FMD outbreak |
| scenario_1a    | Disease-free state with 100% annual vaccination starts from Day 1. |
| scenario_1b    | Disease-endemic state with 100% annual vaccination; Vacc starts from Day 1. |
| scenario_1c    | Baseline vaccination strategy in absence of infections. |
| scenario_2a    | Disease-free state with 100% bi-annual vaccination; Vacc starts from Day 1. |
| scenario_2b    | Disease-endemic state with 100% bi-annual vaccination; Vacc starts from Day 1. |
| scenario_2c    | Disease-free state with 50% bi-annual vaccination; Vacc starts from Day 1. |
| scenario_2d    | Disease-free state with 25% bi-annual vaccination; Vacc starts from Day 1. |
| scenario_3a    | Disease-free state with tri-annual vaccination and 66% coverage; Vacc starts from Day 1. |
| scenario_3b    | Disease-endemic state with tri-annual vaccination and 66% coverage; Vacc starts from Day 1. |
| scenario_3c    | Disease-free state with tri-annual vaccination and 33% coverage; Vacc starts from Day 1. |
| scenario_3d    | Disease-free state with tri-annual vaccination and 16.6% coverage; Vacc starts from Day 1. |
| scenario_4a    | Disease-free state with continuous vaccination and 2/365 daily coverage; Vacc starts from Day 1. |
| scenario_4b    | Disease-endemic state with continuous vaccination and 2/365 daily coverage; Vacc starts from Day 1. |
| scenario_4c    | Disease-free state with continuous vaccination and 1/365 daily coverage; Vacc starts from Day 1. |
| scenario_4d    | Disease-free state with continuous vaccination and 0.5/365 daily coverage; Vacc starts from Day 1. |
| sweep          | Parameter sweep |


### Output
The results of the simulations, including log files and visualisations, will be saved in the *`output`*  directory. Check the relevant subdirectories for logs and plots corresponding to each simulation run.

### Logging
Logs are created in the *`output/logs`* directory with filenames indicating the scenario and timestamp. The logging level is set to INFO by default.

### Contributing
Contributions to improve the model, enhance functionalities, or fix issues are welcome. Please fork the repository and submit a pull request with your changes.

### Licence
This project is licensed under the MIT Licence - see the LICENCE file for details.
