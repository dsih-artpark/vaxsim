Simulation and Usage
====================

Model Components
----------------

Core Functions
~~~~~~~~~~~~~~
.. autofunction:: vaxsim.model.sirsv_model_with_weibull_random
.. autofunction:: vaxsim.model.sirsv_model_with_weibull_targeted

Calibration Tools
~~~~~~~~~~~~~~~~~
.. autofunction:: vaxsim.calibration.loss_function
.. autofunction:: vaxsim.calibration.calibrate_model

Helper Functions
~~~~~~~~~~~~~~~~
.. automodule:: vaxsim.utils
   :members:
   :undoc-members:

Infection Seeding
-----------------

Seeding Functions
~~~~~~~~~~~~~~~~~
.. autofunction:: vaxsim.utils.generate_seeds
.. autofunction:: vaxsim.utils.compute_total_infections

Parameters:
    * seed_size: Number of infected animals
    * seed_interval: Days between seeding events
    * seed_threshold: Infection threshold for reseeding

Outputs & Visualization
-----------------------

Plot Generation
~~~~~~~~~~~~~~~
.. autofunction:: vaxsim.utils.plot_simulation
.. autofunction:: vaxsim.utils.plot_comparison

Data Management
~~~~~~~~~~~~~~~
.. autofunction:: vaxsim.utils.save_results
.. autofunction:: vaxsim.utils.load_results

Analysis Functions
~~~~~~~~~~~~~~~~~~
.. autofunction:: vaxsim.utils.compute_protection
.. autofunction:: vaxsim.utils.compute_total_infections


Examples and Usage
-------------------

.. note:: Launch an IPython shell (by executing `ipython` in your terminal) to interactively run the code.

Loading parameters from YAML file::

    from vaxsim.utils import load_params
    
    parameters = load_params()

Running the baseline simulation::

    from vaxsim.model import sirsv_model_with_weibull_random
    
    params = load_params()['baseline']
    
    S, I, R, V = sirsv_model_with_weibull_random(params)


.. note:: Command line interface::

    python run_vaxsim.py --scenario <scenario_name> --model_type <model_type>

Available options:

* Scenarios

List of provided scenarios are given below:

| Scenario Name | Remarks |
|---------------|---------|
| baseline      | 2021 FMD outbreak |
| scenario_1a   | Disease-free state with 100% annual vaccination starts from Day 1. |
| scenario_1b   | Disease-endemic state with 100% annual vaccination; Vacc starts from Day 1. |
| scenario_1c   | Baseline vaccination strategy in absence of infections. |
| scenario_2a   | Disease-free state with 100% bi-annual vaccination; Vacc starts from Day 1. |
| scenario_2b   | Disease-endemic state with 100% bi-annual vaccination; Vacc starts from Day 1. |
| scenario_2c   | Disease-free state with 50% bi-annual vaccination; Vacc starts from Day 1. |
| scenario_2d   | Disease-free state with 25% bi-annual vaccination; Vacc starts from Day 1. |
| scenario_3a   | Disease-free state with tri-annual vaccination and 66% coverage; Vacc starts from Day 1. |
| scenario_3b   | Disease-endemic state with tri-annual vaccination and 66% coverage; Vacc starts from Day 1. |
| scenario_3c   | Disease-free state with tri-annual vaccination and 33% coverage; Vacc starts from Day 1. |
| scenario_3d   | Disease-free state with tri-annual vaccination and 16.6% coverage; Vacc starts from Day 1. |
| scenario_4a   | Disease-free state with continuous vaccination and 2/365 daily coverage; Vacc starts from Day 1. |
| scenario_4b   | Disease-endemic state with continuous vaccination and 2/365 daily coverage; Vacc starts from Day 1. |
| scenario_4c   | Disease-free state with continuous vaccination and 1/365 daily coverage; Vacc starts from Day 1. |
| scenario_4d   | Disease-free state with continuous vaccination and 0.5/365 daily coverage; Vacc starts from Day 1. |
| parameter_sweep | The parameter sweep evaluates combinations of vaccination rate and period, visualised in a heatmap where colour represents the output analysis function value for each pair. |

* Model types: random, targeted

Example::

    python run_vaxsim.py --scenario baseline --model_type random

Parameter Configuration
------------------------

Define parameters in params.yaml::

    baseline:
        Remarks: 2021 FMD outbreak
        beta: 0.125
        gamma: 0.07
        vax_rate: 0.00833  # 25% coverage
        weibull_shape_vax: 3
        weibull_scale_vax: 220
        weibull_shape_rec: 3
        weibull_scale_rec: 1380
        days: 1095
        seed_rate: 0
        vax_period: 335
        vax_duration: 30
        start_vax_day: 330
        S0: 639996
        I0: 4
        R0: 180000
        V0: 180000

Additional parameter values can be added to the YAML file as needed in the same format.