Simulation and Usage
====================

Model Components
----------------

Core Model
~~~~~~~~~~
.. automodule:: vaxsim.model
   :members:
   :undoc-members:

Calibration
~~~~~~~~~~~
.. automodule:: vaxsim.calibration
   :members:
   :undoc-members:

Utilities
~~~~~~~~~
.. automodule:: vaxsim.utils
   :members:
   :undoc-members:

Visualization
~~~~~~~~~~~~~

The plotting functions are organized into two main categories:

Model Visualization
-------------------
.. autofunction:: vaxsim.plot.plot_model
.. autofunction:: vaxsim.plot.plot_data
.. autofunction:: vaxsim.plot.plot_waning

Analysis and Diagnosis Plots
----------------------------
.. autofunction:: vaxsim.plot.plot_parameter_sweep
.. autofunction:: vaxsim.plot.plot_histogram
.. autofunction:: vaxsim.plot.compare_infections
.. autofunction:: vaxsim.plot.compare_cases_and_infections

All plots are saved in the following directory structure:

.. code-block:: text

    output/
    ├── plots/
    │   ├── scenario_plot_modeltype.png
    │   ├── scenario_waning_modeltype.png
    │   └── infections_comparison_scenario_modeltype.png
    ├── sweep/
    │   └── parameter_sweep_*.png
    └── diagnosis/
        └── scenario/
            └── decay_times_*.png

Usage Examples
--------------

Python Library
~~~~~~~~~~~~~~

.. code-block:: python

    from vaxsim.utils import load_params
    from vaxsim.model import sirsv_model_with_weibull_random

    # Load and select scenario parameters
    params = load_params()['baseline']
    
    # Run simulation
    S, I, R, V = sirsv_model_with_weibull_random(params)

Command Line Interface
~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

    python run_vaxsim.py --scenario baseline --model_type random

Required Arguments:
    * ``--scenario``: Simulation scenario (see :ref:`available-scenarios`)
    * ``--model_type``: Vaccination strategy (see :ref:`vaccination-strategies`)

Configuration
-------------

Model Parameters
~~~~~~~~~~~~~~~~
Parameters are defined in ``params.yaml``:

.. list-table:: Required Parameters
   :header-rows: 1
   :widths: 20 20 60

   * - Parameter
     - Units
     - Description
   * - beta
     - day⁻¹
     - Transmission rate
   * - gamma
     - day⁻¹
     - Recovery rate
   * - vax_rate
     - day⁻¹
     - Vaccination rate
   * - weibull_shape_vax
     - -
     - Vaccine immunity waning shape
   * - weibull_scale_vax
     - days
     - Vaccine immunity waning scale
   * - days
     - days
     - Simulation duration
   * - S0, I0, R0, V0
     - animals
     - Initial population states

Example Configuration
~~~~~~~~~~~~~~~~~~~~~

.. code-block:: yaml

    baseline:
        beta: 0.125          # Transmission rate
        gamma: 0.07          # Recovery rate
        vax_rate: 0.00833    # 25% coverage
        weibull_shape_vax: 3
        weibull_scale_vax: 220
        days: 1095
        S0: 639996
        I0: 4
        R0: 180000
        V0: 180000

.. _available-scenarios:

Available Scenarios
-------------------

.. list-table:: Simulation Scenarios
   :header-rows: 1
   :widths: 20 80

   * - Scenario
     - Description
   * - baseline
     - 2021 FMD outbreak baseline
   * - scenario_1[a-c]
     - Annual vaccination scenarios
   * - scenario_2[a-d]
     - Bi-annual vaccination scenarios
   * - scenario_3[a-d]
     - Tri-annual vaccination scenarios
   * - scenario_4[a-d]
     - Continuous vaccination scenarios
   * - parameter_sweep
     - Parameter sensitivity analysis

.. _vaccination-strategies:

Vaccination Strategies
----------------------

.. list-table::
   :widths: 20 80
   :header-rows: 1

   * - Strategy
     - Description
   * - ``random``
     - Random selection with Weibull immunity waning
   * - ``targeted``
     - Priority-based on immunity levels