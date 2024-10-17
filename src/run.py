import warnings
from pathlib import Path

import numpy as np
import pandas as pd

import yaml
import argparse
import logging
import platform
import sys
from datetime import datetime

from model import sirsv_model_with_weibull_random_vaccination, sirsv_model_with_weibull_targetted_vaccination
from plot import plot_parameter_sweep, plot_waning, plot_model
from utils import analyse_scenarios, run_parameter_sweep

warnings.filterwarnings('ignore')


def load_params():
    params_path = Path(__file__).parent / 'params.yaml'
    with params_path.open('r') as f:
        return yaml.safe_load(f)


def log_system_info():
    """Logs system and environment details."""
    logging.info(f"System: {platform.system()} {platform.release()}")
    logging.info(f"Architecture: {platform.machine()}")
    logging.info(f"Python version: {sys.version}")
    logging.info(f"Numpy version: {np.__version__}")
    logging.info(f"Pandas version: {pd.__version__}")


def main():
    parser = argparse.ArgumentParser(description="Run Discrete SIRSV model simulations.")

    parser.add_argument("--scenario",
                        choices=["baseline", "scenario_1a", "scenario_1b", "scenario_1c", "scenario_2a", "scenario_2b", "scenario_2c", "scenario_2d", "scenario_3a", "scenario_3b", "scenario_3c", "scenario_3d", "scenario_4a", "scenario_4b", "scenario_4c", "scenario_4d", "parameter_sweep", "run_scenarios"],
                        default="baseline",
                        help="Select the scenario to run")

    parser.add_argument("--model_type", nargs='*', choices=["targetted", "random"], default=["random"],
                        help="Select the model type to run. Multiple choices allowed. Default is 'random'.")

    args = parser.parse_args()

    log_filename = f"output/logs/sirsv_model_{args.scenario}_{datetime.now().strftime('%Y%m%d_%H%M')}.log"
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[logging.FileHandler(log_filename)]
    )

    try:
        logging.info(f"Running simulation with scenario: {args.scenario}")
        logging.info(f"Selected model types: {args.model_type}")
        log_system_info()

        param = load_params()

        if args.scenario == "parameter_sweep":

            for model in args.model_type:
                if model == "targetted":
                    sirsv_model = sirsv_model_with_weibull_targetted_vaccination
                elif model == 'random':
                    sirsv_model = sirsv_model_with_weibull_random_vaccination

            base_params = param['sweep']
            vax_rate_range = np.linspace(0.003, 0.033334, 20)
            vax_period_range = np.arange(30, 360, 30)
            results = run_parameter_sweep(sirsv_model, base_params, 'vax_rate', vax_rate_range, 'vax_period', vax_period_range, model_type=args.model_type)
            plot_parameter_sweep(results, 'vax_rate', 'vax_period', model_type=args.model_type)
            logging.info("Parameter sweep completed. Check the output directory for results.")

        elif args.scenario == "run_scenarios":
            output_dir = 'output/scenario_analysis'
            Path(output_dir).mkdir(parents=True, exist_ok=True)

            for model in args.model_type:
                if model == "targetted":
                    sirsv_model = sirsv_model_with_weibull_targetted_vaccination
                elif model == 'random':
                    sirsv_model = sirsv_model_with_weibull_random_vaccination

            analyse_scenarios(sirsv_model, param, output_dir, model_type=args.model_type)
            logging.info(f"Scenario analysis completed. Check the {output_dir} directory for results.")

        else:
            scenario_params = param[args.scenario]
            logging.info(f"Scenario parameters: {scenario_params}")

            for model in args.model_type:
                if model == "targetted":
                    sirsv_model = sirsv_model_with_weibull_targetted_vaccination
                elif model == 'random':
                    sirsv_model = sirsv_model_with_weibull_random_vaccination

                # Run vaccination models and store outputs
                S, I, R, V = sirsv_model(scenario_params, args.scenario, diagnosis=True)

                # Plotting depending on the scenario
                if args.scenario != "baseline":
                    plot_waning(S, I, R, V, scenario_params['days'], scenario=args.scenario, model_type=model)
                else:
                    plot_model(S, I, R, V, scenario_params['days'], scenario=args.scenario, model_type=model)

    except Exception as e:
        logging.error(f"An error occurred: {e}")
        logging.error("Traceback:", exc_info=True)


if __name__ == "__main__":
    main()
