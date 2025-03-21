import argparse
import logging
import os
import platform
import sys
import warnings
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml
from tqdm import tqdm
from vaxsim.model import sirsv_model_with_weibull_random_vaccination, sirsv_model_with_weibull_targeted_vaccination
from vaxsim.plot import plot_model, plot_parameter_sweep, plot_waning
from vaxsim.utils import analyse_scenarios, compute_total_infections, run_parameter_sweep

logger = logging.getLogger("vaxsim.run")
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
    """Entry point for the command line interface."""
    parser = argparse.ArgumentParser(description="Run Discrete SIRSV model simulations.")

    parser.add_argument("--scenario",
                        choices=["baseline", "scenario_1a", "scenario_1b", "scenario_1c",
                                 "scenario_2a", "scenario_2b", "scenario_2c", "scenario_2d",
                                 "scenario_3a", "scenario_3b", "scenario_3c", "scenario_3d",
                                 "scenario_4a", "scenario_4b", "scenario_4c", "scenario_4d", "parameter_sweep", "run_scenarios"],
                        default="baseline",
                        help="Select the scenario to run")

    parser.add_argument("--model_type", choices=["targeted", "random"], default="random",
                        help="Select the model type to run. Default is 'random'.")

    def parse_seed_infection(value):
        try:
            method, rate = value.split(":") if ":" in value else (value, "0")
            rate = int(rate)
            if method not in ["random", "brute", "event_series", "none"]:
                raise ValueError("Invalid seed method.")
            return method, rate
        except ValueError as err:
            raise argparse.ArgumentTypeError("Seed infection format must be 'method:rate', e.g., 'random:10'") from err

    parser.add_argument("--seed_infection", type=parse_seed_infection, default=("none", 0),
                        help="Select the seed method (random, brute, event_series, none) and rate (integer) for importing external infections, formatted as 'method:rate'. Example: 'random:10'. Default is 'none:0'.")

    args = parser.parse_args()

    log_filename = f"output/logs/sirsv_model_{args.scenario}_{datetime.now().strftime('%Y%m%d_%H%M')}.log"
    os.makedirs("output/logs/", exist_ok=True)
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
        seed_method, seed_rate = args.seed_infection

        # Set the model based on model_type selection
        if "targeted" in args.model_type:
            sirsv_model = sirsv_model_with_weibull_targeted_vaccination
        elif "random" in args.model_type:
            sirsv_model = sirsv_model_with_weibull_random_vaccination
        else:
            raise ValueError("Invalid model type specified.")

        if args.scenario == "parameter_sweep":
            base_params = param['sweep']
            vax_rate_range = np.linspace(0.003, 0.033334, 20)
            vax_period_range = np.arange(30, 360, 30)
            results = run_parameter_sweep(sirsv_model, base_params, 'vax_rate', vax_rate_range, 'vax_period', vax_period_range, model_type=args.model_type)
            plot_parameter_sweep(results, 'vax_rate', 'vax_period', model_type=args.model_type)
            logging.info("Parameter sweep completed. Check the output directory for results.")

        elif args.scenario == "run_scenarios":
            output_dir = 'output/scenario_analysis'
            Path(output_dir).mkdir(parents=True, exist_ok=True)
            analyse_scenarios(sirsv_model, param, output_dir, model_type=args.model_type)
            logging.info(f"Scenario analysis completed. Check the {output_dir} directory for results.")

        else:
            scenario_params = param[args.scenario]
            logging.info(f"Scenario parameters: {scenario_params}")

            if seed_rate > 0 and seed_method != "none":
                scenario_params["seed_rate"] = seed_rate

            if seed_method == "none":
                # Run vaccination models and store outputs
                S, I, R, V = sirsv_model(scenario_params, args.scenario, diagnosis=True, seed_method='none')
                if args.scenario != "baseline" and scenario_params['seed_rate'] == 0 and scenario_params['I0'] == 0:
                    plot_waning(S, I, R, V, scenario_params['days'], scenario=args.scenario, model_type=args.model_type)
                else:
                    plot_model(S, I, R, V, scenario_params['days'], scenario=args.scenario, model_type=args.model_type)

            elif seed_method == "random":
                # Run vaccination models with random seeding
                S, I, R, V = sirsv_model(scenario_params, args.scenario, diagnosis=True, seed_method='random')
                plot_model(S, I, R, V, scenario_params['days'], scenario=args.scenario, model_type=args.model_type)

            elif seed_method == "brute":
                logging.warning("Brute force seed infection will be implemented.")
                output_dir = 'output/bruteforceseeding'
                os.makedirs(output_dir, exist_ok=True)
                total_infections_per_seeded_day = []
                normalized_infections_per_seeded_day = []

                for day in tqdm(range(scenario_params['days']), desc="Running Brute Force Simulations"):
                    event_series = [0] * scenario_params['days']
                    event_series[day] = 1
                    S, I, R, V = sirsv_model(scenario_params, args.scenario, diagnosis=False, seed_method='event_series', event_series=event_series, save_variables=False)
                    total_infections = compute_total_infections(I)
                    total_infections_per_seeded_day.append(total_infections)

                for day in range(len(total_infections_per_seeded_day)):
                    remaining_days = scenario_params['days'] - day
                    normalized_infections = total_infections_per_seeded_day[day] / remaining_days
                    normalized_infections_per_seeded_day.append(normalized_infections)

                plt.figure(figsize=(10, 6))
                plt.plot(range(1, scenario_params['days'] + 1), total_infections_per_seeded_day, marker='o', color='red', label='Total Infections')
                plt.plot(range(1, scenario_params['days'] + 1), normalized_infections_per_seeded_day, marker='o', color='orange', label='Normalized Infections')
                plt.title(f'{args.scenario} Total Infections Based on Remaining Days')
                plt.xlabel('Day of Seeding')
                plt.ylabel('Infections')
                plt.xticks(range(1, scenario_params['days'] + 1), rotation=45)
                plt.grid()
                plt.legend()
                plt.savefig(os.path.join(output_dir, f'{args.scenario}_infections_plot.png'))
                plt.close()

            elif seed_method == "event_series":
                logging.warning("Event series seed infection is an upcoming feature where user can manually define the seed schedule.")

            else:
                logging.info("No seed infection will be imported.")

    except Exception as e:
        logging.error(f"An error occurred: {e}")
        logging.error("Traceback:", exc_info=True)


if __name__ == "__main__":
    main()
