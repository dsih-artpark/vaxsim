from numba import jit
from scipy.stats import mode

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import logging
from pathlib import Path
from multiprocessing import Pool, cpu_count
from tqdm import tqdm
from pyDOE import lhs
import yaml
import argparse

from discrete_weibull import sirsv_model_with_weibull_calibration


class ABCCalibration:
    def __init__(self, data, param_bounds, baseline_params, n_samples, epsilon, seed):
        print("Initializing ABCCalibration...")

        print(f"Setting up with {n_samples} samples, epsilon={epsilon}, and seed={seed}")
        self.data = data
        self.param_bounds = param_bounds
        self.baseline_params = baseline_params
        self.n_samples = n_samples
        self.epsilon = epsilon
        self.seed = seed

        print(f"Setting random seed to {self.seed}")
        np.random.seed(self.seed)

        print("Preparing observed data...")
        self.prepare_observed_data()

        # Setup logging
        output_folder = Path('output')
        output_folder.mkdir(parents=True, exist_ok=True)
        log_file = output_folder / 'calibration.log'

        logging.basicConfig(
            filename=log_file,
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s'
        )

        # Log CPU and Pool info
        num_cpus = cpu_count()
        logging.info("ABCCalibration initialized with the following parameters:")
        logging.info(f"Number of samples: {self.n_samples}")
        logging.info(f"Epsilon: {self.epsilon}")
        logging.info(f"Seed: {self.seed}")
        logging.info(f"Available CPU cores: {num_cpus}")
        logging.info(f"Using {min(num_cpus, self.n_samples)} cores for multiprocessing")
        logging.info("Preparing observed data...")

        print("Initialization complete.")

    def prepare_observed_data(self):
        print("Converting date column to datetime...")
        self.data['date'] = pd.to_datetime(self.data['date'])

        print("Dropping NaN values and setting date as index...")
        self.observed = self.data.dropna(subset=['sero', 'diva']).set_index('date')

        print("Calculating time in days since 2020-01-01...")
        self.observed['time'] = (self.observed.index - pd.Timestamp('2020-01-01')).days

        print("Extracting sero and diva values...")
        self.observed_sero = self.observed['sero'].values
        self.observed_diva = self.observed['diva'].values

        print(f"Observed data prepared. Shape: {self.observed.shape}")

    @staticmethod
    @jit(nopython=True)
    def calculate_distance_numba(observed_sero, observed_diva, simulated_sero, simulated_diva):
        distance_sero = np.mean((observed_sero - simulated_sero)**2)
        distance_diva = np.mean((observed_diva - simulated_diva)**2)
        return distance_sero + distance_diva

    def sample_prior_lhs(self):
        lhs_samples = lhs(len(self.param_bounds), samples=1)
        sampled_params = {}
        for i, (param, bounds) in enumerate(self.param_bounds.items()):
            sampled_params[param] = bounds[0] + lhs_samples[0, i] * (bounds[1] - bounds[0])
        return {**self.baseline_params, **sampled_params}

    def simulate_and_calculate_distance(self, params):
        simulated_data = sirsv_model_with_weibull_calibration(params, random_seed=self.seed)

        N = sum(params.get(p, 0) for p in ['S0', 'I0', 'R0', 'V0'])
        simulated_sero = (simulated_data['R'] + simulated_data['V']) / (N - simulated_data['I'])
        simulated_diva = simulated_data['R'] / (N - simulated_data['I'])

        min_length = min(len(self.observed_sero), len(simulated_sero))
        distance = self.calculate_distance_numba(
            self.observed_sero[:min_length],
            self.observed_diva[:min_length],
            simulated_sero[:min_length],
            simulated_diva[:min_length]
        )

        return params, distance

    def run_calibration(self):
        accepted_params = []
        accepted_distances = []
        rejected_distances = []
        total_simulations = 0

        with Pool(processes=cpu_count()) as pool:
            pbar = tqdm(total=self.n_samples, desc='ABC Calibration')

            while len(accepted_params) < self.n_samples:
                batch_size = min(100, self.n_samples - len(accepted_params))
                batch_params = [self.sample_prior_lhs() for _ in range(batch_size)]

                results = pool.map(self.simulate_and_calculate_distance, batch_params)

                for params, distance in results:
                    total_simulations += 1  # Increment total simulations immediately

                    # Log each simulation in real-time
                    logging.info(f"Simulation {total_simulations} - Distance: {distance:.6f} "
                                f"- {'Accepted' if distance < self.epsilon else 'Rejected'}")

                    # Check if the simulation is accepted
                    if distance < self.epsilon:
                        accepted_params.append(params)
                        accepted_distances.append(distance)
                        pbar.update(1)  # Update progress bar only for accepted samples
                    else:
                        rejected_distances.append(distance)

                pbar.set_postfix({'Accepted': len(accepted_params), 'Total Simulations': total_simulations})
                logging.info(f"Batch processed: Total simulations = {total_simulations}, "
                            f"Accepted = {len(accepted_params)}")

            pbar.close()

        acceptance_rate = len(accepted_params) / total_simulations
        logging.info(f"Calibration complete. Acceptance Rate: {acceptance_rate * 100:.2f}%")
        logging.info(f"Total Simulations: {total_simulations}")
        logging.info(f"Accepted Samples: {len(accepted_params)}")
        logging.info(f"Error Distribution: {pd.Series(accepted_distances).describe()}")

        return pd.DataFrame(accepted_params), acceptance_rate, total_simulations, accepted_distances, rejected_distances

    def plot_posterior_distributions(self, accepted_params, param_bounds):
        # Plot posterior distributions
        bounded_params = list(param_bounds.keys())
        n_params = len(bounded_params)
        n_cols = 3
        n_rows = (n_params + n_cols - 1) // n_cols

        fig, axes = plt.subplots(n_rows, n_cols, figsize=(15, 5 * n_rows))
        axes = axes.flatten()

        for i, param in enumerate(bounded_params):
            ax = axes[i]
            values = accepted_params[param]
            ax.hist(values, bins=30, edgecolor='black')

            # Find the peak (mode) of the distribution
            peak = mode(values)[0][0]
            ax.axvline(peak, color='r', linestyle='--', label='Peak')
            ax.set_title(param)
            ax.legend()

        for j in range(i + 1, len(axes)):
            fig.delaxes(axes[j])

        plt.tight_layout()
        plt.savefig('output/posterior_distribution.png', dpi=300)
        print("Posterior distribution plots saved.")


def load_params():
    params_path = Path(__file__).parent / 'params.yaml'
    print(f"Loading parameters from {params_path}")
    with params_path.open('r') as f:
        params = yaml.safe_load(f)
    print("Parameters loaded successfully.")
    return params


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='ABC calibration for SIRSV model')
    parser.add_argument('--samples', type=int, default=10000, help='Number of samples')
    parser.add_argument('--epsilon', type=float, default=0.2, help='Acceptance threshold')
    parser.add_argument('--seed', type=int, default=42, help='Random seed for reproducibility')

    args = parser.parse_args()

    print("Loading data...")
    data = pd.read_csv('data.csv', parse_dates=['date'])
    print(f"Data loaded. Shape: {data.shape}")

    print("Loading parameters...")
    param = load_params()
    param_bounds = param['bounds']
    baseline_params = param['baseline']
    print("Parameters loaded.")

    print("Creating ABCCalibration instance...")
    abc = ABCCalibration(data, param_bounds, baseline_params, args.samples, args.epsilon, args.seed)
    print("ABCCalibration instance created.")

    print("Starting calibration...")
    accepted_params, acceptance_rate, total_simulations, accepted_distances, rejected_distances = abc.run_calibration()
    print(f"Calibration completed with acceptance rate: {acceptance_rate * 100:.2f}%")

    print("Saving accepted parameters to CSV...")
    accepted_params.to_csv('output/accepted_params.csv', index=False)
    print("Accepted parameters saved to output/accepted_params.csv.")

    print("Plotting posterior distributions...")
    abc.plot_posterior_distributions(accepted_params, param_bounds)
    print("Posterior distributions plotted.")

    # Plot histograms for accepted and rejected distances
    print("Plotting histograms for accepted and rejected distances...")
    fig, ax = plt.subplots(figsize=(10, 6))

    ax.hist(accepted_distances, bins=30, color='g', alpha=0.7, label='Accepted')
    ax.hist(rejected_distances, bins=30, color='r', alpha=0.7, label='Rejected')
    ax.set_xlabel('Distance')
    ax.set_ylabel('Frequency')
    ax.set_title('Histogram of Distances for Accepted and Rejected Samples')
    ax.legend()

    plt.tight_layout()
    plt.savefig('output/distance_histogram.png', dpi=300)
    print("Distance histogram saved.")
