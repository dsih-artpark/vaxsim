"""
Calibration module using Sequential Monte Carlo ABC sampling.
"""

import time
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import gaussian_kde

from vaxsim.model import sirsv_model_with_weibull_calibration
from vaxsim.utils import model_loss, load_params


def log_results(params, loss, iteration, log_file):
    """Log calibration results to a CSV file."""
    log_data = {'iteration': iteration, 'loss': loss, **params}
    log_df = pd.DataFrame([log_data])
    log_df.to_csv(log_file, mode='a', header=not log_file.exists(), index=False)

def loss_function(params, bounds_keys, baseline, data):
    """
    Compute loss for given parameters.

    Parameters
    ----------
    params : array-like
        Parameter values (including scale_diva).
    bounds_keys : list
        List of parameter names.
    baseline : dict
        Baseline parameter values.
    data : pd.DataFrame
        Input data with seromonitoring and diva columns.

    Returns
    -------
    float
       Computed loss value.
    """
    # Reconstruct parameter dictionary and extract scale_diva
    param_dict = {key: val for key, val in zip(bounds_keys, params)}
    scale_diva = param_dict.pop('scale_diva', 0.5)
    param_dict = {**baseline, **param_dict}

    S, I, R, V = sirsv_model_with_weibull_calibration(param_dict)
    loss = model_loss(S, I, R, V, data, scale_diva)
    print(f"Loss: {loss:.4f}")
    return loss

def plot_parameter_distributions(samples, param_names, output_dir):
    """
    Plot KDE distributions for each parameter with the mode indicated by a vertical dashed line.
    
    The mode value is annotated with 90° rotated text, positioned inside the plot area.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    for i, param in enumerate(param_names):
        plt.figure(figsize=(6, 4))
        kde = gaussian_kde(samples[:, i])
        xs = np.linspace(samples[:, i].min(), samples[:, i].max(), 1000)
        ys = kde(xs)
        mode_value = xs[np.argmax(ys)]
        
        sns.kdeplot(samples[:, i], fill=True, color='lightblue', linewidth=2)
        plt.axvline(mode_value, color='r', linestyle='dashed', linewidth=1)
        
        # Get current axis limits for better placement
        ax = plt.gca()
        ylim = ax.get_ylim()[1]
        # Place text with a slight horizontal offset and inside the plotting area
        plt.text(mode_value, ylim * 0.95, f'{mode_value:.2f}', color='r',
                 ha='center', va='top', fontsize=10, rotation=90, alpha=0.7,
                 bbox=dict(facecolor='white', edgecolor='none', pad=1.0, alpha=0.5))
        
        plt.title(f"Distribution for {param}")
        plt.xlabel(param)
        plt.ylabel("Density")
        plt.tight_layout()
        plt.savefig(output_dir / f"{param}_distribution.png")
        plt.close()

def smc_abc_sampling(num_particles=200, num_generations=5, initial_epsilon=1.0, final_epsilon=0.1):
    """
    Perform Sequential Monte Carlo ABC sampling with profiling.

    Parameters
    ----------
    num_particles : int
        Number of particles per generation.
    num_generations : int
        Total number of generations.
    initial_epsilon : float
        Epsilon threshold for the first generation.
    final_epsilon : float
        Epsilon threshold for the final generation.

    Returns
    -------
    np.ndarray
        Final set of accepted particles.
    """
    overall_start = time.time()
    params_dict = load_params()
    bounds_dict = params_dict['bounds']
    baseline = params_dict['baseline']
    bounds_keys = list(bounds_dict.keys())
    bounds = np.array([bounds_dict[key] for key in bounds_keys])

    data_path = Path(__file__).parent.parent.parent / 'data copy.csv'
    data = pd.read_csv(data_path, parse_dates=['date'], index_col='date')

    # Create a unique run identifier based on timestamp and input settings
    run_id = f"Run_{time.strftime('%Y%m%d_%H%M%S')}_P{num_particles}_G{num_generations}_E{initial_epsilon}-{final_epsilon}"
    # Use run_id for log file and plots folder
    log_file = Path(__file__).parent.parent.parent / 'output' / 'calibration' / f"calibration_log_{run_id}.csv"
    plots_dir = Path(__file__).parent.parent.parent / 'output' / 'calibration' / "plots" / run_id

    epsilons = np.linspace(initial_epsilon, final_epsilon, num_generations)

    # Generation 0: Uniform sampling within bounds
    particles = []
    for i in range(num_particles):
        sampled_params = np.random.uniform(bounds[:, 0], bounds[:, 1])
        loss = loss_function(sampled_params, bounds_keys, baseline, data)
        if loss <= epsilons[0]:
            particles.append(sampled_params)
            log_results(dict(zip(bounds_keys, sampled_params)), loss, iteration=f"Gen0_{i}", log_file=log_file)
            print(f"Gen0_{i}: Loss = {loss:.4f}")
    particles = np.array(particles)
    print(f"Generation 0: Accepted {particles.shape[0]} out of {num_particles} particles.")

    # Sequential generations: perturb accepted particles
    for gen in range(1, num_generations):
        gen_start = time.time()
        new_particles = []
        current_epsilon = epsilons[gen]
        print(f"\nStarting Generation {gen} with epsilon = {current_epsilon}")

        if particles.shape[0] == 0:
            print("No accepted particles in previous generation; terminating SMC ABC sampling.")
            break

        cov = np.cov(particles, rowvar=False) + 1e-6 * np.eye(particles.shape[1])
        iteration_counter = 0
        while len(new_particles) < num_particles:
            idx = np.random.choice(len(particles))
            base_particle = particles[idx]
            perturbed = np.random.multivariate_normal(base_particle, cov)
            if not all(bounds[j, 0] <= perturbed[j] <= bounds[j, 1] for j in range(len(bounds_keys))):
                continue
            loss = loss_function(perturbed, bounds_keys, baseline, data)
            iteration_id = f"Gen{gen}_{iteration_counter}"
            print(f"{iteration_id}: Loss = {loss:.4f}")
            if loss <= current_epsilon:
                new_particles.append(perturbed)
                log_results(dict(zip(bounds_keys, perturbed)), loss, iteration=iteration_id, log_file=log_file)
            iteration_counter += 1
        particles = np.array(new_particles)
        print(f"Generation {gen}: Accepted {particles.shape[0]} particles in {time.time() - gen_start:.2f} seconds")

    overall_time = time.time() - overall_start
    acceptance_rate = particles.shape[0] / (num_particles * num_generations)
    print(f"\nFinal acceptance rate: {acceptance_rate:.4f}")
    print(f"Total SMC ABC sampling time: {overall_time:.2f} seconds")

    plot_parameter_distributions(particles, bounds_keys, plots_dir)

    return particles

if __name__ == "__main__":
    final_particles = smc_abc_sampling(num_particles=10, num_generations=10,
                                       initial_epsilon=0.7, final_epsilon=0.07)
