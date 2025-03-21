import os
import sys

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from vaxsim.utils import auc_below_threshold

plt.rcParams.update({
    'font.size': 14, 
    'axes.labelsize': 14,
    'axes.titlesize': 16,
    'legend.fontsize': 12,
    'xtick.labelsize': 12,
    'ytick.labelsize': 12
})


def plot_histogram(decay_times_vax, decay_times_rec, scenario, round_counter, start=True):
    """Plot histogram of immunity decay times for vaccinated and recovered populations.

    Parameters
    ----------
    decay_times_vax : array-like
        Decay times for vaccinated population
    decay_times_rec : array-like
        Decay times for recovered population
    scenario : str
        Name of simulation scenario
    round_counter : int
        Current vaccination round number
    start : bool, optional
        If True, plots beginning of round, if False, end of round

    Notes
    -----
    Saves plot to: output/diagnosis/{scenario}/decay_times_{scenario}_round_{round}_[begin|end].png
    """
    output_dir = f'output/diagnosis/{scenario}'
    os.makedirs(output_dir, exist_ok=True)

    plt.figure(figsize=(10, 6))
    plt.hist(decay_times_vax, bins=30, alpha=0.5, label='Vaccinated')
    plt.hist(decay_times_rec, bins=30, alpha=0.5, label='Recovered')
    plt.xlabel('Decay Time')
    plt.ylabel('Frequency')
    plt.title(f'Decay Times for {scenario.capitalize()} - Round {round_counter} {"Beginning" if start else "End"}')
    plt.legend()

    plt.text(0.95, 0.05, f"Vaccinated Count: {len(decay_times_vax)}", transform=plt.gca().transAxes, fontsize=10, horizontalalignment='center', bbox=dict(facecolor='white', alpha=0.5))
    file_path = f'{output_dir}/decay_times_{scenario}_round_{round_counter}_{"begin" if start else "end"}.png'
    plt.savefig(file_path)
    plt.close()


def plot_model(S, I, R, V, days, scenario, model_type, output_dir='output/plots'):
    """Plot SIRSV model simulation results showing protected and recovered fractions.

    Parameters
    ----------
    S, I, R, V : numpy.ndarray
        Time series of population in each compartment
    days : int
        Simulation duration in days
    scenario : str
        Name of simulation scenario
    model_type : str
        Type of vaccination strategy ('random' or 'targeted')
    output_dir : str, optional
        Output directory for plots

    Notes
    -----
    Saves plot to: output/plots
    """
    os.makedirs(output_dir, exist_ok=True)

    data = pd.read_csv('data copy.csv', parse_dates=['date'], index_col='date')
    data['month'] = (data.index - data.index[0]).days / 30

    t = np.arange(days) / 30  # Convert days to months
    N = S + I + R + V

    fig, axs = plt.subplots(1, 2, figsize=(14, 6), sharey=True)

    protected = (R + V) / (N - I)
    axs[0].plot(t, protected, color='#1b9e77', linestyle='-', linewidth=2, label='Protected Fraction (Baseline)')
    axs[0].set_ylim(0, 1)
    axs[0].set_ylabel("Fraction of Population")
    axs[0].set_xlabel("Months since start of simulation")

    axs[1].plot(t, R / (N - I), color='#7570b3', linestyle='-', linewidth=2, label='Recovered Fraction (Baseline)')
    axs[1].set_xlabel("Months since start of simulation")

    if scenario == 'baseline':
        plot_data(axs, data)

    for ax in axs:
        ax.grid(True, linestyle='--', linewidth=0.5, color='#DDDDDD', alpha=0.7)
        ax.set_xlim(0, t[-1])
        ax.legend(loc='upper left', frameon=False)

    # fig.suptitle("Discrete SIRSV Model Simulation", fontsize=16)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'{scenario}_plot_{model_type}.png'), dpi=300, bbox_inches='tight')
    plt.close()


def plot_data(axs, data):
    """Add observed data points to model comparison plots.

    Parameters
    ----------
    axs : list of matplotlib.axes.Axes
        List of two axes for plotting protected and recovered fractions
    data : pandas.DataFrame
        Data frame containing:
        - sero_eff : Seromonitoring effectiveness data
        - diva : DIVA test results
        - month : Time points in months since start

    Notes
    -----
    Adds two types of data points:
    - Harmonised Seromonitoring Data to first subplot
    - Serosurveillance (DIVA) Data to second subplot
    
    Also adds shaded regions between consecutive vaccination rounds
    when seromonitoring data points are within ~1 month of each other.
    """
    # Seromonitoring data
    sero_data = data.dropna(subset=['sero_eff'])
    if not sero_data.empty:
        sero_eff = sero_data['sero_eff'].values
        sero_months = sero_data['month'].values
        axs[0].errorbar(sero_months, sero_eff, fmt='o', color='#d95f02', capsize=6, label='Harmonised Seromonitoring Data', alpha=0.7)

        # Add vertical shading only between consecutive seromonitoring dates (consecutive months)
        for i in range(len(sero_months) - 1):
            # Only add shading if the next month is within ~1 month (allowing a small margin)
            if (sero_months[i+1] - sero_months[i]) <= 1.1:
                # Label only once for the legend
                label = 'Vaccination Round' if i == 0 else None
                axs[0].axvspan(sero_months[i], sero_months[i+1], color='gray', alpha=0.2, label=label)

    # DIVA data
    diva_data = data.dropna(subset=['diva'])
    if not diva_data.empty:
        diva_true = diva_data['diva'].values
        diva_months = diva_data['month'].values
        axs[1].errorbar(diva_months, diva_true, fmt='s', color='#e7298a', capsize=6, label='Serosurveillance (DIVA) Data', alpha=0.7)


def plot_model_v0(S, I, R, V, days, scenario, model_type, trajectories=None, output_dir='output/plots'):
    os.makedirs(output_dir, exist_ok=True)

    data = pd.read_csv('data copy.csv', parse_dates=['date'], index_col='date')
    data['month'] = (data.index - data.index[0]).days / 30

    t = np.arange(days) / 30  # Convert days to months
    N = S + I + R + V

    fig, axs = plt.subplots(4, 1, figsize=(10, 10))

    if trajectories is not None:
        for traj in trajectories:
            St, It, Rt, Vt = traj
            axs[0].plot(t, (Rt + Vt) / (St + It + Rt + Vt - It), 'b-', alpha=0.2)
            axs[1].plot(t, It / (St + It + Rt + Vt), 'r-', alpha=0.2)
            axs[2].plot(t, Vt / (St + It + Rt + Vt), 'g-', alpha=0.2)
            axs[3].plot(t, Rt / (St + It + Rt + Vt), 'g-', alpha=0.2)

    protected = (R + V) / (N - I)
    axs[0].plot(t, protected, 'b-', linewidth=2, label='Protected (Fit)')
    axs[0].set_ylim(0, 1)

    infected_fraction = I / N
    axs[1].plot(t, infected_fraction, 'r-', linewidth=2, label='Infected')
    axs[1].set_ylim(0, max(0.1, min(np.max(infected_fraction) * 1.1, 1)))

    axs[2].plot(t, V / N, 'g-', linewidth=2, label='Vaccinated')
    axs[2].set_ylim(0, 1)

    # axs[3].plot(t, R / N, 'g-', linewidth=2, label='Recovered')
    axs[3].plot(t, R / (N - I), 'c-', linewidth=2, label='Recovered (fit)')
    axs[3].set_ylim(0, 1)

    if scenario == 'baseline':
        plot_data(axs, data)

    # fig.suptitle('Discrete SIRSV Model', fontsize=16)
    fig.text(0.5, 0.04, "Months since start of simulation", ha='center', fontsize=12)
    fig.text(0.04, 0.5, "Fraction of population", va='center', rotation='vertical', fontsize=12)

    for ax in axs:
        ax.grid(True, linestyle='--', alpha=0.7)
        ax.set_xlim(0, t[-1])

    axs[0].legend(loc='upper left')
    axs[1].legend(loc='upper left')
    axs[2].legend(loc='upper left')
    axs[3].legend(loc='upper left')

    plt.tight_layout(rect=[0.05, 0.05, 0.95, 0.95])
    plt.savefig(os.path.join(output_dir, f'{scenario}_plot_{model_type}.png'), dpi=300, bbox_inches='tight')
    plt.close()


def plot_data_v0(ax, data):
    # Seromonitoring data
    sero_data = data.dropna(subset=['sero_eff'])
    if not sero_data.empty:
        sero_eff = sero_data['sero_eff'].values
        sero_months = sero_data['month'].values
        sero_error = sero_eff * 0.1  # 10% error

        ax[0].errorbar(sero_months, sero_eff, yerr=sero_error, fmt='bo', capsize=5, label='Effective Protection\n(Seromonitoring and Vaccination Coverage Data)')

    # Confirmed cases with smoothing
    # inf_data = data.dropna(subset=['inf_obs'])
    # if not inf_data.empty:
    #     inf_true = inf_data['inf_obs'].rolling(window=3, min_periods=1).mean().values
    #     inf_months = inf_data['month'].values

    #     ax2 = ax[1].twinx()
    #     ax2.plot(inf_months, inf_true, 'ko-', markersize=3, label='Smoothed FMD Cases')
    #     ax2.set_ylabel('Number of Cases', color='k')
    #     ax2.tick_params(axis='y', labelcolor='k')

    # DIVA data
    diva_data = data.dropna(subset=['diva'])
    if not diva_data.empty:
        diva_true = diva_data['diva'].values
        diva_months = diva_data['month'].values
        diva_error = diva_true * 0.1  # 10% error

        ax[1].errorbar(diva_months, diva_true, yerr=diva_error, fmt='mo', capsize=5, label='DIVA Data')


def plot_waning(S, I, R, V, days, scenario, model_type, output_dir='output/plots', herd_threshold=0.416):
    """Plot immunity waning analysis showing protected fraction and vulnerability regions.

    Parameters
    ----------
    S, I, R, V : numpy.ndarray
        Time series of population in each compartment
    days : int
        Simulation duration in days
    scenario : str
        Name of simulation scenario
    model_type : str
        Type of vaccination strategy
    output_dir : str, optional
        Output directory for plots
    herd_threshold : float, optional
        Herd immunity threshold (default: 0.416)

    Notes
    -----
    Saves plot to: output/plots
    Shows regions where protection falls below herd immunity threshold.
    """
    os.makedirs(output_dir, exist_ok=True)

    t = np.arange(days) / 30  # Convert days to months
    N = S + I + R + V

    plt.figure(figsize=(10, 6))

    protected = (R + V) / (N - I)

    # area under the curve
    auc = auc_below_threshold(S, I, R, V, days, herd_threshold)

    plt.plot(t, protected, 'b', label='Protected')

    plt.axhline(y=herd_threshold, color='r', linestyle='--', label='Herd Immunity Threshold')

    plt.fill_between(t, protected, herd_threshold, where=(protected < herd_threshold),
                     color='red', alpha=0.3, interpolate=True, label='Region of vulnerability')

    plt.text(0.05, 0.05, f'Cumulative Vulnerability: {auc:.4f}',
             transform=plt.gca().transAxes, fontsize=12, bbox=dict(facecolor='white', alpha=0.5))
    # plt.title('SIRSV Model (Immunity waning)', fontsize=16)
    plt.xlabel("Months since start of simulation", fontsize=12)
    plt.ylabel("Fraction of population", fontsize=12)
    plt.ylim(0, 1)
    plt.legend(loc='upper right')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'{scenario}_waning_{model_type}.png'), dpi=300, bbox_inches='tight')
    plt.close()


def plot_parameter_sweep(results, param1_name, param2_name, output_variable='protected', 
                        vaccine_efficacy=1, herd_threshold=0.416, model_type='random'):
    """Plot parameter sweep results as a heatmap with threshold indicators.

    Parameters
    ----------
    results : list of dict
        List of simulation results with parameter combinations
    param1_name, param2_name : str
        Names of parameters being swept
    output_variable : str, optional
        Variable to plot (default: 'protected')
    vaccine_efficacy : float, optional
        Vaccine efficacy factor (default: 1)
    herd_threshold : float, optional
        Herd immunity threshold (default: 0.416)
    model_type : str, optional
        Type of vaccination strategy (default: 'random')

    Notes
    -----
    Saves plot to: output/sweep/parameter_sweep_{param1}_{param2}_{output_variable}_{efficacy}_{model_type}.png
    """
    param1_values = sorted(set(result[param1_name] for result in results if result is not None))
    param2_values = sorted(set(result[param2_name] for result in results if result is not None))

    output_grid = np.full((len(param1_values), len(param2_values)), np.nan)
    for result in results:
        if result is None:
            continue

        i = param1_values.index(result[param1_name])
        j = param2_values.index(result[param2_name])
        output_grid[i, j] = result[output_variable] * vaccine_efficacy

    plt.figure(figsize=(12, 10))
    vmin, vmax = np.nanmin(output_grid), np.nanmax(output_grid)
    im = plt.imshow(output_grid, origin='lower', aspect='auto',
                    extent=[min(param2_values), max(param2_values),
                            min(param1_values), max(param1_values)],
                    vmin=vmin, vmax=vmax, cmap='viridis')

    below_threshold = output_grid < herd_threshold
    x_values, y_values = np.meshgrid(param2_values, param1_values)
    plt.scatter(x_values[below_threshold], y_values[below_threshold], color='red', marker='x', label='Below Herd Threshold')

    cbar_label = 'Min protected fraction'
    cbar = plt.colorbar(im)
    cbar.set_label(cbar_label)
    plt.xlabel(f'{param2_name} units')
    plt.ylabel(f'{param1_name} units')
    plt.title(f'Minimum {output_variable} for different {param1_name} and {param2_name}')
    plt.legend()
    plt.savefig(f'output/sweep/parameter_sweep_{param1_name}_{param2_name}_{output_variable}_{vaccine_efficacy}_{model_type}.png', dpi=300, bbox_inches='tight')
    plt.close()


def compare_infections(scenario, model_type='random', output_dir='output/plots'):
    """Compare infection dynamics between baseline and scenario simulations.

    Parameters
    ----------
    scenario : str
        Name of scenario to compare with baseline
    model_type : str, optional
        Type of vaccination strategy (default: 'random')
    output_dir : str, optional
        Output directory for plots

    Notes
    -----
    Saves plot to: {output_dir}/infections_comparison_{scenario}_{model_type}.png
    Loads data from: output/saved_variables/{model_type}_vaccination/
    """
    input_dir = 'output/saved_variables'
    os.makedirs(output_dir, exist_ok=True)
    baseline_inf = np.load(f'{input_dir}/{model_type}_vaccination/baseline/baseline_simulation_results.npz')
    scenario_inf = np.load(f'{input_dir}/{model_type}_vaccination/{scenario}/{scenario}_simulation_results.npz')

    # Assume simulation time is in days; convert to dates starting 2020-01-01
    num_days = baseline_inf['I'].shape[0]
    start_date = pd.to_datetime("2020-01-01")
    dates = start_date + pd.to_timedelta(np.arange(num_days), unit='D')

    plt.figure(figsize=(16, 6))
    plt.plot(dates, baseline_inf['I'], label='Baseline', color='red', linestyle='-')
    plt.plot(dates, scenario_inf['I'], label='Bi-annual', color='orange', linestyle='--')

    plt.xlabel('Date')
    plt.ylabel('Number of infected individuals')
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.legend(loc='upper left')

    # Set x-axis ticks to every 3 months and rotate labels 45°.
    ax = plt.gca()
    ax.xaxis.set_major_locator(mdates.MonthLocator(interval=3))
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m'))
    plt.xticks(rotation=45)
    plt.savefig(f'{output_dir}/infections_comparison_{scenario}_{model_type}.png', dpi=300, bbox_inches='tight')
    plt.close()


def compare_cases_and_infections(scenario, model_type='random', output_dir='output/plots'):
    """Compare observed FMD cases with simulated infections.

    Parameters
    ----------
    scenario : str
        Name of scenario for comparison
    model_type : str, optional
        Type of vaccination strategy (default: 'random')
    output_dir : str, optional
        Output directory for plots

    Notes
    -----
    Saves plot to: {output_dir}/cases_and_infections_{scenario}_{model_type}.png
    Loads data from: 
    - data copy.csv (observed cases)
    - output/saved_variables/{model_type}_vaccination/baseline/
    """
    # Load observed data from CSV
    data = pd.read_csv('data copy.csv', parse_dates=['date'], index_col='date')
    data['month'] = (data.index - data.index[0]).days / 30.0
    inf_data = data.dropna(subset=['inf_obs'])
    if not inf_data.empty:
        # Unsmoothed observed values as dots
        obs_values = inf_data['inf_obs'].values
        inf_months = inf_data['month'].values
        # Convert month values to dates (approximate by multiplying by 30 days)
        start_date = pd.to_datetime("2020-01-01")
        obs_dates = start_date + pd.to_timedelta(inf_months * 30, unit='D')
        # Smooth confirmed cases with a rolling window (window=3, min_periods=1)
        smoothed = inf_data['inf_obs'].rolling(window=3, min_periods=1).mean().values
    else:
        obs_dates = None
        obs_values = None
        smoothed = None

    # Load baseline simulation results
    input_dir = 'output/saved_variables'
    baseline_inf = np.load(f'{input_dir}/{model_type}_vaccination/baseline/baseline_simulation_results.npz')
    num_days = baseline_inf['I'].shape[0]
    start_date = pd.to_datetime("2020-01-01")
    sim_dates = start_date + pd.to_timedelta(np.arange(num_days), unit='D')

    # Create figure with two stacked panels
    fig, axs = plt.subplots(2, 1, figsize=(16, 10), sharex=True)

    # Top panel: Observed FMD Cases (dots and smoothed line)
    if obs_values is not None:
        # Plot the unsmoothed observed data as dots
        axs[0].plot(obs_dates, obs_values, 'ko', markersize=4, label='Observed Cases')
        # Plot the smoothed line on top of the dots
        axs[0].plot(obs_dates, smoothed, 'k-', linewidth=2, label='Smoothed Cases (window=3)')
        axs[0].set_ylabel('Number of Cases')
        axs[0].tick_params(axis='y')
        axs[0].legend(loc='upper left')
        axs[0].grid(True, linestyle='--', alpha=0.7)
    else:
        axs[0].text(0.5, 0.5, "No observed case data", transform=axs[0].transAxes,
                    ha='center', fontsize=14)

    # Bottom panel: Baseline Simulation Infections
    axs[1].plot(sim_dates, baseline_inf['I'], label='Baseline',
                color='red', linestyle='-')
    axs[1].set_ylabel('Number of Infected Individuals')
    axs[1].tick_params(axis='y')
    axs[1].legend(loc='upper left')
    axs[1].grid(True, linestyle='--', alpha=0.7)

    # X-axis formatting: major ticks every 3 months, 45° rotation, date formatter.
    axs[1].set_xlabel('Date')
    axs[1].xaxis.set_major_locator(mdates.MonthLocator(interval=3))
    axs[1].xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m'))
    plt.xticks(rotation=45)

    # plt.suptitle(f"Observed Cases and Baseline Infections for {scenario}", fontsize=16)
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])

    os.makedirs(output_dir, exist_ok=True)
    plt.savefig(f'{output_dir}/cases_and_infections_{scenario}_{model_type}.png', dpi=300, bbox_inches='tight')
    plt.close()


if __name__ == "__main__":
    # Usage: python plot.py <scenario> [<model_type>]
    scenario = sys.argv[1] if len(sys.argv) > 1 else 'scenario_2b'
    model_type = sys.argv[2] if len(sys.argv) > 2 else 'random'
    compare_infections(scenario, model_type=model_type)
    # compare_cases_and_infections(scenario, model_type=model_type)
