import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from vaxsim.utils import auc_below_threshold


def plot_histogram(decay_times_vax, decay_times_rec, scenario, round_counter, start=True):
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


def plot_model(S, I, R, V, days, scenario, model_type, trajectories=None, output_dir='output/plots'):
    os.makedirs(output_dir, exist_ok=True)

    data = pd.read_csv('data.csv', parse_dates=['date'], index_col='date')
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
    axs[1].set_ylim(0, min(np.max(infected_fraction) * 1.1, 1))
    axs[1].set_ylabel('Number of Confirmed Cases', color='k')

    axs[2].plot(t, V / N, 'g-', linewidth=2, label='Vaccinated')
    axs[2].set_ylim(0, 1)

    axs[3].plot(t, R / N, 'g-', linewidth=2, label='Recovered')
    axs[3].plot(t, R / (N - I), 'c-', linewidth=2, label='Recovered (DIVA)')
    axs[3].set_ylim(0, 1)

    plot_data(axs, data)

    fig.suptitle('Discrete SIRSV Model', fontsize=16)
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


def plot_data(ax, data):
    # Seromonitoring data
    sero_data = data.dropna(subset=['sero'])
    sero_true = sero_data['sero'].values
    sero_months = sero_data['month'].values
    sero_error = sero_true * 0.1  # 10% error

    ax[0].errorbar(sero_months, sero_true, yerr=sero_error, fmt='bo', capsize=5, label='Seromonitoring Data')

    # Lab confirmed cases
    inf_data = data.dropna(subset=['inf_obs'])
    inf_true = inf_data['inf_obs'].values
    inf_months = inf_data['month'].values

    ax2 = ax[1].twinx()
    ax2.plot(inf_months, inf_true, 'ko', markersize=3, label='Confirmed cases')
    ax2.tick_params(axis='y', labelcolor='k')

    # DIVA data
    diva_data = data.dropna(subset=['diva'])
    diva_true = diva_data['diva'].values
    diva_months = diva_data['month'].values
    diva_error = diva_true * 0.1  # 10% error

    ax[3].errorbar(diva_months, diva_true, yerr=diva_error, fmt='mo', capsize=5, label='DIVA Data')
    ax2.legend(loc='upper right')


def plot_waning(S, I, R, V, days, scenario, model_type, output_dir='output/plots', herd_threshold=0.416):
    os.makedirs(output_dir, exist_ok=True)

    t = np.arange(days) / 30  # Convert days to months
    N = S + I + R + V

    plt.figure(figsize=(10, 6))

    protected = (R + V) / (N - I)

    # area under the curve
    auc = auc_below_threshold(S, I, R, V, days, herd_threshold)

    plt.plot(t, protected, 'b', label='Protected')

    plt.axhline(y=herd_threshold, color='r', linestyle='--', label='Herd Immunity Threshold')

    plt.fill_between(t, protected, 1, where=(protected < herd_threshold),
                     color='red', alpha=0.3, interpolate=True, label='Region of vulnerability')

    plt.text(0.05, 0.05, f'Cumulative Vulnerability: {auc:.4f}',
             transform=plt.gca().transAxes, fontsize=12, bbox=dict(facecolor='white', alpha=0.5))
    plt.title('SIRSV Model (Immunity waning)', fontsize=16)
    plt.xlabel("Months since start of simulation", fontsize=12)
    plt.ylabel("Fraction of population", fontsize=12)
    plt.ylim(0, 1)
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'{scenario}_waning_{model_type}.png'), dpi=300, bbox_inches='tight')
    plt.close()


def plot_parameter_sweep(results, param1_name, param2_name, output_variable='protected', vaccine_efficacy=1, herd_threshold=0.416, model_type='random'):
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
