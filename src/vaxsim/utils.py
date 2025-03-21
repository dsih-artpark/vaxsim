import logging
import random
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.integrate import simpson
from tqdm import tqdm
import yaml


def load_params():
    """Load model parameters from params.yaml."""
    params_path = Path(__file__).parent.parent.parent / 'params.yaml'
    with params_path.open('r') as f:
        return yaml.safe_load(f)


def auc_below_threshold(S, I, R, V, days, herd_threshold=0.416):
    """
    Helper function to compute the normalized area under the curve below a threshold value
    for the protected fraction, and normalize by total days.

    Args:
    S: Susceptible population array over time.
    I: Infected population array over time.
    R: Recovered population array over time.
    V: Vaccinated population array over time.
    days: Total number of days.
    herd_threshold: Threshold for the protected fraction.

    Returns:
    auc_normalised: Area under the curve normalized by the total number of days.
    """
    N = S + I + R + V
    protected = (R + V) / (N - I)
    t = np.arange(days) / 30

    area_below_threshold = np.maximum(herd_threshold - protected, 0) * (protected < herd_threshold)
    shaded_area = simpson(area_below_threshold, x=t)

    total_area = simpson(np.ones_like(protected), x=t)  # Area of the box from 0 to 1
    percentage_shaded = (shaded_area / total_area) * 100

    return percentage_shaded


def equilibrium_min_protected_fraction(S, I, R, V):
    """Helper function to compute minimum long-term protected fraction."""
    N = S + I + R + V
    protected_fraction = np.minimum((R + V) / (N - I), 1)
    return np.min(protected_fraction[-365:])


def compute_total_infections(I):
    """Helper function to compute total infections."""
    daily_infections = np.maximum(0, np.diff(I))
    return np.sum(daily_infections)


def model_loss(S, I, R, V, data, scale_diva=0.5):
    """Calculate loss between model predictions and observed data.

    This function computes the sum of squared errors between model predictions
    and seromonitoring/DIVA data.

    Args
    ----
    S : array-like
        Susceptible population over time
    I : array-like
        Infected population over time
    R : array-like
        Recovered population over time
    V : array-like
        Vaccinated population over time
    data : pd.DataFrame
        Data containing 'sero_eff' and 'diva' columns
    scale_diva : float, optional
        Scaling factor for DIVA predictions, by default 0.5

    Returns
    -------
    float
        Total loss value (or 1e6 if error occurs)

    Note
    ----
    Predictions are calculated as:
    - sero_pred = (R + V) / (N - I)
    - diva_pred = R / (N - I)
    where N = S + I + R + V
    """
    try:
        start_date = pd.to_datetime('2020-01-01')

        valid_sero = data['sero_eff'].dropna()
        valid_diva = data['diva'].dropna()
        sero_days = (valid_sero.index - start_date).days.astype(int)
        diva_days = (valid_diva.index - start_date).days.astype(int)

        max_day = max(max(sero_days), max(diva_days))
        if len(S) > max_day + 1:
            S = S[:max_day + 1]
            I = I[:max_day + 1]
            R = R[:max_day + 1] 
            V = V[:max_day + 1]

        N = S + I + R + V
        sero_pred = (R + V) / (N - I)
        diva_pred = R / (N - I)

        matched_sero = sero_pred[sero_days]
        matched_diva = diva_pred[diva_days]

        sero_error = ((valid_sero.values - matched_sero) ** 2).sum()
        diva_error = ((scale_diva * valid_diva.values - matched_diva) ** 2).sum()

        total_error = sero_error + diva_error
        print(f"Current loss: {total_error:.6f}")
        return total_error

    except (ValueError, ZeroDivisionError, IndexError) as e:
        print(f"Error in loss calculation: {e!s}")
        return 1e6


def analyse_scenarios(sirsv_model, params, output_dir='output', model_type='random'):
    """
    Analyses each scenario to compute total infections, percentage of infections averted, and protected fraction.
    Saves the results as both a CSV file and a PNG image.

    Parameters:
        params (dict): Dictionary containing parameters for each scenario.
        output_dir (str): Directory to save the output files (CSV and PNG).
    """
    results = []

    # Run simulations for each scenario
    filtered_scenarios = {k: v for k, v in params.items() if k == 'baseline' or k.startswith('scenario_')}

    # Run simulations for each scenario
    simulation_results = {}
    for scenario, scenario_params in filtered_scenarios.items():
        try:
            # Run model with current scenario parameters
            logging.info(f"Running scenario: {scenario} with params: {scenario_params}")
            S, I, R, V = sirsv_model(scenario_params, scenario)

            # Store simulation results
            simulation_results[scenario] = {'S': S, 'I': I, 'R': R, 'V': V}

        except Exception as e:
            logging.error(f"Error running scenario {scenario}: {e}")

    # Baseline scenario results for comparison
    baseline_data = simulation_results.get('baseline', {})
    baseline_infections = compute_total_infections(np.array(baseline_data.get('I', [])))

    # Analyze and compile results
    for scenario, data in simulation_results.items():
        I = np.array(data.get('I', []))
        R = np.array(data.get('R', []))
        V = np.array(data.get('V', []))
        S = np.array(data.get('S', []))

        # Calculate total infections
        total_infections = compute_total_infections(I)

        # Calculate percentage of infections averted
        I0 = params[scenario].get('I0', 'NA')
        if I0 == 'NA' or I0 == 0:
            infections_averted = 'NA'
        else:
            infections_averted = (
                (baseline_infections - total_infections) / baseline_infections * 100
                if baseline_infections > total_infections and scenario != 'baseline' else 'NA'
            )

        # Calculate protected fraction using the last 365 days
        last_365_days = slice(-365, None)
        total_population = S[last_365_days] + I[last_365_days] + R[last_365_days] + V[last_365_days]
        protected_fraction = np.min((R[last_365_days] + V[last_365_days]) / total_population) if total_population.any() else 'NA'
        cumulative_vulnerability = auc_below_threshold(S, I, R, V, len(I), herd_threshold=0.416)

        # Append results for this scenario
        vax_rate = params[scenario].get('vax_rate', 'NA')
        vax_period = params[scenario].get('vax_period', 'NA')
        remarks = params[scenario].get('Remarks', 'NA')
        results.append([
            scenario,
            round(vax_rate, 4) if vax_rate != 'NA' else vax_rate,
            round(vax_period, 4) if vax_period != 'NA' else vax_period,
            round(I0, 2) if I0 != 'NA' else I0,
            remarks,
            round(total_infections, 4),
            round(infections_averted, 4) if infections_averted != 'NA' else infections_averted,
            round(protected_fraction, 4) if protected_fraction != 'NA' else protected_fraction,
            round(cumulative_vulnerability, 4)
        ])

    # Create DataFrame for the results
    columns = [
        'scenario.ID',
        'Vax Rate (day^-1)',
        'Vax Period (days)',
        'I0',
        'Remarks',
        'Total Infections',
        'Percentage of Infections Averted',
        'Protected Fraction',
        'Cumulative Vulnerability'
    ]
    df_results = pd.DataFrame(results, columns=columns)
    csv_path = Path(output_dir) / 'scenario_analysis.csv'
    df_results.to_csv(csv_path, index=False)

    fig, ax = plt.subplots(figsize=(12, 6))
    ax.axis('tight')
    ax.axis('off')

    col_widths = [0.3, 0.2, 0.2, 0.2, 1.0, 0.3, 0.3, 0.3]
    row_height = 0.1

    table = ax.table(cellText=df_results.values, colLabels=df_results.columns, cellLoc='center', loc='center')

    for i, col_width in enumerate(col_widths):
        for j in range(len(df_results) + 1):
            cell = table[j, i]
            cell.set_width(col_width)
            cell.set_height(row_height)

    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1, 1.5)

    for key, cell in table.get_celld().items():
        cell.set_edgecolor('black')
        cell.set_linewidth(0.5)
        if key[1] == 4:
            cell.set_text_props(fontfamily='sans-serif', wrap=True)
            cell.set_text_props(visible=True)

    png_path = Path(output_dir) / f'scenario_analysis_{model_type}.png'
    plt.savefig(png_path, bbox_inches='tight', dpi=300)
    plt.close()

    latex_path = Path(output_dir) / f'scenario_analysis_{model_type}.tex'
    try:
        with open(latex_path, 'w') as f:
            latex_code = df_results.to_latex(index=False, column_format='|l|' + 'r|'*7, na_rep='NA')
            f.write(latex_code)
        logging.info(f"Results saved to LaTeX: {latex_path}")
    except Exception as e:
        logging.error(f"Error saving LaTeX file: {e}")

    logging.info(f"Results saved to CSV, PNG, and LaTeX files in {output_dir}.")


def run_parameter_sweep(sirsv_model, base_params, param1_name, param1_range, param2_name, param2_range,
                        diagonal=False, analysis_function='auc', model_type='random'):
    results = []

    # Function mapping for analysis type
    analysis_func = auc_below_threshold if analysis_function == 'auc' else equilibrium_min_protected_fraction

    if diagonal:
        # Ensure param1_range and param2_range are of equal length for diagonal sweep
        param_range = np.minimum(param1_range, param2_range)
        total_combinations = len(param_range)
    else:
        total_combinations = len(param1_range) * len(param2_range)

    with tqdm(total=total_combinations, desc="Parameter sweep progress") as pbar:
        if diagonal:
            # Diagonal sweep (where param1_value == param2_value)
            for value in param_range:
                current_params = base_params.copy()
                current_params[param1_name] = value
                current_params[param2_name] = value
                logging.info(f"Running params: {param1_name}, {param2_name}")

                try:
                    S, I, R, V = sirsv_model(current_params, "parameter_sweep")

                    # Check for NaN values in results
                    if np.isnan(S).any() or np.isnan(I).any() or np.isnan(R).any() or np.isnan(V).any():
                        logging.warning(f"NaN detected for param1={value} and param2={value}. Skipping this run.")
                        pbar.update(1)
                        continue

                    # Perform analysis using the selected function
                    analysis_result = analysis_func(S, I, R, V)

                    equilibrium_values = {
                        param1_name: value,
                        param2_name: value,
                        analysis_function: analysis_result
                    }
                    results.append(equilibrium_values)

                except Exception as e:
                    logging.error(f"Error occurred for param1={value} and param2={value}: {e}")

                pbar.update(1)

        else:
            # Regular sweep with independent param1 and param2 ranges
            for param1_value in param1_range:
                for param2_value in param2_range:
                    current_params = base_params.copy()
                    current_params[param1_name] = param1_value
                    current_params[param2_name] = param2_value

                    try:
                        S, I, R, V = sirsv_model(current_params, "parameter_sweep")

                        # Check for NaN values in results
                        if np.isnan(S).any() or np.isnan(I).any() or np.isnan(R).any() or np.isnan(V).any():
                            logging.warning(f"NaN detected for param1={param1_value} and param2={param2_value}. Skipping this run.")
                            pbar.update(1)
                            continue

                        # Perform analysis using the selected function
                        analysis_result = analysis_func(S, I, R, V)

                        equilibrium_values = {
                            param1_name: param1_value,
                            param2_name: param2_value,
                            analysis_function: analysis_result
                        }
                        results.append(equilibrium_values)

                    except Exception as e:
                        logging.error(f"Error occurred for param1={param1_value} and param2={param2_value}: {e}")

                    pbar.update(1)

    return results


def seed_infection(t, seed_schedule, seed_rate=0):
    """
    Returns the number of new seeds (infections) at time t based on a predefined schedule.

    Args:
        t: Current time step (day).
        seed_schedule: A list of days when seeding should occur.
        seed_rate: Number of seeds (infections) to introduce at the scheduled times.

    Returns:
        new_seeds: Number of new seeds (infections) at time t.
    """
    if t in seed_schedule:
        return seed_rate
    return 0


def generate_seed_schedule(method, min_day=None, max_day=None, days=None, num_seeds=None, event_series=None):
    """
    Generates a seed schedule based on the specified method.

    Args:
        method (str): The method for generating the schedule, either "random" or "event_series".
        min_day (int, optional): The minimum day for seeding in the "random" method.
        max_day (int, optional): The maximum day for seeding in the "random" method.
        days (int, optional): Total length of the binary schedule for "random" method.
        num_seeds (int, optional): Number of seeding events for the "random" method.
        event_series (list of int, optional): A manually defined binary series indicating seeding days 
                                              for the "event_series" method.

    Returns:
        list: A list of indices indicating the seeding days.

    Raises:
        ValueError: If required parameters are missing or if an invalid method is provided.
    """
    if method == "random":
        if min_day is None or max_day is None or days is None or num_seeds is None:
            raise ValueError("All parameters min_day, max_day, days, and num_seeds must be specified for 'random' method.")
        seed_schedule = [0] * days
        seed_days = random.sample(range(min_day, max_day), num_seeds)
        for day in seed_days:
            seed_schedule[day] = 1
        return seed_schedule
    elif method == "event_series":
        if event_series is None:
            raise ValueError("For 'event_series' method, event_series must be specified.")
        seed_schedule = [i for i in range(len(event_series)) if event_series[i] == 1]
        return seed_schedule
    else:
        raise ValueError("Invalid method. Choose 'random' or 'event_series'.")


def find_local_minima(data, days):
    """
    Finds local minima in a list of data points and returns a binary event series.

    Args:
        data: A list or array of data points.
        days: Total number of days (length of data).

    Returns:
        event_series: A binary list of length `days` with `1` at local minima indices and `0` elsewhere.
    """
    event_series = [0] * days
    for i in range(1, days - 1):
        if data[i - 1] > data[i] < data[i + 1]:
            event_series[i] = 1
    return event_series

