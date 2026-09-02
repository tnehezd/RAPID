import os
import numpy as np
from scipy.special import ive
import yaml

def load_config():
    yaml_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "ring_test.yaml"))
    if not os.path.exists(yaml_path):
        yaml_path = "ring_test.yaml"
    with open(yaml_path, "r") as f:
        config = yaml.safe_load(f)
    return config


def analytical_viscous_ring(r, r_peak, ring_width, ring_mass, T):
    """
    Lynden-Bell & Pringle (1974) viscous ring analytical solution.
    Fully compatible with the C simulation, even with Gaussian initial condition.
    """

    # --- 1. Gaussian initial mass ---
    M0 = ring_mass

    # --- 2. LB&P normalization ---
    sigma0 = M0 / (2.0 * np.pi * r_peak)

    # --- 3. Dimensionless radius ---
    r_safe = np.where(r == 0.0, 1e-10, r)
    x = r_safe / r_peak

    # --- 4. T = 0 case: pure Gaussian ---
    if T == 0.0:
        return M0/(2.0 * np.pi * r_peak) * np.exp(-((r_safe - r_peak)**2) / (2.0 * ring_width**2))

    # --- 5. LB&P Green's function ---
    prefactor = sigma0 / T
    expo = np.exp(-(x - 1.0)**2 / T)
    bessel = ive(0.25, 2.0 * x / T)

    sigma = prefactor * x**(-0.25) * expo * bessel

    return sigma


if __name__ == "__main__":
    config = load_config()

    r0 = config["test"]["ring_center_au"]
    ring_width = config["test"]["ring_width_au"]
    grid_points = config["disk_parameters"]["number_of_grid_points"]
    r_min = config["disk_parameters"]["inner_radius_au"]
    r_max = config["disk_parameters"]["outer_radius_au"]
    ring_mass = config["test"]["ring_mass"]

    r_grid = np.linspace(r_min, r_max, grid_points)

    out_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "analytical_outputs"))
    os.makedirs(out_dir, exist_ok=True)

    dimensionless_times = [0.0, 0.004, 0.008, 0.016, 0.032, 0.064, 0.128, 0.256]

    for T in dimensionless_times:
        sigma_an = analytical_viscous_ring(r_grid, r0, ring_width, ring_mass, T)

        t_str = f"{T:.4f}".replace(".", "_")
        filename = f"analytical_T_{t_str}.dat"
        filepath = os.path.join(out_dir, filename)

        data_to_save = np.column_stack((r_grid, sigma_an))
        np.savetxt(filepath, data_to_save,
                header=f"Analytical Reference (C-matched) for T={T}",
                fmt="%.6e")

        print(f"Saved: {filepath}")
