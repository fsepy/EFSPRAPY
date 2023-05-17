# Calculation of probability of exceeding limiting radiation level at a receiver
from efsprapy.func.functions import *


def main(
    N: int,  # Number of simulations
    fire_duration: float,  # fire duration [min]
    time_step: float,  # time step [sec]
    h_opening: float,  # Height of opening [m]
    area_opening: float,  # Area of opening [m²]
    area_floor: float,  # Floor area [m²]
    area_total: float,  # Total compartment area [m²]
    fire_load_density_mean: float,  # Mean fire load density [MJ/m²]
    comb_eff: float,  # Combustion efficiency
    lining_rho: float,  # Density of compartment lining [kg/m³]
    lining_cp: float,  # Specific heat of compartment [J/kgK]
    lining_k: float,  # Thermal conductivity of compartment [W/mK]
    limiting_time: float,   # Time to reach maximum temperature for fuel controlled fire [min]
    panel_width: float,  # Width of radiating panel [m]
    panel_height: float,  # Height of radiating panel [m]
    boundary_distance: float,  # Emitter to boundary [m]):
    heat_flux_critical: float,  # Critical heat flux of receiver material [kW/m²]
    n_factor: float,  # FTP index [n_factor = 1 (for thermally thin); n_factor = 2 (for thermally thick); n_factor = 1.5 (for intermediate)]
):
    time = list(np.arange(0, fire_duration * 60 + time_step, time_step, dtype=float))  # Time in [sec]
    fire_load_density_sd = 126  # Standard deviation of fire load density [MJ/m²]
    design_fire_load_density = gumbel_dist(fire_load_density_mean, fire_load_density_sd, N)

    # w_opening = area_opening / h_opening

    emissivity = np.random.uniform(0.7, 1.0, N)

    view_factor = view_factor_parallel(panel_width, panel_height, boundary_distance)

    heat_flux_s = np.zeros((N, len(time)), dtype=float)
    heat_flux_r = np.zeros((N, len(time)), dtype=float)
    ftp = np.zeros((N, len(time)), dtype=float)

    for i in range(0, N):
        fire_temp = en_parametric_fire(
            fire_duration,
            time_step,
            h_opening,
            area_opening,
            area_floor,
            area_total,
            design_fire_load_density[i],
            comb_eff,
            lining_rho,
            lining_cp,
            lining_k,
            limiting_time
        )
        for j in range(len(time)):
            heat_flux_s[i, j] = heat_flux_source(fire_temp[j], emissivity[i])  # Heat flux radiated from source [kW]
            heat_flux_r[i, j] = view_factor * heat_flux_s[i, j]  # Heat flux received at boundary due to radiation from source [kW]
        ftp[i, :] = flux_time_product(heat_flux_r[i, :], heat_flux_critical, time, time_step, n_factor)  # Flux time product

    return ftp

