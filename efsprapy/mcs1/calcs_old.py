import numpy as np
from fsetools.lib.fse_activation_hd import heat_detector_temperature_pd7974
from fsetools.lib.fse_din_en_1991_1_2_parametric_fire import temperature as din_param_temperature
from fsetools.lib.fse_thermal_radiation import phi_parallel_any_br187


def calculate_incident_heat_flux_from_controlled_fire_simple(
        fire_mode: int,
        t_arr: np.ndarray,

        A_t: float,
        A_f: float,
        H: float,
        b: float,

        W_v: float,
        H_v: float,

        q_x_d: float,
        hrr_density_kWm2: float,
        alpha: float,
        t_alpha: float,

        H_d: float,
        R_d: float,
        RTI: float,
        C: float,

        T_act: float,
        S: float,

        C_conv: float = 0.7,
        sigma: float = 5.67e-8,
        epsilon_s: float = 1.0,
        epsilon_f: float = 1.0,
):
    # A potential issue with fire models is the possibility of unrealistically high heat release rates (HRR) when the
    # fire is detected. To address this, we propose a simple radial spread fire model that replaces the continuous
    # t-square fire model. This new model caps the peak HRR based on the heat release rate per unit area (HRRPUA) and
    # the fuel density.
    #
    # The fire's area is computed in two parts:
    #
    #   Fire Area 1: Represents the total area that the fire has affected. It's calculated as:
    #       fire_area_1 = pi * (fire_travel_speed * time) ** 2
    #
    #   Fire Area 2: Represents the area where the fire has exhausted all fuel. It's computed as:
    #       fire_area_2 = pi * (fire_travel_speed * max(time - fuel_density / hrr_density, 0)) ** 2
    #
    # The overall fire area is then the difference between Fire Area 1 and Fire Area 2:
    # fire_area = fire_area_1 - fire_area_2
    #
    # The fire's travel speed is computed with the following formula derived from equating the total heat content in
    # Fire Area 1 with the energy growth of a t-square fire:
    #   fire_travel_speed = (alpha * 1e3 / hrr_density / pi) ** 0.5
    #
    # Here, alpha is the fire growth coefficient of the t-square fire model.
    fire_travel_speed = (alpha / hrr_density_kWm2 / np.pi) ** 0.5
    fire_area_1 = np.where((_ := np.pi * (fire_travel_speed * t_arr) ** 2) > A_f, A_f, _)
    fire_area_2 = np.pi * (fire_travel_speed * np.where((_ := t_arr - q_x_d / (hrr_density_kWm2 * 1e3)) < 0, 0, _)) ** 2
    fire_area = np.where((_ := (fire_area_1 - fire_area_2)) < 0, 0, _)
    fire_hrr_kW = fire_area * hrr_density_kWm2

    # ==============================================================
    # calculate fire HRR and jet temperature at sprinkler activation
    # ==============================================================
    *_, T_jet, T_d = heat_detector_temperature_pd7974(
        fire_time=t_arr,
        fire_hrr_kW=fire_hrr_kW,
        detector_to_fire_vertical_distance=H_d,
        detector_to_fire_horizontal_distance=R_d,
        detector_response_time_index=RTI,
        detector_conduction_factor=C,
        fire_hrr_density_kWm2=hrr_density_kWm2,
        fire_conv_frac=C_conv,
        activation_temperature=T_act,
    )
    if not (np.nanmin(T_d) < T_act < np.nanmax(T_d)):
        # no detection, no treatment to HRR and smoke temperature
        t_d = np.inf
        T_jet[np.isnan(T_jet)] = 0.
    elif fire_mode == 1:
        # fire capped at sprinkler activation
        t_d = t_arr[len(T_d)]
        T_jet = np.append(T_jet, np.full((len(t_arr) - len(T_jet),), T_jet[-1]))
        fire_hrr_kW[len(T_d):] = fire_hrr_kW[len(T_d)]
    elif fire_mode == 2:
        # fire decay at sprinkler activation
        t_d = t_arr[len(T_d)]
        res = dict(t_2_x=-1, t_3_x=-1, Q_2=-1, T_2_x=-1, T_3_x=-1, )
        din_param_temperature(
            t=t_arr, A_w=W_v * H_v, h_w=H_v, A_t=A_t, A_f=A_f, t_alpha=t_alpha, b=b, q_x_d=q_x_d, outputs=res
        )
        t_2_x, t_3_x, Q_2, T_2_x, T_3_x = res['t_2_x'], res['t_3_x'], res['Q_2'], res['T_2_x'], res['T_3_x']
        fire_hrr_kW[len(T_d):] = fire_hrr_kW[len(T_d)] + (np.arange(len(t_arr) - len(T_d), dtype=float)) * (
                0 - Q_2 * 1e3) / (t_3_x - t_2_x)
        fire_hrr_kW[fire_hrr_kW < 0] = 0

        T_jet = np.append(
            T_jet,
            T_jet[-1] + (np.arange(len(t_arr) - len(T_d), dtype=float)) * (T_3_x - T_2_x) / (t_3_x - t_2_x)
        )
        T_jet[T_jet < T_jet[0]] = T_jet[0]
    else:
        raise ValueError(f'Unknown `fire_mode` {fire_mode}.')

    # ================================
    # Calculate radiation from hot gas
    # ================================
    q_s = sigma * epsilon_s * T_jet ** 4  # Incident radiation at source due to smoke

    # view factor from hot gas
    phi_s = phi_parallel_any_br187(W_m=W_v, H_m=H_v, w_m=W_v / 2, h_m=H_v / 2, S_m=S, )

    # ==============================
    # Calculate radiation from flame
    # ==============================
    # flame type 1
    q_f_1 = (1 - C_conv) * fire_hrr_kW / (4 * np.pi * np.power(S, 2))

    # flame type 2
    D_f = (4 * (fire_hrr_kW / hrr_density_kWm2) / np.pi) ** 0.5
    H_f = 0.235 * np.power(fire_hrr_kW, 2. / 5.) - 1.02 * D_f
    W_e = np.where(H_f <= H, D_f, np.maximum(2 * 0.95 * (H_f - H), D_f))
    W_e[W_e > W_v] = W_v
    H_e = np.where(H_f > H, H, H_f)
    q_f_2_mask = np.logical_and(W_e > 0, H_e > 0)
    q_f_2 = np.zeros_like(fire_hrr_kW)
    q_f_2[q_f_2_mask] = ((1 - C_conv) * fire_hrr_kW[q_f_2_mask] / 2.) / (W_e[q_f_2_mask] * H_e[q_f_2_mask])

    # combined flame type 1 and flame type 2
    q_f = np.zeros_like(fire_hrr_kW)
    q_f_mask = D_f > 0
    q_f[q_f_mask] = epsilon_f * np.where(S / D_f[q_f_mask] > 2.5, q_f_1[q_f_mask], q_f_2[q_f_mask])
    q_f *= epsilon_f * 1e3  # Radiation at receiver due to flame [kw/m²] and convert to W/m²

    # view factor
    phi_f_mask = np.logical_and(q_f_mask, q_f_2_mask)

    phi_f_2 = list()
    for i in np.arange(len(phi_f_mask))[phi_f_mask]:
        phi_f_2.append(phi_parallel_any_br187(W_m=W_e[i], H_m=H_e[i], w_m=W_e[i] / 2, h_m=H_e[i] / 2, S_m=S, ))
    phi_f_2 = np.array(phi_f_2)

    phi_f = np.zeros_like(fire_hrr_kW)
    phi_f[phi_f_mask] = np.where(
        (S / D_f[phi_f_mask]) > 2.5,
        np.full_like(q_f_1[phi_f_mask], phi_parallel_any_br187(W_m=W_v, H_m=H_v, w_m=W_v / 2, h_m=H_v / 2, S_m=S, )),
        phi_f_2
    )

    return q_f, phi_f, q_s, phi_s, t_d  # fire hrr is in kW but need to be SI unit when returned
