__all__ = 'heat_flux_from_controlled_fire',

import numpy as np

from fsetools.lib.fse_activation_hd import heat_detector_temperature_pd7974


def heat_flux_from_controlled_fire(
        t: np.ndarray,
        fire_hrr_density_kWm2: float,
        fire_alpha,
        detector_to_fire_vertical_distance,
        detector_to_fire_horizontal_distance,
        detector_response_time_index,
        detector_conduction_factor,
        detector_act_temp,
        H,
        W_o,
        S,
        fire_conv_frac=0.7,
        sigma=5.67e-8,
        epsilon_s=1.0,
        epsilon_f=1.0,
):
    fire_hrr_kW = fire_alpha * (t) ** 2

    # ==============================================================
    # calculate fire HRR and jet temperature at sprinkler activation
    # ==============================================================
    *_, jet_temperature, detector_temperature = heat_detector_temperature_pd7974(
        fire_time=t,
        fire_hrr_kW=fire_hrr_kW,
        detector_to_fire_vertical_distance=detector_to_fire_vertical_distance,
        detector_to_fire_horizontal_distance=detector_to_fire_horizontal_distance,
        detector_response_time_index=detector_response_time_index,
        detector_conduction_factor=detector_conduction_factor,
        fire_hrr_density_kWm2=fire_hrr_density_kWm2,
        fire_conv_frac=fire_conv_frac,
        activation_temperature=detector_act_temp,
    )
    if not (min(detector_temperature) < detector_act_temp < max(detector_temperature)):
        raise ValueError(
            f'Heat detector not activated, {min(detector_temperature)} < {detector_act_temp} < {max(detector_temperature)} not satisfied'
        )
    t_act = t[np.argmin(np.abs(detector_temperature - detector_act_temp))]
    jet_temperature_act = jet_temperature[detector_temperature > detector_act_temp][0]
    fire_hrr_act_kW = fire_hrr_kW[detector_temperature > detector_act_temp][0]

    jet_temperature[t > t_act] = jet_temperature_act
    fire_hrr_kW[t > t_act] = fire_hrr_act_kW

    # ================================
    # Calculate radiation from hot gas
    # ================================
    q_s = sigma * epsilon_s * jet_temperature ** 4  # Incident radiation at source due to smoke

    # ==============================
    # Calculate radiation from flame
    # ==============================
    # flame type 1
    q_f_1 = (1 - fire_conv_frac) * fire_hrr_kW / (4 * np.pi * np.power(S, 2))

    # flame type 2
    D_f = (4 * (fire_hrr_kW / fire_hrr_density_kWm2) / np.pi) ** 0.5
    H_f = 0.235 * np.power(fire_hrr_kW, 2. / 5.) - 1.02 * D_f
    W_e = np.where(H_f <= H, D_f, np.maximum(2 * 0.95 * (H_f - H), D_f))
    W_e[W_e > W_o] = W_o
    H_e = np.where(H_f > H, H, H_f)
    q_f_2 = ((1 - fire_conv_frac) * fire_hrr_kW / 2.) / (W_e * H_e)

    # combined flame type 1 and flame type 2
    q_f = np.where(S / D_f > 2.5, q_f_1, q_f_2)
    q_f *= epsilon_f  # Radiation at receiver due to flame [kw/m²]

    return q_f, q_s
