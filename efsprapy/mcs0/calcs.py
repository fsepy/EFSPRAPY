__all__ = (
    'main',
)

import os
import tempfile
from typing import Tuple

import numpy as np
from fsetools.lib.fse_activation_hd import heat_detector_temperature_pd7974
from fsetools.lib.fse_bs_en_1991_1_2_parametric_fire import temperature as para_fire
from fsetools.lib.fse_din_en_1991_1_2_parametric_fire import temperature as din_param_temperature
from fsetools.lib.fse_thermal_radiation import phi_parallel_any_br187

from ..safir.test_files import therm1d_hf_ft_20
from ..safir.therm2d import Run, PPXML

with open(therm1d_hf_ft_20, 'r') as f:
    therm1d_hf_ft_20_s = f.read()


def calculate_incident_heat_flux_from_controlled_fire(
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


def calculate_incident_heat_flux(
        fire_mode: int,

        t: np.ndarray,
        A_t: float,
        A_f: float,
        A_v: float,
        h_eq: float,
        q_fd: float,
        lbd: float,
        rho: float,
        c: float,
        t_lim: float,

        W_e: float,
        H_e: float,
        S: float,

        hrr_density_kWm2: float,
        alpha: float,
        H_d: float,
        R_d: float,
        RTI: float,
        C: float,
        T_act: float,
        H: float,
        b: float,
        C_conv: float,
        t_alpha: float,
):
    # calculate the view factors
    if fire_mode == 0:
        T = para_fire(
            t=t,
            A_t=A_t,
            A_f=A_f,
            A_v=A_v,
            h_eq=h_eq,
            q_fd=q_fd,
            lbd=lbd,
            rho=rho,
            c=c,
            t_lim=t_lim,
        )
        q_1 = 5.67e-8 * 1.0 * T ** 4
        q_2 = 0
        phi_1 = phi_parallel_any_br187(W_m=W_e, H_m=H_e, w_m=W_e / 2, h_m=H_e / 2, S_m=S, )
        phi_2 = 0
        t_d = np.nan
    elif fire_mode == 1 or fire_mode == 2:
        q_1, phi_1, q_2, phi_2, t_d = calculate_incident_heat_flux_from_controlled_fire(
            fire_mode=fire_mode,
            t_arr=t,
            hrr_density_kWm2=hrr_density_kWm2,
            alpha=alpha,
            H_d=H_d,
            R_d=R_d,
            RTI=RTI,
            C=C,
            T_act=T_act,
            H=H,
            S=S,
            A_f=A_f,
            A_t=A_t,
            H_v=h_eq,
            W_v=A_v / h_eq,
            b=b,
            q_x_d=q_fd,
            C_conv=C_conv,
            t_alpha=t_alpha,
        )
    else:
        raise ValueError('Unknown `fire_mode`')

    return q_1, phi_1, q_2, phi_2, t_d


def calculate_ignition_time_ftp(
        t: np.ndarray,
        q_inc: np.ndarray,
        ftp_chf: float,
        ftp_index: float,
        ftp_target: float,
):
    ftp = np.zeros_like(t)
    ftp_i_diff = (q_inc[:-1] + q_inc[1:]) * 0.5
    ftp_i_diff[ftp_i_diff < ftp_chf] = ftp_chf
    ftp_i = ((ftp_i_diff * 1e-3 - ftp_chf * 1e-3) ** ftp_index) * (t[1:] - t[:-1])
    ftp[1:] = np.cumsum(ftp_i)
    try:
        if ftp_target <= np.amax(ftp):
            t_ig = t[np.argmin(np.abs(ftp - ftp_target))]
        else:
            t_ig = np.inf
    except ValueError:
        t_ig = np.nan
    return t_ig, ftp


def calculate_ignition_time_temperature(t: np.ndarray, q_inc: np.ndarray, T_ig: float, safir_time_step: float = 5., ):
    # First to check , under ideal conditions, if the receiver temperature could reach ignition temperature,  which is
    # the highest possible temperature the receiver could attain through radiation following the Stefan-Boltzmann law:
    #
    #     heat_flux = sigma * epsilon * (emitter_temperature ** 4 - 293.15 ** 4)
    #
    # The above equation is rearranged to solve for the receiver temperature:
    #
    #     emitter_temperature = (heat_flux / (sigma * epsilon) + 293.15 ** 4) ** (1/4)
    #
    # The receiver temperature will not exceed emitter temperature. Thus, if the calculated emitter temperature is not
    # greater than the ignition temperature, the function should return immediately, avoiding further heat transfer.
    if (np.amax(q_inc) / (5.67e-8 * 1.0) + 293.15 ** 4) ** 0.25 < T_ig:
        return np.nan, np.nan, np.nan

    t_end = np.amax(t[q_inc > np.amin(q_inc)])

    with tempfile.TemporaryDirectory() as dir_work:
        fn_bc = f'b'
        with open(os.path.join(dir_work, fn_bc), 'w+') as f_bc:
            f_bc.write('\n'.join(f'{t[i]:g}, {q_inc[i]:.3f}' for i in range(len(t))))

        fn = f'a'
        fn_safir_in = f'{fn}.in'
        fp_safir_in = os.path.join(dir_work, fn_safir_in)
        with open(fp_safir_in, 'w+') as f_in:
            f_in.write(therm1d_hf_ft_20_s.format(
                fn_bc=fn_bc,
                materials=f'WOODEC5\n    450. 0 25 9 0.8 1.2 0.0 0.0 1.0',
                t_step=safir_time_step,
                t_end=t_end,
            ))
        Run().run(fp_safir_in)

        try:
            with open(os.path.join(dir_work, f'{fn}.XML')) as f_xml:
                pp = PPXML(xml=f_xml.read())
            T_1 = np.interp(t, pp.t, pp.get_nodes_temp(np.array([401]))[0, :] + 273.15)
        except ValueError:
            T_1 = None

    if T_1 is not None:
        T_max = np.amax(T_1)
        t_max = t[np.argmax(T_1)]
        try:
            t_ig = np.amin(t[T_1 >= T_ig])
        except ValueError:
            t_ig = np.inf
    else:
        return np.nan, np.nan, np.nan

    return t_ig, t_max, T_max


def main(
        t_end: float,
        t_step: float,

        opening_width: float,
        opening_height: float,

        room_height: float,
        room_floor_area: float,
        room_total_surface_area: float,

        fire_mode: int,
        fire_fuel_density_MJm2: float,
        fire_hrr_density_kWm2: float,
        fire_growth_factor: float,
        fire_t_lim: float,
        fire_din_growth_factor: float,
        fire_convection_factor: float,

        detector_to_fire_vertical_distance: float,
        detector_to_fire_horizontal_distance: float,
        detector_act_temp: float,
        detector_response_time_index: float,
        detector_conduction_factor: float,

        receiver_ignition_temperature: float,
        receiver_separation: float,

        lining_rho: float,
        lining_c: float,
        lining_k: float,

        ftp_chf: float,
        ftp_index: float,
        ftp_target: float,
) -> Tuple[float, float, float, float, float, float, float, float, float, float]:
    """Calculates flux-time product based on PD 7974-1 Clause 8.2.2

    :param t_end: [s], end time
    :param t_step: [s], time step
    :param opening_width: [m]
    :param opening_height: [m], ventilation opening height
    :param room_height: [m]
    :param room_floor_area: [m^2] room floor area
    :param room_total_surface_area: [m^2] room total internal surface area, including ventilation openings
    :param fire_mode: [1]
    :param fire_fuel_density_MJm2: [MJ/m^2]
    :param fire_convection_factor: [1]
    :param fire_din_growth_factor: [1]
    :param fire_growth_factor: [kW/s2]
    :param fire_hrr_density_kWm2:  [kW/m^2]
    :param fire_t_lim: [s]
    :param receiver_ignition_temperature:
    :param detector_conduction_factor:
    :param detector_response_time_index:
    :param detector_to_fire_horizontal_distance:
    :param detector_act_temp:
    :param detector_to_fire_vertical_distance:
    :param lining_rho: [kg/m^3] lining density
    :param lining_c: [J/K/kg] lining specific heat capacity
    :param lining_k: [W/m/K] lining thermal conductivity
    :param receiver_emissivity: [1], emissivity of the emitter, e.g., from the opening
    :param receiver_separation: [m], distance between emitter and receiver, used to calculate the view factor
    :param ftp_chf: [W], critical heat flux of the receiver surface
    :param ftp_index: [1], FTP index, 1 for thermally thin; 2 for thermally thick; 1.5 for intermediate
    :param ftp_target: [1], FTP target for ignition
    :return:
    """
    fire_fuel_density = fire_fuel_density_MJm2 * 1e6
    t_arr = np.arange(0, t_end + t_step / 2., t_step, dtype=float)

    q_1, phi_1, q_2, phi_2, t_d = calculate_incident_heat_flux(
        fire_mode=fire_mode,

        t=t_arr,
        A_t=room_total_surface_area,
        A_f=room_floor_area,
        A_v=opening_width * opening_height,
        h_eq=opening_height,
        q_fd=fire_fuel_density,
        lbd=lining_k,
        rho=lining_rho,
        c=lining_c,
        t_lim=fire_t_lim,

        W_e=opening_width,
        H_e=opening_height,
        S=receiver_separation,

        hrr_density_kWm2=fire_hrr_density_kWm2,
        alpha=fire_growth_factor,
        H_d=detector_to_fire_vertical_distance,
        R_d=detector_to_fire_horizontal_distance,
        RTI=detector_response_time_index,
        C=detector_conduction_factor,
        T_act=detector_act_temp,
        H=room_height,
        b=(lining_k * lining_rho * lining_c) ** 0.5,
        C_conv=fire_convection_factor,
        t_alpha=fire_din_growth_factor,
    )
    q_inc = q_1 * phi_1 + q_2 * phi_2

    # ==============================
    # Calculate ignition temperature
    # ==============================

    # flux-time product
    t_ig_ftp, ftp = calculate_ignition_time_ftp(
        t=t_arr,
        q_inc=q_1 * phi_1 + q_2 * phi_2,
        ftp_chf=ftp_chf,
        ftp_index=ftp_index,
        ftp_target=ftp_target,
    )

    # surface temperature
    try:
        q_2 = 0 if q_2 is np.nan else q_2
        phi_2 = 0 if phi_2 is np.nan else phi_2
        t_ig_safir, t_max_safir, T_max_safir = calculate_ignition_time_temperature(
            t=t_arr, q_inc=q_1 * phi_1 + q_2 * phi_2, T_ig=receiver_ignition_temperature,
        )
    except ValueError:
        t_ig_safir, t_max_safir, T_max_safir = np.nan, np.nan, np.nan

    if isinstance(phi_1, np.ndarray):
        phi_1 = np.average(phi_1)

    return phi_1, phi_2, ftp[-1], t_ig_ftp, t_ig_safir, t_max_safir, T_max_safir, fire_mode, t_d, np.amax(q_inc)
