import os
import sys
import traceback
from tempfile import TemporaryDirectory
from typing import Optional

import numpy as np
from fsetools.lib.fse_bs_en_1991_1_2_parametric_fire import temperature as para_fire
from fsetools.lib.fse_thermal_radiation import phi_parallel_any_br187

from ..cfast.cfast import Run as RunCFAST
from ..cfast.test_files import simple
from ..project_info import logger
from ..safir.therm2d import Run, PPXML


def calculate_hrr_and_smoke_with_sprinkler_suppression_cfast(
        t_arr: np.ndarray,

        W: float,
        D: float,
        H: float,

        q_fd: float,
        hrr_density_kWm2: float,
        alpha_kWs2: float,

        H_d: float,
        R_d: float,
        RTI: float,
        T_act: float,

        dir_temp: Optional[str] = None,
):
    """

    :param t_arr: [s], an array represents the time
    :param W: [m], room floor width
    :param D: [m], room floor depth
    :param H: [m], room height
    :param q_fd: [W/m2], fuel load density
    :param hrr_density_kWm2: [kW/m2], heat release rate density
    :param alpha_kWs2: [kW/s2], fire growth rate
    :param H_d: [m], sprinkler / heat detector height (vertical distance above the fire)
    :param R_d: [m], sprinkler / heat detector radial distance from the fire
    :param RTI: [m0.5 s0.5], sprinkler / heat detector response time index
    :param T_act: [K], sprinkler / heat detector activation temperature
    :return:
    """
    A_t = W * D
    A_f = 2 * (A_t + D * H + H * W)
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
    t_arr_ = np.arange(0, t_arr[-1] + 1, 30., dtype=float)
    fire_travel_speed = (alpha_kWs2 / hrr_density_kWm2 / np.pi) ** 0.5
    fire_area_1 = np.where((_ := np.pi * (fire_travel_speed * t_arr_) ** 2) > A_f, A_f, _)
    fire_area_2 = np.pi * (
            fire_travel_speed * np.where((_ := t_arr_ - q_fd / (hrr_density_kWm2 * 1e3)) < 0, 0, _)) ** 2
    fire_area = np.where((_ := (fire_area_1 - fire_area_2)) < 0, 0, _)
    fire_hrr_kW = fire_area * hrr_density_kWm2

    # the calculated fire hrr to be converted into following format:
    fire_hrr_curve_tabl = (
        "&TABL ID = 'Constant Fire' "
        "LABELS = 'TIME','HRR','HEIGHT','AREA','CO_YIELD','SOOT_YIELD','HCN_YIELD','HCL_YIELD','TRACE_YIELD' /\n"
    )
    # &TABL ID = 'Constant Fire', DATA = 0,    0,   0, 0.01, 0, 0, 0, 0, 0 /
    # &TABL ID = 'Constant Fire', DATA = 10,   100, 0, 0.01, 0, 0, 0, 0, 0 /
    # &TABL ID = 'Constant Fire', DATA = 990,  100, 0, 0.01, 0, 0, 0, 0, 0 /
    # &TABL ID = 'Constant Fire', DATA = 1000, 0,   0, 0.01, 0, 0, 0, 0, 0 /
    fire_hrr_curve_tabl += "&TABL ID = 'Constant Fire', DATA = 0, 0, 0, 0, 0, 0.07, 0, 0, 0 /"
    fire_hrr_curve_tabl += '\n'.join(filter(
        None,
        [
            f"&TABL "
            f"ID = 'Constant Fire', "
            f"DATA = {t_arr_[i]:.0f}, {fire_hrr_kW[i]:.1f}, 0, {fire_area[i]:.2f}, 0, 0.07, 0, 0, 0 /"
            if fire_hrr_kW[i - 1] != fire_hrr_kW[i] != fire_hrr_kW[i + 1] else ''
            for i in range(1, len(t_arr_) - 1)
        ]
    ))

    t_end = np.amax(t_arr_[fire_hrr_kW > np.amin(fire_hrr_kW)])

    # calculate sprinkler location
    if ((W ** 2 + D ** 2) ** 0.5 / 4.) < R_d:
        sprinkler_loc_x = W / 2
        sprinkler_loc_y = D / 2
    else:
        # x / y = W / D
        # x = W / D * y
        # (x**2 + y**2) ** 0.5 = R_d
        # ((W / D * x)**2 + x**2) ** 0.5 = R_d
        sprinkler_loc_x = (R_d * D) / ((D ** 2 + W ** 2) ** 0.5)
        sprinkler_loc_y = sprinkler_loc_x / W * D

    # ===================================================================================
    # calculate smoke layer temperature and actual fire hrr considering sprinkler effects
    # ===================================================================================
    with TemporaryDirectory(dir=dir_temp) as dir_work:
        fn = f'a'
        fn_cfast_in = f'{fn}.in'
        fp_cfast_in = os.path.join(dir_work, fn_cfast_in)
        with open(fp_cfast_in, 'w+') as f_in:
            f_in.write(simple.format(
                t_end=t_end,
                t_step=t_arr[1] - t_arr[0],
                room_width=W,
                room_depth=D,
                room_height=H,
                opening_width=max(2., W / 2.),
                opening_height=0.5,
                fire_hrr_curve_tabl=fire_hrr_curve_tabl,
                sprinkler_loc_x=sprinkler_loc_x,
                sprinkler_loc_y=sprinkler_loc_y,
                sprinkler_loc_z=H_d - 0.02,
                sprinkler_activation_temperature=T_act - 273.15,
                sprinkler_rti=RTI,
            ))
        _ = RunCFAST().run(fp_cfast_in)

        t_, ult_, llt_, lh_, hrr_, spt_ = _.read_outputs()

    ult_ = np.interp(t_arr, t_, ult_) + 273.15
    # llt_ = np.interp(t_arr, t_, llt_)
    # lh_ = np.interp(t_arr, t_, lh_)
    hrr_ = np.interp(t_arr, t_, hrr_)
    t_d = np.interp(T_act - 273.15, spt_, t_)

    return hrr_, ult_, t_d


def calculate_incident_heat_flux_from_sprinkler_suppressed_fire(
        t_arr: np.ndarray,
        fire_hrr_kW: np.ndarray,
        T_smoke: np.ndarray,
        hrr_density_kWm2: float,

        H: float,
        W_o: float,
        H_o: float,
        S: float,

        C_conv: float = 0.7,
        sigma: float = 5.67e-8,
        epsilon_s: float = 1.0,
        epsilon_f: float = 1.0,

        t_cap: float = 0.,
):
    # if `t_d` is provided, then cap fire heat release rate and smoke temperature at `t_d`. namely, steady-state after
    # activation time `t_d`.
    if t_cap > 0.:
        fire_hrr_kW = np.where(t_arr > t_cap, np.interp(t_cap, t_arr, fire_hrr_kW), fire_hrr_kW)
        T_smoke = np.where(t_arr > t_cap, np.interp(t_cap, t_arr, T_smoke), T_smoke)

    # ================================
    # Calculate radiation from hot gas
    # ================================
    q_s = sigma * epsilon_s * T_smoke ** 4  # Incident radiation at source due to smoke

    # view factor from hot gas
    phi_s = phi_parallel_any_br187(W_m=W_o, H_m=H_o, w_m=W_o / 2, h_m=H_o / 2, S_m=S, )

    # ==============================
    # Calculate radiation from flame
    # ==============================
    # flame type 1
    q_f_1 = (1 - C_conv) * fire_hrr_kW / (4 * np.pi * np.power(S, 2))

    # flame type 2
    D_f = (4 * (fire_hrr_kW / hrr_density_kWm2) / np.pi) ** 0.5
    H_f = 0.235 * np.power(fire_hrr_kW, 2. / 5.) - 1.02 * D_f
    W_e = np.where(H_f <= H, D_f, np.maximum(2 * 0.95 * (H_f - H), D_f))
    W_e[W_e > W_o] = W_o
    H_e = np.where(H_f > H_o, H_o, H_f)
    q_f_2_mask = np.logical_and(W_e > 0, H_e > 0)
    q_f_2 = np.zeros_like(fire_hrr_kW)
    q_f_2[q_f_2_mask] = ((1 - C_conv) * fire_hrr_kW[q_f_2_mask] / 4.) / (W_e[q_f_2_mask] * H_e[q_f_2_mask])

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
        np.full_like(q_f_1[phi_f_mask], phi_parallel_any_br187(W_m=W_o, H_m=H_o, w_m=W_o / 2, h_m=H_o / 2, S_m=S, )),
        phi_f_2
    )
    return q_f, phi_f, q_s, phi_s


def calculate_incident_heat_flux_from_parametric_fire(
        t_arr: np.ndarray,
        W: float,
        D: float,
        H: float,
        q_fd: float,
        lbd: float,
        rho: float,
        c: float,
        t_lim: float,

        W_o: float,
        H_o: float,
        S: float,
):
    # calculate the view factors
    T = para_fire(
        t=t_arr,
        A_t=2 * (W * D + D * H + H * W),
        A_f=W * D,
        A_v=W_o * H_o,
        h_eq=H_o,
        q_fd=q_fd,
        lbd=lbd,
        rho=rho,
        c=c,
        t_lim=t_lim,
    )
    q = 5.67e-8 * 1.0 * T ** 4
    phi = phi_parallel_any_br187(W_m=W_o, H_m=H_o, w_m=W_o / 2, h_m=H_o / 2, S_m=S, )

    return q, phi


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


def calculate_ignition_time_from_temperature(t_arr: np.ndarray, q_inc: np.ndarray, T_ig: float, safir_in_s: str,
                                             t_step: float = 5., dir_temp: str = None):
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
    if T_ig < 0:
        return np.nan, np.nan, np.nan
    if (np.amax(q_inc) / (5.67e-8 * 1.0) + 293.15 ** 4) ** 0.25 < T_ig:
        return np.nan, np.nan, np.nan

    t_end = np.amax(t_arr[q_inc > np.amin(q_inc)])
    fn_bc = f'b'
    fn = f'i'
    fn_safir_in = f'{fn}.in'
    safir_in_formatted_s = safir_in_s.format(fn_bc=fn_bc, t_step=t_step, t_end=t_end, )

    with TemporaryDirectory(dir=dir_temp) as dir_work:
        with open(os.path.join(dir_work, fn_bc), 'w+') as f_bc:
            f_bc.write('\n'.join(f'{t_arr[i]:g}, {q_inc[i]:.3f}' for i in range(len(t_arr))))

        fp_safir_in = os.path.join(dir_work, fn_safir_in)
        with open(fp_safir_in, 'w+') as f_in:
            f_in.write(safir_in_formatted_s)
        Run().run(fp_safir_in)

        try:
            with open(os.path.join(dir_work, f'{fn}.XML')) as f_xml:
                pp = PPXML(xml=f_xml.read())
            T_1 = np.interp(t_arr, pp.t, pp.get_nodes_temp(np.array([1]))[0, :] + 273.15)
        except ValueError:
            T_1 = None

    if T_1 is not None:
        T_max = np.amax(T_1)
        t_max = t_arr[np.argmax(T_1)]
        try:
            t_ig = np.amin(t_arr[T_1 >= T_ig])
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
        opening_ventilation_fraction: float,

        room_height: float,
        room_width: float,
        room_width_depth_ratio: float,

        fire_mode: int,
        fire_fuel_density_MJm2: float,
        fire_combustion_efficiency: float,
        fire_hrr_density_kWm2: float,
        fire_growth_factor: float,
        fire_t_lim: float,

        detector_to_fire_vertical_distance: float,
        detector_to_fire_horizontal_distance: float,
        detector_act_temp: float,
        detector_response_time_index: float,

        receiver_separation: float,

        lining_rho: float,
        lining_c: float,
        lining_k: float,

        ftp_chf: Optional[float] = None,
        ftp_index: Optional[float] = None,
        ftp_target: Optional[float] = None,

        receiver_ignition_temperature: Optional[float] = -1,
        safir_input_file_s: Optional[str] = None,

        dir_temp: Optional[str] = None,
):
    """Calculates flux-time product based on PD 7974-1 Clause 8.2.2

    :param t_end: [s], end time
    :param t_step: [s], time step
    :param opening_width: [m]
    :param opening_height: [m], ventilation opening height
    :param room_height: [m]
    :param fire_mode: [1]
    :param fire_fuel_density_MJm2: [MJ/m^2]
    :param fire_growth_factor: [kW/s2]
    :param fire_hrr_density_kWm2:  [kW/m^2]
    :param fire_t_lim: [s]
    :param receiver_ignition_temperature:
    :param detector_response_time_index:
    :param detector_to_fire_horizontal_distance:
    :param detector_act_temp:
    :param detector_to_fire_vertical_distance:
    :param lining_rho: [kg/m^3] lining density
    :param lining_c: [J/K/kg] lining specific heat capacity
    :param lining_k: [W/m/K] lining thermal conductivity
    :param receiver_separation: [m], distance between emitter and receiver, used to calculate the view factor
    :param ftp_chf: [W], critical heat flux of the receiver surface
    :param ftp_index: [1], FTP index, 1 for thermally thin; 2 for thermally thick; 1.5 for intermediate
    :param ftp_target: [1], FTP target for ignition
    :param safir_input_file_s: [1], safir input file in string format
    :return:
    """
    fire_fuel_density = fire_fuel_density_MJm2 * 1e6 * fire_combustion_efficiency
    t_arr = np.arange(0, t_end + t_step / 2., t_step, dtype=float)
    room_depth = room_width / room_width_depth_ratio

    # ============================
    # Calculate incident heat flux
    # ============================
    if fire_mode == 0:
        q_f, phi_f = calculate_incident_heat_flux_from_parametric_fire(
            t_arr=t_arr,
            W=room_width, D=room_depth, H=room_height,
            q_fd=fire_fuel_density,
            lbd=lining_k, rho=lining_rho, c=lining_c,
            t_lim=fire_t_lim,
            W_o=opening_width, H_o=opening_height * opening_ventilation_fraction,
            S=receiver_separation,
        )
        q_inc = q_f * phi_f
        t_d = np.inf
    elif fire_mode == 1 or fire_mode == 2:
        hrr, T_smoke, t_d = calculate_hrr_and_smoke_with_sprinkler_suppression_cfast(
            t_arr=t_arr,
            W=room_width, D=room_depth, H=room_height,
            q_fd=fire_fuel_density, hrr_density_kWm2=fire_hrr_density_kWm2, alpha_kWs2=fire_growth_factor,
            H_d=detector_to_fire_vertical_distance, R_d=detector_to_fire_horizontal_distance,
            RTI=detector_response_time_index, T_act=detector_act_temp,
            dir_temp=dir_temp
        )

        q_f, phi_f, q_s, phi_s = calculate_incident_heat_flux_from_sprinkler_suppressed_fire(
            t_arr=t_arr, fire_hrr_kW=hrr * 1e-3,
            T_smoke=T_smoke, hrr_density_kWm2=fire_hrr_density_kWm2,
            H=room_height,
            W_o=opening_width, H_o=opening_height,
            S=receiver_separation,
            t_cap=0. if fire_mode == 1 else t_d
        )
        q_inc = q_f * phi_f
    else:
        logger.debug('{}\n{}\n{}'.format(*sys.exc_info()[:2], traceback.format_exc()))
        return [np.nan] * 9

    # except Exception as e:
    #     logger.debug('{}\n{}\n{}'.format(*sys.exc_info()[:2], traceback.format_exc()))
    #     return [np.nan] * 9

    # ==============================
    # Calculate ignition temperature
    # ==============================
    # flux-time product
    try:
        if ftp_target <= 0 or ftp_chf <= 0 or ftp_index <= 0:
            raise ValueError
        t_ig_ftp, ftp = calculate_ignition_time_ftp(
            t=t_arr, q_inc=q_inc, ftp_chf=ftp_chf, ftp_index=ftp_index, ftp_target=ftp_target,
        )
    except Exception:
        t_ig_ftp, ftp = np.nan, np.nan

    # surface temperature
    try:
        if receiver_ignition_temperature <= 0:
            raise ValueError
        t_ig_safir, t_max_safir, T_max_safir = calculate_ignition_time_from_temperature(
            t_arr=t_arr, q_inc=q_inc, T_ig=receiver_ignition_temperature, safir_in_s=safir_input_file_s,
            dir_temp=dir_temp
        )
    except Exception:
        t_ig_safir, t_max_safir, T_max_safir = np.nan, np.nan, np.nan

    # reduce vector to scalar for outputs
    return tuple(np.average(i) if isinstance(i, np.ndarray) else i for i in (
        q_inc, t_ig_ftp, ftp[-1], t_ig_safir, t_max_safir, T_max_safir, fire_mode, t_d, fire_fuel_density,
        fire_hrr_density_kWm2
    ))
