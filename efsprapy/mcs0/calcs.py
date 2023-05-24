__all__ = (
    'main',
)

import numpy as np
from fsetools.lib.fse_bs_en_1991_1_2_parametric_fire import temperature as param_fire
from fsetools.lib.fse_thermal_radiation import phi_parallel_any_br187


def main(
        t_end: float,
        t_step: float,
        vent_width: float,
        vent_height: float,
        room_floor_area: float,
        room_total_surface_area: float,
        fuel_density: float,
        rho: float,
        c: float,
        k: float,
        t_lim: float,
        emissivity: float,
        emitter_width: float,
        emitter_height: float,
        emitter_receiver_separation: float,
        chf: float,
        ftp_index: float,
        ftp_target: float,
) -> tuple:
    """Calculates flux-time product based on PD 7974-1 Clause 8.2.2

    :param vent_width:
    :param t_end: [s], end time
    :param t_step: [s], time step
    :param vent_height: [m], ventilation opening height
    :param vent_area: [m], ventilation opening area
    :param room_floor_area: [m^2] room floor area
    :param room_total_surface_area: [m^2] room total internal surface area, including ventilation openings
    :param fuel_density: [MJ/m^2] fue load density
    :param rho: [kg/m^3] lining density
    :param c: [??] lining specific heat capacity
    :param k: [??] lining thermal conductivity
    :param t_lim:
    :param emissivity: [1], emissivity of the emitter, e.g., from the opening
    :param emitter_width: [m], emitter width, used to calculate the view factor
    :param emitter_height:  [m], emitter height, used to calculate the view factor
    :param emitter_receiver_separation: [m], distance between emitter and receiver, used to calculate the view factor
    :param chf: [W], critical heat flux of the receiver surface
    :param ftp_index: [1], FTP index, 1 for thermally thin; 2 for thermally thick; 1.5 for intermediate
    :param ftp_target: [1], FTP index, 1 for thermally thin; 2 for thermally thick; 1.5 for intermediate
    :return:
    """
    fuel_density *= 1e6

    # prepare time array
    t_arr = np.arange(0, t_end * 60 + t_step / 2., t_step, dtype=float)

    # calculate emitter temperature
    emitter_temperature = param_fire(
        t=t_arr,
        A_t=room_total_surface_area,
        A_f=room_floor_area,
        A_v=vent_width * vent_height,
        h_eq=vent_height,
        q_fd=fuel_density,
        lbd=k,
        rho=rho,
        c=c,
        t_lim=t_lim
    )

    # calculate the view factor between the emitter and receiver
    phi = phi_parallel_any_br187(
        emitter_width,
        emitter_height,
        emitter_width / 2,
        emitter_height / 2,
        emitter_receiver_separation
    )

    # calculate imposed (radiation) heat flux at the receiver
    hf = phi * 5.78e-8 * emissivity * (emitter_temperature ** 4 - 293.15 ** 4)

    # calculate flux-time product
    ftp = np.zeros_like(t_arr)
    ftp_i = ((((hf[:-1] + hf[1:]) * 1e-3) * 0.5 - chf * 1e-3) ** ftp_index) * (t_arr[1:] - t_arr[:-1])
    ftp[1:] = np.cumsum(np.where(ftp_i < 0, 0., ftp_i))

    try:
        if ftp_target <= np.amax(ftp):
            t_ig = t_arr[np.argmin(np.abs(ftp - ftp_target))]
        else:
            t_ig = np.inf
    except:
        t_ig = np.nan

    return t_ig, phi, ftp[-1]


def calc_receiver_temperature(
        fire_type: int, t: np.ndarray, A_t: float, A_f: float, A_v: float, h_eq: float, q_fd: float, lbd: float,
        rho: float, c: float, t_lim: float,
        emitter_width: float,
        emitter_height: float,
        emitter_receiver_separation: float,
        emissivity: float,
) -> np.ndarray:
    if fire_type == 0:
        temperature = param_fire(
            t=t, A_t=A_t, A_f=A_f, A_v=A_v, h_eq=h_eq, q_fd=q_fd, lbd=lbd, rho=rho, c=c, t_lim=t_lim
        )

        # calculate the view factor between the emitter and receiver
        phi = phi_parallel_any_br187(
            emitter_width,
            emitter_height,
            emitter_width / 2,
            emitter_height / 2,
            emitter_receiver_separation
        )

        # calculate imposed (radiation) heat flux at the receiver
        hf = phi * 5.78e-8 * emissivity * (temperature ** 4 - 293.15 ** 4)

    elif fire_type == 1:
        # todo
        # radiative_heat_flux_from_fire(
        #     t=t,
        #     fire_hrr_density_kWm2=0,
        # )
        temperature = None
    else:
        raise ValueError(f'Unknown `fire_type` {fire_type}')

    return temperature
