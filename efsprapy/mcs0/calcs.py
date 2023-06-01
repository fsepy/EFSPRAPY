__all__ = (
    'main',
)

import os
import tempfile
from typing import Optional

import numpy as np
from fsetools.lib.fse_bs_en_1991_1_2_parametric_fire import temperature as param_temperature
from fsetools.lib.fse_thermal_radiation import phi_parallel_any_br187

from efsprapy.func.controlled_fire import heat_flux_from_controlled_fire
from efsprapy.safir.test_files import therm1d_hf_50
from efsprapy.safir.therm2d import Run, PPXML


def main(
        t_end: float,
        t_step: float,
        vent_width: float,
        vent_height: float,
        room_height: float,
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
        fire_mode: int,
        fire_hrr_density_kWm2: float,
        fire_alpha: float,
        detector_to_fire_vertical_distance: float,
        detector_act_temp: float,
        detector_to_fire_horizontal_distance: float,
        detector_response_time_index: float,
        detector_conduction_factor: float,
        T_ig: float,
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

    # calculate incident heat flux at the receiver
    q_1, q_2 = calculate_incident_heat_flux(
        t=t_arr,
        fire_mode=fire_mode,
        A_t=room_total_surface_area,
        A_f=room_floor_area,
        A_v=vent_width * vent_height,
        h_eq=vent_height,
        q_fd=fuel_density,
        lbd=k,
        rho=rho,
        c=c,
        t_lim=t_lim,
        fire_hrr_density_kWm2=fire_hrr_density_kWm2,
        fire_alpha=fire_alpha,
        H=room_height,
        S=emitter_receiver_separation,
        W_o=vent_width,
        detector_to_fire_vertical_distance=detector_to_fire_vertical_distance,
        detector_act_temp=detector_act_temp,
        detector_to_fire_horizontal_distance=detector_to_fire_horizontal_distance,
        detector_response_time_index=detector_response_time_index,
        detector_conduction_factor=detector_conduction_factor,
    )

    # calculate the view factors
    if fire_mode == 0:
        phi_1 = phi_parallel_any_br187(
            emitter_width,
            emitter_height,
            emitter_width / 2,
            emitter_height / 2,
            emitter_receiver_separation
        )
        phi_2 = np.nan
    elif fire_mode == 1:
        phi_1 = phi_parallel_any_br187(
            emitter_width,
            emitter_height,
            emitter_width / 2,
            emitter_height / 2,
            emitter_receiver_separation
        )
        phi_2 = phi_1
    else:
        raise ValueError('Unknown `fire_mode`')

    t_ig_2, solver_iter = IgnitionTimeSolverIncidentHeatFluxSafir(
        t=t_arr,
        q_1=q_1,
        phi_1=phi_1,
        epsilon_1=1.,
        q_2=q_2,
        phi_2=phi_2,
        epsilon_2=1.,

    ).solve(T_ig=T_ig)

    # calculate flux-time product
    t_ig_1, ftp = calculate_t_ig_ftp(
        t=t_arr,
        q_1=q_1,
        phi_1=phi_1,
        emissivity=emissivity,
        chf=chf,
        ftp_index=ftp_index,
        ftp_target=ftp_target,
        q_2=q_2,
        phi_2=phi_2,
    )

    return t_ig_1, t_ig_2, phi_1, phi_2, ftp[-1], solver_iter,


def calculate_t_ig_ftp(t, q_1, phi_1, emissivity, chf, ftp_index, ftp_target, q_2=None, phi_2=None):
    hf = phi_1 * (q_1 - 5.67e-8 * emissivity * 293.15 ** 4)
    if q_2 is not None:
        hf += phi_2 * (q_1 - 5.67e-8 * emissivity * 293.15 ** 4)
    ftp = np.zeros_like(t)
    ftp_i_diff = (hf[:-1] + hf[1:]) * 0.5
    ftp_i_diff[ftp_i_diff < chf] = chf
    ftp_i = ((ftp_i_diff * 1e-3 - chf * 1e-3) ** ftp_index) * (t[1:] - t[:-1])
    ftp[1:] = np.cumsum(ftp_i)
    try:
        if ftp_target <= np.amax(ftp):
            t_ig = t[np.argmin(np.abs(ftp - ftp_target))]
        else:
            t_ig = np.inf
    except:
        t_ig = np.nan
    return t_ig, ftp


def calculate_incident_heat_flux(
        t: np.ndarray,
        fire_mode: int,
        A_t: float,
        A_f: float,
        A_v: float,
        h_eq: float,
        q_fd: float,
        lbd: float,
        rho: float,
        c: float,
        t_lim: float,
        fire_hrr_density_kWm2: float,
        fire_alpha: float,
        H: float,
        S: float,
        W_o: float,
        detector_to_fire_vertical_distance: float,
        detector_act_temp: float,
        detector_to_fire_horizontal_distance: float,
        detector_response_time_index: float,
        detector_conduction_factor: float,
        epsilon_1: float = .99
):
    if fire_mode == 0:
        T = param_temperature(
            A_t=A_t, A_f=A_f, A_v=A_v, h_eq=h_eq, q_fd=q_fd, lbd=lbd, rho=rho, c=c, t_lim=t_lim, t=t,
        )
        q_1 = 5.67e-8 * epsilon_1 * T ** 4
        q_2 = None
    elif fire_mode == 1:
        q_1, q_2 = heat_flux_from_controlled_fire(
            t=t,
            fire_hrr_density_kWm2=fire_hrr_density_kWm2,
            fire_alpha=fire_alpha,
            H=H,
            S=S,
            W_o=W_o,
            detector_to_fire_vertical_distance=detector_to_fire_vertical_distance,
            detector_act_temp=detector_act_temp,
            detector_to_fire_horizontal_distance=detector_to_fire_horizontal_distance,
            detector_response_time_index=detector_response_time_index,
            detector_conduction_factor=detector_conduction_factor,
        )
    else:
        raise ValueError(f'Unknown `fire_mode')
    return q_1, q_2


class IgnitionTimeSolverIncidentHeatFluxSafir:
    def __init__(
            self,
            t: np.ndarray,
            q_1: np.ndarray,
            phi_1: float,
            epsilon_1: float,
            q_2: Optional[np.ndarray] = None,
            phi_2: Optional[float] = None,
            epsilon_2: Optional[float] = None,
            safir_input: Optional[str] = None,
            dir_work: Optional[str] = None,
    ):
        if safir_input is None:
            with open(therm1d_hf_50, 'r') as f:
                self.safir_in_0_s = f.read()
        else:
            self.safir_in_0_s = safir_input

        if dir_work is not None:
            assert os.path.isdir(dir_work) and os.path.exists(dir_work)
        else:
            self.__dir_work = dir_work

        assert t.shape == q_1.shape
        self.t = t
        self.q_1 = q_1
        self.phi_1 = phi_1
        self.epsilon_1 = epsilon_1
        self.q_2 = q_2
        self.phi_2 = phi_2
        self.epsilon_2 = epsilon_2

    def solve(
            self,
            T_ig: float,
            solver_iter_cap: int = 20,
            solver_t_ig_tol: float = 30,
    ):
        if self.__dir_work is not None:
            return self.__solve(T_ig=T_ig, solver_iter_cap=solver_iter_cap, solver_t_ig_tol=solver_t_ig_tol)
        else:
            with tempfile.TemporaryDirectory() as dir_work:
                print(dir_work)
                self.__dir_work = dir_work
                return self.__solve(T_ig=T_ig, solver_iter_cap=solver_iter_cap, solver_t_ig_tol=solver_t_ig_tol)

    def __solve(
            self,
            T_ig: float,
            solver_iter_cap: int = 20,
            solver_t_ig_tol: float = 30,
    ):
        t_step = solver_t_ig_tol
        t_end = np.amax(self.t[self.q_1 > np.amin(self.q_1)])

        # initial heat flux
        hf = self.phi_1 * (self.q_1 - self.epsilon_1 * 5.67e-8 * np.full_like(self.t, 293.15) ** 4)
        if self.q_2 is not None:
            hf += self.phi_2 * (self.q_2 - self.epsilon_2 * 5.67e-8 * np.full_like(self.t, 293.15) ** 4)
        hf[hf < 0] = 0

        t_ig_arr = list()
        n_iter = 1
        t_ig = -1
        non_ignition_streak: int = 0
        while n_iter < solver_iter_cap:
            print(f'Iteration {n_iter}', end=': ')
            T_1 = self.__safir_receiver_temperature(n_iter=0, t=self.t, t_end=t_end, t_step=t_step, hf=hf)
            if T_1 is not None:
                t_ig = self.__get_ignition_time(self.t, T_1, T_ig)
            else:
                print('Fail')
                return np.nan, n_iter

            if t_ig == np.inf:
                if n_iter == 1:
                    return t_ig, n_iter
                if non_ignition_streak >= 1:
                    return np.inf, n_iter
                non_ignition_streak += 1
                print('N/A')
            else:
                t_ig_arr.append(t_ig)
                non_ignition_streak = 0
                print(t_ig)

            # check if convergence is sought
            if len(t_ig_arr) >= 2 and abs(t_ig_arr[-1] - t_ig_arr[-2]) <= solver_t_ig_tol:
                t_ig = t_ig_arr[-1]
                break

            hf = self.phi_1 * (self.q_1 - self.epsilon_1 * 5.67e-8 * T_1 ** 4)
            if self.q_2 is not None:
                hf += self.phi_2 * (self.q_2 - self.epsilon_2 * 5.67e-8 * T_1 ** 4)
            hf[hf < 0] = 0

            n_iter += 1

        return t_ig, n_iter

    def __safir_receiver_temperature(
            self,
            n_iter: int,
            t: np.ndarray,
            hf: np.ndarray,
            t_step: float,
            t_end: float,
    ):
        """"""

        fn_bc = f'{n_iter}_hf'
        with open(os.path.join(self.__dir_work, fn_bc), 'w+') as f:
            f.write('\n'.join(f'{t[i]:g}, {hf[i]:.3f}' for i in range(len(t))))

        fn = f'i{n_iter}'
        fn_safir_in = f'{fn}.in'
        fp_safir_in = os.path.join(self.__dir_work, fn_safir_in)
        with open(fp_safir_in, 'w+') as f:
            f.write(self.safir_in_0_s.format(
                fn_bc=fn_bc,
                materials=f'WOODEC5\n    450. 0 25 9 0.8 1.2 0.0 0.0 1.0',
                t_step=t_step,
                t_end=t_end,
            ))

        Run().run(fp_safir_in)

        try:
            with open(os.path.join(self.__dir_work, f'{fn}.XML')) as f:
                pp = PPXML(xml=f.read())

            return np.interp(t, pp.t, pp.get_nodes_temp(nodes=np.array([1, ]))[0, :] + 273.15)

        except:
            return None

    @staticmethod
    def __get_ignition_time(t, T, T_ig):
        try:
            return np.amin(t[T >= T_ig])
        except ValueError:
            return np.inf
