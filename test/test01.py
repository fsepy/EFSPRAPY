import os

import matplotlib.pyplot as plt
import numpy as np
from fsetools.lib.fse_bs_en_1991_1_2_parametric_fire import temperature as param_temperature

from efsprapy.safir.test_files import therm1d_hf_50 as fp_safir_in
from efsprapy.safir.therm2d import Run, PPXML

# fig, axes = plt.subplots(ncols=2, nrows=10)
fig, ax = plt.subplots()


class IgnitionTimeSolver:
    def __init__(self, fp_safir_input: str, dir_work: str):
        assert os.path.isfile(fp_safir_input) and os.path.exists(fp_safir_input)
        self.__fp_safir_input_0 = fp_safir_input
        with open(self.__fp_safir_input_0, 'r') as f:
            self.safir_in_0_s = f.read()

        assert os.path.isdir(dir_work) and os.path.exists(dir_work)
        self.__dir_work = dir_work

        self.__safir_runner = Run()
        self.__p = PPXML()

        self.__safir_input_kwargs = dict()

    def run(self, t_2: np.ndarray, T_2: np.ndarray, T_ig: float, solver_iter_cap: int = 10,
            solver_t_ig_tol: float = 10):
        epsilon = 1.0
        sigma = 5.67e-8
        phi = 0.1
        T_1 = np.full_like(t_2, 293.15)
        t_end = np.amax(t_2[T_2 > np.amin(T_2)])
        solver_t_ig_arr = np.full((solver_iter_cap,), np.inf, float)

        T_1 = self.__run_single_worker(
            n_iter=0, t=t_2, T_1=T_1, T_2=T_2, t_end=t_end, epsilon=epsilon, sigma=sigma, phi=phi,
        )
        if np.amax(T_1) > T_ig:
            solver_t_ig_arr[0] = np.amin(t_2[T_1 >= T_ig])
        else:
            return np.inf

        T_1_ = T_1.copy()

        n_iter_count = 1
        while n_iter_count < solver_iter_cap:
            print(f'running {n_iter_count}')
            T_1_ave = (T_1 + T_1_) * 0.5
            T_1_ = T_1
            T_1 = self.__run_single_worker(
                n_iter=n_iter_count, t=t_2, T_1=T_1_ave, T_2=T_2, t_end=t_end, epsilon=epsilon, sigma=sigma, phi=phi,
            )

            # note ignition time
            T_1_max = np.amax(T_1)
            if n_iter_count > 0 and T_1_max > T_ig:
                solver_t_ig_arr[n_iter_count] = np.min(t_2[T_1 >= T_ig])
            else:
                solver_t_ig_arr[n_iter_count] = np.inf
            print(solver_t_ig_arr)

            # update heat flux
            ax.plot(t_2, T_1, label=f'{n_iter_count}')

            n_iter_count += 1

        ax.legend().set_visible(True)

        plt.show()

    def __run_single_worker(
            self, n_iter: int, t: np.ndarray, T_1: np.ndarray, T_2: np.ndarray, epsilon: float,
            sigma: float, phi: float, t_end: float
    ):
        dir_work = self.__dir_work

        hf = phi * epsilon * sigma * (T_2 ** 4 - T_1 ** 4)
        hf[hf < 0] = 0

        fn_bc = f'{n_iter}_hf'
        with open(os.path.join(dir_work, fn_bc), 'w+') as f:
            f.write('\n'.join(f'{t[i]:g}, {hf[i]:.3f}' for i in range(len(t))))

        fn = f'{n_iter}_{os.path.splitext(os.path.basename(self.__fp_safir_input_0))[0]}'
        fn_safir_in = f'{fn}.in'
        fp_safir_in = os.path.join(dir_work, fn_safir_in)
        with open(fp_safir_in, 'w+') as f:
            f.write(self.safir_in_0_s.format(
                fn_bc=fn_bc,
                materials=f'WOODEC5\n    450. 0 25 9 0.8 1.2 0.0 0.0 1.0',
                t_step=5,
                t_end=t_end,
            ))

        self.__safir_runner.run(fp_safir_in)

        with open(os.path.join(dir_work, f'{fn}.XML')) as f:
            self.__p.xml = f.read()
        T_1 = np.interp(t, self.__p.t, self.__p.get_nodes_temp(nodes=np.array([1, ]))[0, :] + 273.15)

        return T_1

    def __update_safir_input_kwargs(self, fn_flux_bc: str, time_step: int, time_end: int):
        self.__safir_input_kwargs = dict(fn_flux_bc=fn_flux_bc, time_step=time_step, time_end=time_end)


if __name__ == '__main__':
    run = IgnitionTimeSolver(fp_safir_input=fp_safir_in, dir_work=r'C:\Users\IanFu\Desktop\New folder (2)')

    _t_2_ = np.arange(0, 10800 + 0.5, 10)
    _T_2_ = param_temperature(
        A_t=360, A_f=100, A_v=36.1, h_eq=1, q_fd=600e6, lbd=1, rho=1, c=2250000, t_lim=20 * 60, t=_t_2_, T_0=293.15
    )

    # prepare boundary condition
    run.run(t_2=_t_2_, T_2=_T_2_, T_ig=300 + 273.15)
