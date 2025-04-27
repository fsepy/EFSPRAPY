import shutil
import subprocess
import time
from os import path, environ, getcwd
from subprocess import Popen, PIPE
from tempfile import TemporaryDirectory
from time import time
from typing import Callable, Optional, Tuple

import numpy as np

from ..cfast.test_files import simple

try:
    from subprocess import CREATE_NO_WINDOW
except ImportError:
    CREATE_NO_WINDOW = 0


class CFastError(Exception):
    """Base class for CFAST-related errors."""
    pass


def detect_binary(fp=None):
    """
    Detect CFAST binary path.

    Args:
        fp: Optional explicit path to binary

    Returns:
        Path to CFAST binary

    Raises:
        FileNotFoundError: If binary cannot be found
    """
    from os import name as platform_name
    if platform_name == 'nt':
        executable = "cfast.exe"
    else:
        # For Linux/Mac support, you'd need to add the correct executable name
        raise NotImplementedError("Only Windows is currently supported")

    env_home = environ.get("cfast7.7.3")
    env_path = shutil.which(executable)

    possible_paths = [
        fp,
        path.join(getcwd(), executable),
        path.join('C:', 'Program Files', 'firemodels', 'cfast7', executable),
        path.join(env_home, executable) if env_home is not None else None,
        env_path if env_path is not None else None
    ]

    for fp_ in possible_paths:
        if fp_ is not None and path.isfile(fp_):
            return fp_

    raise FileNotFoundError(f"Unable to find {executable}")


class Run:
    """Class to run CFAST simulations with improved support for multiprocessing."""

    FP_BIN = None  # Will be set on first initialization

    def __init__(self, binary_path: Optional[str] = None):
        """
        Initialize a CFAST runner.

        Args:
            binary_path: Optional explicit path to CFAST binary
        """
        self.__fp_in = None

        # Initialize the binary path only once for all instances
        if Run.FP_BIN is None:
            try:
                Run.FP_BIN = binary_path if binary_path else detect_binary(
                    path.join(path.dirname(__file__), 'cfast7.7.3')
                )
            except FileNotFoundError:
                Run.FP_BIN = detect_binary()

        self.__fp_cfast_exe = Run.FP_BIN
        self.__run_stdout = ''

    def get_stdout(self) -> str:
        """Get the standard output from the last run."""
        return self.__run_stdout

    def run(self, fp: str, timeout: int = 1800, print_time: Callable = None) -> 'Run':
        """
        Run CFAST simulation.

        Args:
            fp: Path to input file
            timeout: Timeout in seconds
            print_time: Optional callback for time updates

        Returns:
            Self for method chaining
        """
        self.set_fp_in(fp)
        try:
            self.__run_stdout = self.__run_worker(
                exe=self.__fp_cfast_exe,
                fp_in=self.__fp_in,
                fp_stdout=f'{path.splitext(self.__fp_in)[0]}.stdout',
                print_time=print_time,
                timeout=timeout
            )
            return self
        except subprocess.TimeoutExpired:
            # Silently handle timeout
            if print_time:
                print_time(0)
            return self
        except Exception:
            # Silently handle other errors
            return self

    def read_outputs(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """
        Read CFAST output files.

        Returns:
            Tuple containing:
            - time array
            - upper layer temperature
            - lower layer temperature
            - layer height
            - actual heat release rate
            - sprinkler temperature
        """
        fn = path.basename(self.__fp_in)
        case_name = path.splitext(fn)[0]
        dir_name = path.dirname(self.__fp_in)

        compartments_file = path.join(dir_name, f'{case_name}_compartments.csv')
        devices_file = path.join(dir_name, f'{case_name}_devices.csv')

        # Read compartments data
        compartment_data = np.genfromtxt(compartments_file, delimiter=',', skip_header=4)
        t = compartment_data[:, 0]
        upper_layer_temperature = compartment_data[:, 1]
        lower_layer_temperature = compartment_data[:, 2]
        layer_height = compartment_data[:, 3]
        hrr_actual = compartment_data[:, 32]

        # Read devices data
        device_data = np.genfromtxt(devices_file, delimiter=',', skip_header=4)
        sprinkler_temperature = device_data[:, 1]

        return t, upper_layer_temperature, lower_layer_temperature, layer_height, hrr_actual, sprinkler_temperature

    def set_fp_in(self, fp_in: str) -> None:
        """
        Set input file path.

        Args:
            fp_in: Path to input file
        """
        fp_in = path.realpath(fp_in)
        if path.isfile(fp_in):
            self.__fp_in = fp_in
        else:
            raise FileNotFoundError(f'File does not exist {fp_in}')

    @staticmethod
    def __run_worker(exe: str, fp_in: str, timeout: int = 5, fp_stdout: str = None,
                     print_time: Callable = None) -> str:
        """
        Run CFAST process with improved handling for multiprocessing.

        Args:
            exe: Path to executable
            fp_in: Path to input file
            timeout: Timeout in seconds
            fp_stdout: Path to save stdout
            print_time: Optional callback for time updates

        Returns:
            Standard output as string
        """
        fn_in = path.basename(fp_in)
        if fn_in.endswith('.in'):
            fn_in = fn_in[:-3]

        cwd = path.dirname(fp_in)
        command = f'{exe} {fn_in} -O:CD'

        proc = Popen(
            args=command,
            stdout=PIPE,
            stdin=PIPE,
            stderr=PIPE,
            creationflags=CREATE_NO_WINDOW,
            cwd=cwd,
            universal_newlines=True,
            encoding='utf-8',
        )

        lines = []
        start_time = time()
        try:
            while time() - start_time < timeout:
                line = proc.stdout.readline()
                if line == '' and proc.poll() is not None:
                    break
                if line:
                    lines.append(line.strip())

            # Check if process is still running after the loop
            if proc.poll() is None:
                # Process timed out - kill it
                if print_time:
                    print_time(0)
                lines.append(f'Process timed out after {timeout:g} seconds.')
                proc.kill()
                proc.communicate()
                raise subprocess.TimeoutExpired(command, timeout)

        except subprocess.TimeoutExpired:
            if print_time:
                print_time(0)
            lines.append(f'Process timed out after {timeout:g} seconds.')
            proc.kill()
            proc.communicate()
            raise
        except Exception:
            proc.kill()
            proc.communicate()
            raise
        finally:
            if fp_stdout is not None:
                try:
                    with open(fp_stdout, 'w+') as f:
                        f.write('\n'.join(lines))
                except Exception:
                    pass

        return '\n'.join(lines)


def calculate_hrr_and_smoke_with_sprinkler_suppression_cfast(
        t_arr: np.ndarray,
        fire_mode: int,

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

    if fire_mode == 1:
        t_end = min(7200., t_end)
    if fire_mode == 2:
        t_end = min(3600., t_end)

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
        fp_cfast_in = path.join(dir_work, fn_cfast_in)
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
        _ = Run().run(fp_cfast_in)

        t_, ult_, llt_, lh_, hrr_, spt_ = _.read_outputs()

    ult_ = np.interp(t_arr, t_, ult_) + 273.15
    # llt_ = np.interp(t_arr, t_, llt_)
    # lh_ = np.interp(t_arr, t_, lh_)
    hrr_ = np.interp(t_arr, t_, hrr_)
    t_d = np.interp(T_act - 273.15, spt_, t_)

    return hrr_, ult_, t_d
