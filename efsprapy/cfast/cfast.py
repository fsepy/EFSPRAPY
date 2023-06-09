import logging
import shutil
import subprocess
import time
from os import path, environ, getcwd
from subprocess import Popen, PIPE
from time import time
from typing import Callable

import numpy as np

try:
    from subprocess import CREATE_NO_WINDOW
except ImportError:
    CREATE_NO_WINDOW = 0

logger = logging.getLogger('gui')


def detect_binary(fp=None):
    from os import name as platform_name
    if platform_name == 'nt':
        executable = "cfast.exe"
    else:
        raise NotImplementedError

    env_home = environ.get("cfast")
    env_path = shutil.which(executable)

    possible_paths = [
        fp,
        path.join(getcwd(), executable),
        path.join('C:', 'Program', 'Files', 'firemodels', 'cfast7', executable),
        path.join(env_home, executable) if env_home is not None else None,
        path.join(env_path) if env_path is not None else None
    ]

    for fp_ in possible_paths:
        if fp_ is not None and path.isfile(fp_):
            return fp_

    raise FileNotFoundError(f"Unable to find {executable}")


class Run:
    FP_CFAST_EXE = detect_binary()

    def __init__(self):
        self.__fp_in = None
        self.__fp_cfast_exe = Run.FP_CFAST_EXE if Run.FP_CFAST_EXE else detect_binary()
        self.__run_stdout = ''

    def get_stdout(self):
        return self.__run_stdout

    def run(self, fp: str, timeout: int = 1800, print_time: Callable = None):
        self.set_fp_in(fp)
        self.__run_stdout = self.__run_worker(
            exe=self.__fp_cfast_exe,
            fp_in=self.__fp_in,
            fp_stdout=f'{path.splitext(self.__fp_in)[0]}.stdout',
            print_time=print_time,
            timeout=timeout
        )
        return self

    def read_outputs(self):
        fn = path.basename(self.__fp_in)
        case_name = path.splitext(fn)[0]
        dir_name = path.dirname(self.__fp_in)

        _ = np.genfromtxt(path.join(dir_name, f'{case_name}_compartments.csv'), delimiter=',', skip_header=4)
        t = _[:, 0]
        upper_layer_temperature = _[:, 1]
        lower_layer_temperature = _[:, 2]
        layer_height = _[:, 3]
        hrr_actual = _[:, 32]

        _ = np.genfromtxt(path.join(dir_name, f'{case_name}_devices.csv'), delimiter=',', skip_header=4)
        sprinkler_temperature = _[:, 1]

        return t, upper_layer_temperature, lower_layer_temperature, layer_height, hrr_actual, sprinkler_temperature

    def set_fp_in(self, fp_in: str):
        fp_in = path.realpath(fp_in)
        if path.isfile(fp_in):
            self.__fp_in = fp_in
        else:
            raise FileNotFoundError(f'File does not exist {fp_in}')

    @staticmethod
    def __run_worker(exe, fp_in, timeout: int = 5, fp_stdout: str = None, print_time: Callable = None):
        """"""
        fn_in = path.basename(fp_in)
        if fn_in.endswith('.in'):
            fn_in = fn_in[:-3]

        proc = Popen(
            args=f'{exe} {fn_in} -O:CD',
            stdout=PIPE,
            stdin=PIPE,
            stderr=PIPE,
            creationflags=CREATE_NO_WINDOW,
            cwd=path.dirname(fp_in),
            universal_newlines=True,
            encoding='utf-8',
        )

        lines = list()
        time_0 = time()
        while time() - time_0 < timeout:
            line = proc.stdout.readline()
            if line == '' and proc.poll() is not None:
                break
            if line:
                if fp_stdout is not None:
                    lines.append(line.strip())
        try:
            proc.communicate(timeout=0.1)
        except subprocess.TimeoutExpired as e:
            if print_time:
                print_time(0)
            lines.append(f'Process timed out after {timeout:g} seconds.')
            proc.kill()
            proc.communicate()
            raise e
        finally:
            if fp_stdout is not None:
                try:
                    with open(fp_stdout, 'w+') as f:
                        f.write('\n'.join(lines))
                except Exception as e2:
                    logger.warning(f'Unable to write to {fp_stdout}. {type(e2).__name__}.')
        return '\n'.join(lines)
