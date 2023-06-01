import logging
import re
import shutil
import subprocess
import time
from os import path, environ, getcwd, devnull
from subprocess import Popen, PIPE
from time import time
from typing import List, Dict, Callable, Union, Optional

import numpy as np

from .re_patterns import stdout_last_time

try:
    from subprocess import CREATE_NO_WINDOW
except ImportError:
    CREATE_NO_WINDOW = 0

logger = logging.getLogger('gui')

__all__ = 'PPXML', 'Run', 'safir_single_run', 'batch_run'


class PPXML:
    """
    SAFIR Therm2D Post Processor (XML)
    """

    def __init__(self, xml: Optional[str] = None):
        self.__data: Union[None, dict] = None  # dict parsed from xml
        self.__ns: Union[None, np.ndarray] = None  # node indexes where temperatures are measured `ns=[n1, n2, ...]`
        self.__xs: Union[None, np.ndarray] = None  # x coordinates [x1, x2, ...]
        self.__ys: Union[None, np.ndarray] = None  # y coordinates [y1, y2, ...]
        self.__xys: Union[None, np.ndarray] = None  # np.array(list(zip(x, y)))
        self.__ts: Union[None, np.ndarray] = None  # time steps [30, 60, ...]
        self.__Ts: Union[None, np.ndarray] = None  # temperatures [T1, T2, ...], T1 = [T11, T12, ...]
        self.__xml = None  # xml raw
        self.__xml_changed = True  # whether `self.__xml` has changed
        self.__ns_xys_changed = True  # whether `self.__ns` or `self.__xys` haved changed
        self.__func_xy_T: dict = dict()  # time: interp2d temperature at a given x, y

        self.xml = xml

    def get_nodes_temp(self, nodes: np.ndarray) -> np.ndarray:
        self.process_xml()
        temps = self.nodes2temp(Ts=self.__Ts, nodes=np.array(nodes))
        return temps

    @staticmethod
    def __get_line_xy(x1: float, y1: float, x2: float, y2: float, n: int, inclusive: bool = True) -> np.ndarray:
        """Get an array of coordinates [[xi1, yi1], [xi2, yi2], ...] evenly spaced between vectors (x1, y1) and
        (x2, y2). The separation `di=((xi2-xi1)**2+(yi2-yi1)**2)**0.5` is less or equal to `d`.

        :param x1:
        :param y1:
        :param x2:
        :param y2:
        :param n:
        :param inclusive:
        :return:
        """

        # make [[xi1, yi1], [xi2, yi2], ...]
        arr = np.stack((np.linspace(x1, x2, n), np.linspace(y1, y2, n)), axis=1)

        # exclude first and last item in `arr` if `inclusive` is True
        if not inclusive:
            return (arr[1:] + arr[:-1]) / 2
        return arr

    def __get_line_temp_single_t(
            self, t_i: int, x1: float, y1: float, x2: float, y2: float, n: int, inclusive: bool = True
    ) -> np.ndarray:
        from scipy.interpolate import LinearNDInterpolator
        """To get temperature along a line (x1, y1) -> (x2, y2) at a given index corresponding time step.

        Make an interpolation function that interpolate temperature `T` for given location(s) `x, y`. This interp
        function is time-step-index `i` specific and stores instantiated interp object to a dict.

        :param t_i:
        :param x1:
        :param y1:
        :param x2:
        :param y2:
        :param n:
        :param inclusive:
        :return:
        """

        # Make an interpolation function that interpolate temperature `T` for given location(s) `x, y`. This interp
        # function is time-step-index `i` specific and stores instantiated interp object to a dict.
        if t_i in self.__func_xy_T:
            l_interp = self.__func_xy_T[t_i]
        else:
            l_interp = LinearNDInterpolator(np.column_stack((self.x, self.y)), self.T[t_i, :])
            self.__func_xy_T[t_i] = l_interp

        xy = self.__get_line_xy(x1=x1, y1=y1, x2=x2, y2=y2, n=n, inclusive=inclusive)
        return l_interp(xy)

    def __get_line_isotherm_single_t(self, t_i, T_iso, x1, y1, x2, y2, n):
        """Make an interpolation function `f(x, y)=T` at a given time index `t_i` that yield temperature at a given
        x, y location."""
        from scipy.interpolate import interp1d

        # workout number of sampling points based on resolution `r`
        line_xy = self.__get_line_xy(x1=x1, x2=x2, y1=y1, y2=y2, n=n)
        line_length = np.linspace(0, ((x2 - x1) ** 2 + (y2 - y1) ** 2) ** 0.5, line_xy.shape[0])

        temp = self.__get_line_temp_single_t(t_i=t_i, x1=x1, y1=y1, x2=x2, y2=y2, n=n)

        if np.amax(temp) < T_iso:
            return None, None, np.amin(line_length)
        if np.amin(temp) > T_iso:
            return None, None, np.amax(line_length)

        isotherm = interp1d(temp, line_length)(T_iso)
        return None, None, isotherm

    def get_line_temp_isotherm(self, T_iso, x1, y1, x2, y2, n):
        self.process_xml()
        isotherm_depth = np.zeros_like(self.__ts)
        for i, t_ in enumerate(self.__ts):
            isotherm_depth[i] = self.__get_line_isotherm_single_t(i, T_iso, x1, y1, x2, y2, n)[-1]
        return isotherm_depth

    def get_line_temp(self, x1, y1, x2, y2, n, inclusive: bool = True):
        """Produce a mxn array with m represents the depth from (x1, y1) -> (x2, y2) and n represents time step"""
        self.process_xml()
        self.__get_line_xy(x1=x1, y1=y1, x2=x2, y2=y2, n=n)
        temp = np.zeros(
            shape=(self.__get_line_xy(x1=x1, y1=y1, x2=x2, y2=y2, n=n).shape[0], self.__ts.shape[0]),
            dtype=float
        )
        for i, t_ in enumerate(self.__ts):
            temp[:, i] = self.__get_line_temp_single_t(t_i=i, x1=x1, y1=y1, x2=x2, y2=y2, n=n, inclusive=inclusive)
        return temp

    def get_nodes_temp_ave(
            self,
            ns: np.ndarray = None,
            ws: np.ndarray = None,
            mode: str = 'ws',
            temp_ns: np.ndarray = None,
            fp_save_num: str = None,
    ):
        """
        Get average temperature based on the defined nodes :code:`ns` or temperatures :code:`T_ns`.

        :param ns:              A list/tuple of integers describing node index. This is not required if :code:`T_ns` is
                                provided.
                                :code:`ns = [n1, n2, ...]`
        :param ws:              A list/tuple of floats describing weighting of the defined nodes when averaging
                                temperatures.
                                :code:`ws = [w1, w2, ...]`
        :param mode:            `ws` for weighted average, k_y_theta for
        :param temp_ns:            An array of temperatures. This will be derived from the xml if not provided.
                                :code:`T_ns = [T_n1, T_n2, ...]`
                                :code:`T_n1 = [float, float, ...]`
        :param fp_save_num:
        :return:                [C] Averaged temperatures

        The below condition should be satisfied.

        :code:`len(ns) == len(ws) == T_ns.shape[0]`
        """
        self.process_xml()

        if temp_ns is None:
            if ns is None:
                raise ValueError('Either `nodes` or `temperature_nodes` must be defined')
            temp_ns = self.get_nodes_temp(ns)

        if mode == 'ws':
            T_ave = self.__T_ave(Ts=temp_ns, weights=ws)
        elif mode == 'k_y_theta':
            T_ave = self.__T_ave_k_y_theta(Ts=temp_ns, weights=ws)
        else:
            raise ValueError(f'`mode` can be either `ws` or `k_y_theta`, `{mode}` is provided')

        if fp_save_num is not None:
            try:
                np.savetxt(
                    fname=fp_save_num,
                    X=np.vstack((self.__ts, temp_ns, T_ave)).transpose(),
                    fmt='%.3f',
                    header=f'TIME,' + ','.join([f'node {n}' for n in ns]) + ',MEAN',
                    delimiter=',',
                    comments='',
                )
            except Exception as e:
                logger.warning(f'Failed to numerical data to {fp_save_num}, {e}')

        return T_ave

    def get_nodes_from_xy(self, xys: Union[tuple, list]) -> np.ndarray:
        self.process_xml()
        return self.xys2nodes(self.__xs, self.__ys, xys)

    def get_xys_temp(self, xys: Union[tuple, list]):
        return self.get_nodes_temp(nodes=self.get_nodes_from_xy(xys))

    def get_xys_temp_ave(self, xys: Union[tuple, list], *args, **kwargs):
        ns = self.get_nodes_from_xy(xys=xys)
        return self.get_nodes_temp_ave(ns=ns, *args, **kwargs)

    def process_xml(self):
        # only process if the xml is changed
        if self.__xml_changed is True:
            # self.xml = self.xml  # assign xml
            self.__data = self.xml2dict(self.xml)  # convert xml to dict
            self.__xs, self.__ys = self.dict2xys(self.__data)  # obtain all node x and y coordinates
            self.__xys = np.column_stack((self.__xs, self.__ys))  # nx2 array concatenating x and y coors
            self.__ts, self.__Ts = self.dict2tsTs(self.__data)  # obtain time steps and temperatures
            self.__func_xy_T = dict()
            self.__xml_changed = False

    @property
    def t(self):
        self.process_xml()
        return self.__ts

    @property
    def T(self) -> np.ndarray:
        self.process_xml()
        return self.__Ts

    @property
    def x(self):
        self.process_xml()
        return self.__xs

    @property
    def y(self):
        self.process_xml()
        return self.__ys

    @property
    def xml(self):
        return self.__xml

    @xml.setter
    def xml(self, xml: str):
        self.__xml = xml
        self.__xml_changed = True

    @staticmethod
    def xml2dict(xml: str) -> dict:
        import xmltodict
        data_dict = xmltodict.parse(xml)
        return data_dict

    @staticmethod
    def dict2xys(data: dict) -> tuple:
        xys = data['SAFIR_RESULTS']['NODES']['N']
        xs, ys = list(), list()
        for xy in xys:
            xs.append(float(xy['P2']))
            ys.append(float(xy['P1']))
        return np.array(xs), np.array(ys)

    @staticmethod
    def dict2tsTs(data: dict) -> tuple:
        ts, Ts = list(), list()
        steps = data['SAFIR_RESULTS']['STEP']
        for step in steps:
            ts.append(step['TIME']['#text'])
            Ts.append(step['TEMPERATURES']['T'])
        ts.insert(0, 0)
        Ts.insert(0, [20] * len(Ts[0]))  # todo: initial temperature can be found in *.in safir input file
        # safir xml outputs last line repeating the same time step
        return np.array(ts, dtype=float)[:-1], np.array(Ts, dtype=float)[:-1]

    @staticmethod
    def xys2nodes(xs, ys, xys) -> np.ndarray:
        nodes = list()
        xys1 = np.stack((xs, ys), axis=0)
        for xy in xys:
            i = np.argmin(((xys1[0, :] - xy[0]) ** 2) + ((xys1[1, :] - xy[1]) ** 2))
            i += 1
            nodes.append(i)

        return np.array(nodes)

    @staticmethod
    def nodes2temp(Ts: np.ndarray, nodes: np.ndarray) -> np.ndarray:
        nodes_safir2py = [i - 1 for i in nodes]
        Ts_nodes = list()
        for node in nodes_safir2py:
            Ts_nodes.append([T[node] for T in Ts])
        return np.array(Ts_nodes)

    @staticmethod
    def __T_ave(Ts: np.ndarray, weights: np.ndarray) -> np.ndarray:
        # populate `weights` to match the shape of `Ts_nodes`
        Ts_weights = np.repeat(weights, Ts.shape[1])
        Ts_weights = np.reshape(Ts_weights, (Ts.shape[1], len(weights)), order='F')
        Ts_weights = Ts_weights.transpose()

        # calculate the summed weights for all time steps
        Ts_weights_sum = np.sum(Ts_weights, axis=0)

        # calculate normalised weights
        Ts_weights_normalised = np.zeros_like(Ts)
        for i in range(Ts_weights_normalised.shape[0]):
            Ts_weights_normalised[i, :] = np.divide(Ts_weights[i, :], Ts_weights_sum[:])

        return np.sum(Ts_weights_normalised * Ts, axis=0)

    @staticmethod
    def __T_ave_k_y_theta(Ts: np.ndarray, weights: np.ndarray) -> np.ndarray:
        from fsetools.libstd.bs_en_1993_1_2_2005_k_y_theta import (
            clause_3_2_1_1_k_y_theta_mod, clause_3_2_1_1_k_y_theta_mod_reversed
        )

        T_ave = np.zeros_like(Ts[0, :])
        for i in range(len(T_ave)):
            k_y_theta_i = clause_3_2_1_1_k_y_theta_mod(Ts[:, i] + 273.15)
            k_y_theta_mean = np.sum(np.multiply(k_y_theta_i, weights)) / np.sum(weights)
            T_ave_ = clause_3_2_1_1_k_y_theta_mod_reversed(k_y_theta_mean) - 273.15
            T_ave[i] = T_ave_

        return T_ave


class Run:
    def __init__(self):
        self.__fp_in = None
        self.__fp_safir_exe = self.detect_binary()

    def run(self, fp: str, timeout: int = 1800, print_time: Callable = None):
        self.set_fp_in(fp)
        return self.__run_worker(
            exe=self.__fp_safir_exe,
            fp_in=self.__fp_in,
            fp_stdout=f'{path.splitext(self.__fp_in)[0]}.stdout',
            print_time=print_time,
            timeout=timeout
        )

    def set_fp_in(self, fp_in: str):
        fp_in = path.realpath(fp_in)
        if path.isfile(fp_in):
            self.__fp_in = fp_in
        else:
            raise FileNotFoundError(f'File does not exist {fp_in}')

    @staticmethod
    def detect_binary(fp=None):
        from os import name as platform_name
        if platform_name == 'nt':
            executable = "safir.exe"
        else:
            raise NotImplementedError

        env_home = environ.get("safir")
        env_path = shutil.which(executable)

        possible_paths = [
            fp,
            path.join(getcwd(), executable),
            path.join("C:\work\\fem\\safir", executable),
            path.join(env_home, executable) if env_home is not None else None,
            path.join(env_path) if env_path is not None else None
        ]

        for fp_ in possible_paths:
            if fp_ is not None and path.isfile(fp_):
                return fp_

        raise FileNotFoundError(f"Unable to find {executable}")

    @staticmethod
    def __run_worker(exe, fp_in, timeout: int = 1800, fp_stdout: str = None, print_time: Callable = None):
        """"""
        fn_in = path.basename(fp_in)
        if fn_in.endswith('.in'):
            fn_in = fn_in[:-3]

        proc = Popen(
            args=f'{exe} {fn_in}',
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
                    print_time(float(re.findall(stdout_last_time, line.strip())[-1])) if print_time else None
                except Exception:
                    pass

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


class Therm2D(Run, PPXML):
    def __init__(self):
        super().__init__()


def safir_single_run(fp_safir_exe: str, fp_in: str, fp_stdout: str = None):
    dir_in = path.dirname(fp_in)
    fn_in = path.basename(fp_in)

    if fn_in.endswith('.in'):
        fn_in = fn_in[:-3]

    # construct command
    cmd = f'{fp_safir_exe} {fn_in.rstrip(".in")}'

    if fp_stdout is not None:
        subprocess.run(args=cmd, cwd=dir_in, stdout=open(fp_stdout, 'w+'))
    else:
        subprocess.run(args=cmd, cwd=dir_in, stdout=open(devnull, 'w'))


def batch_run_worker(args: List) -> List:
    def worker(
            cmd: str,
            cwd: str,
            fp_stdout: str = None,
            timeout_seconds: int = 1 * 60,
    ) -> List:
        try:
            if fp_stdout:
                subprocess.call(cmd, cwd=cwd, timeout=timeout_seconds, stdout=open(fp_stdout, 'w+'))
            else:
                subprocess.call(cmd, cwd=cwd, timeout=timeout_seconds, stdout=open(devnull, 'w'))
            return [cmd, 'Success']
        except subprocess.TimeoutExpired:
            return [cmd, 'Timed Out']

    kwargs, q = args
    result = worker(**kwargs)
    q.put(1)
    return result


def batch_run(
        list_kwargs_in: List[Dict],
        func_mp: Callable = batch_run_worker,
        n_proc: int = 1,
        dir_work: str = None,
        qt_progress_signal=None
):
    # ------------------------------------------
    # prepare variables used for multiprocessing
    # ------------------------------------------
    import multiprocessing as mp
    m, p = mp.Manager(), mp.Pool(n_proc, maxtasksperchild=1000)
    q = m.Queue()
    jobs = p.map_async(func_mp, [(dict_, q) for dict_ in list_kwargs_in])
    n_simulations = len(list_kwargs_in)

    # ---------------------
    # multiprocessing start
    # ---------------------
    while True:
        if jobs.ready():
            if qt_progress_signal:
                qt_progress_signal.emit(100)
            break  # complete
        else:
            if qt_progress_signal:
                qt_progress_signal.emit(int(q.qsize() / n_simulations * 100))
            time.sleep(1)  # in progress

    # --------------------------------------------
    # pull results and close multiprocess pipeline
    # --------------------------------------------
    p.close()
    p.join()
    mp_out = jobs.get()
    time.sleep(0.5)

    # ----------------------
    # save and print summary
    # ----------------------
    if dir_work:
        out = mp_out
        len_1 = int(max([len(' '.join(i[0])) for i in out]))
        summary = '\n'.join([f'{" ".join(i[0]):<{len_1}} - {i[1]:<{len_1}}' for i in out])
        print(summary)
        with open(path.join(dir_work, 'summary.txt'), 'w+') as f:
            f.write(summary)

    return mp_out
