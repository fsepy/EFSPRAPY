__all__ = (
    'MCS0',
    'main',
    'EXAMPLE_INPUT', 'EXAMPLE_INPUT_DETERMINISTIC'
)

import os
from typing import Callable

import numpy as np
from sfeprapy.mcs import MCSSingle, MCS

from .calcs import main
from .inputs import EXAMPLE_INPUT_DETERMINISTIC, EXAMPLE_INPUT


class MCS0Single(MCSSingle):
    OUTPUT_KEYS = (
        'phi_1', 'phi_2', 'ftp', 't_ig_ftp', 't_ig_safir', 't_max_safir', 'T_max_safir'
    )

    def __init__(self, name, n_simulations, sim_kwargs, save_dir):
        super().__init__(name=name, n_simulations=n_simulations, sim_kwargs=sim_kwargs, save_dir=save_dir)

    @property
    def worker(self) -> Callable:
        return main

    def get_pdf(self, bin_width: float = 0.2) -> (np.ndarray, np.ndarray, np.ndarray):
        ftp: np.ndarray = None
        for i in range(len(self.output_keys)):
            if self.output_keys[i] == 'ftp':
                ftp = self.output[:, i]
        return MCS0Single.make_pdf(ftp, bin_width=bin_width)

    def get_cdf(self, bin_width: float = 0.2):
        x, y_pdf = self.get_pdf(bin_width=bin_width)
        return x, np.cumsum(y_pdf)

    @property
    def output_keys(self) -> tuple:
        return MCS0Single.OUTPUT_KEYS


class MCS0(MCS):
    def __getitem__(self, item) -> MCS0Single:
        return self.mcs_cases[item]

    @property
    def new_mcs_case(self):
        return MCS0Single


def cli_main(fp_mcs_in: str, n_threads: int = 1):
    fp_mcs_in = os.path.realpath(fp_mcs_in)

    mcs = MCS0()
    mcs.set_inputs_file_path(fp_mcs_in)
    mcs.run(n_proc=n_threads)
    mcs.save_all(True)
