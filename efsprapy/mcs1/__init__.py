__all__ = (
    'main',

    'EXAMPLE_INPUT', 'EXAMPLE_INPUT_DETERMINISTIC',

    'hf_ft_40mm_ft', 'hf_ft_20mm_ft'
)

from os import path

from .calcs import main
from .inputs import EXAMPLE_INPUT_DETERMINISTIC, EXAMPLE_INPUT
from .safir_input_files import hf_ft_40mm_ft, hf_ft_20mm_ft


def cli_main(fp_mcs_in: str, n_threads: int = 1):
    fp_mcs_in = path.realpath(fp_mcs_in)

    mcs = MCS1()
    mcs.set_inputs_file_path(fp_mcs_in)
    mcs.run(n_proc=n_threads)
    mcs.save_all(True)
