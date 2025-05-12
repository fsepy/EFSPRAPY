import pathlib
import shutil

import numpy as np

from efsprapy.mcs.parser import InputParser
from efsprapy.mcs.xlsx import dict_to_xlsx
from efsprapy.mcs1 import EXAMPLE_INPUT
from efsprapy.prepare_inputs import convert_input_file_xlsx_to_cases
from efsprapy.run_analysis import process_multiple_cases

if __name__ == '__main__':
    dir_project = pathlib.Path(r'C:\Users\IanFu\Desktop\~1_CURRENT\efsprapy\01-analysis\case_study_solve_separation')
    case_name = 'sprinklered_cfast'

    dir_case = dir_project / case_name
    shutil.rmtree(dir_case, ignore_errors=True)
    dir_case.mkdir(parents=True, exist_ok=True)
    fp_xlsx = dir_case / f'{case_name}.xlsx'

    kwargs_ = dict()

    # for receiver_separation in np.concatenate((
    #         np.arange(0.01, 0.15, 0.04), np.arange(0.15, 0.26 + 1e-5, 0.01), np.arange(0.26, 0.4, 0.02),
    #         np.arange(0.4, 1, 0.2), np.arange(1, 6, 1)
    # )):
    for receiver_separation in np.linspace(0, 15, 16)[1:-1]:
        kwargs_[f'{case_name}-{receiver_separation:06.3f}'] = EXAMPLE_INPUT | dict(
            n_simulations=1000,

            fire_fuel_density_MJm2=dict(dist="gumbel_r_", lbound=10, ubound=3000, mean=780, sd=234),
            fire_hrr_density_kWm2=dict(dist="uniform_", lbound=0.15e3, ubound=0.65e3),
            fire_mode=dict(dist='discrete_', values='0,1,2', weights='0.07,0.59,0.34', lbound=None, ubound=None),
            receiver_separation=receiver_separation,
        )

    dict_to_xlsx({k_: InputParser.flatten_dict(v_) for k_, v_ in kwargs_.items()}, fp_xlsx.as_posix())
    convert_input_file_xlsx_to_cases(fp_xlsx)

    folders = [item for item in dir_case.iterdir() if item.is_dir()]
    folders.sort()
    process_multiple_cases(folders, n_proc=20)
