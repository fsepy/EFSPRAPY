import pathlib
import shutil
from os import chdir

import numpy as np

from efsprapy.mcs.parser import InputParser
from efsprapy.mcs.xlsx import dict_to_xlsx
from efsprapy.mcs1 import EXAMPLE_INPUT
from efsprapy.prepare_inputs import convert_input_file_xlsx_to_cases
from efsprapy.run_analysis import process_multiple_cases

if __name__ == '__main__':
    chdir(r'C:\Users\IanFu\Desktop\~1_CURRENT\efsprapy\01-analysis\case_study_solve_upa')
    name = 'sprinklered_cfast'
    dir_project = pathlib.Path(r'C:\Users\IanFu\Desktop\~1_CURRENT\efsprapy\01-analysis\case_study_solve_upa')
    case_name = 'sprinklered_cfast'

    dir_case = dir_project / case_name
    shutil.rmtree(dir_case, ignore_errors=True)
    dir_case.mkdir(parents=True, exist_ok=True)
    fp_xlsx = dir_case / f'{case_name}.xlsx'

    kwargs_ = dict()
    for unprotected_area_percentage in np.linspace(0.1, 1, 10):
        kwargs_[f'{case_name}-{unprotected_area_percentage:06.6f}'] = EXAMPLE_INPUT | dict(
            n_simulations=1000,
            receiver_separation=4.5,
            fire_heat_flux_reduction_factor=unprotected_area_percentage,

            fire_combustion_efficiency=1, receiver_ignition_temperature=-1, safir_input_file_s=None,

            fire_fuel_density_MJm2=dict(dist="gumbel_r_", lbound=10, ubound=3000, mean=780, sd=234),
            fire_hrr_density_kWm2=dict(dist="uniform_", lbound=0.15e3, ubound=0.65e3),
            fire_mode=dict(dist='discrete_', values='0,1,2', weights='0.07,0.59,0.34', lbound=None, ubound=None),
        )

    kwargs_input = {k_: InputParser.flatten_dict(v_) for k_, v_ in kwargs_.items()}
    dict_to_xlsx(kwargs_input, fp_xlsx.as_posix())
    convert_input_file_xlsx_to_cases(fp_xlsx)

    folders = [item for item in dir_case.iterdir() if item.is_dir()]
    folders.sort()
    process_multiple_cases(folders, n_proc=20)
