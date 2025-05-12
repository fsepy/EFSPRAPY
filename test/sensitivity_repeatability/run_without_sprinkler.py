import pathlib
import shutil

from efsprapy.mcs.parser import InputParser
from efsprapy.mcs.xlsx import dict_to_xlsx
from efsprapy.mcs1 import EXAMPLE_INPUT
from efsprapy.prepare_inputs import convert_input_file_xlsx_to_cases
from efsprapy.run_analysis import process_multiple_cases_2

if __name__ == '__main__':
    dir_project = pathlib.Path(r'C:\Users\IanFu\Desktop\~1_CURRENT\efsprapy\01-analysis\sensitivity_repeatability')
    case_name = 'without_sprinkler'

    dir_case = dir_project / case_name
    shutil.rmtree(dir_case, ignore_errors=True)
    dir_case.mkdir(parents=True, exist_ok=True)
    fp_xlsx = dir_case / f'{case_name}.xlsx'

    input_data_dict = dict()
    for n_sim in (100, 200, 500, 1000, 2000, 5000, 10000):
        for i in range(100):
            input_data_dict[f'n{n_sim}_{i}'] = EXAMPLE_INPUT | dict(
                n_simulations=n_sim,
                fire_combustion_efficiency=1, receiver_ignition_temperature=-1, safir_input_file_s=None,
                fire_fuel_density_MJm2=dict(dist="gumbel_r_", lbound=10, ubound=3000, mean=780, sd=234),
                fire_hrr_density_kWm2=dict(dist="uniform_", lbound=0.15e3, ubound=0.65e3),
                fire_mode=dict(dist='constant_', lbound=0, ubound=0, values=None, weights=None),
            )

    dict_to_xlsx({k_: InputParser.flatten_dict(v_) for k_, v_ in input_data_dict.items()}, fp_xlsx.as_posix())
    convert_input_file_xlsx_to_cases(fp_xlsx)

    folders = [item for item in dir_case.iterdir() if item.is_dir()]
    folders.sort()
    process_multiple_cases_2(folders, n_proc=20)
