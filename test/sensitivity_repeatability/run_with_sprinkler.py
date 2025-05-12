import pathlib
from os import getcwd, chdir

from efsprapy.mcs.parser import InputParser
from efsprapy.mcs.xlsx import dict_to_xlsx, xlsx_to_dict
from efsprapy.mcs1 import EXAMPLE_INPUT
from efsprapy.prepare_inputs import convert_input_file_xlsx_to_cases
from efsprapy.run_analysis import process_multiple_cases

if __name__ == '__main__':
    chdir(r'C:\Users\IanFu\Desktop\~1_CURRENT\efsprapy\01-analysis\sensitivity_repeatability')
    name = 'with_sprinkler'

    path_work = pathlib.Path(getcwd()) / name
    path_work.mkdir(parents=True, exist_ok=True)
    fp_xlsx = path_work / f'{name}.xlsx'

    input_data_dict = dict()
    for n_sim in (100, 200, 500, 1000, 2000, 5000,):
        for i in range(100):
            input_data_dict[f'n{n_sim}_{i}'] = EXAMPLE_INPUT | dict(
                n_simulations=n_sim,
                fire_combustion_efficiency=1, receiver_ignition_temperature=-1, safir_input_file_s=None,
                fire_fuel_density_MJm2=dict(dist="gumbel_r_", lbound=10, ubound=3000, mean=780, sd=234),
                fire_hrr_density_kWm2=dict(dist="uniform_", lbound=0.15e3, ubound=0.65e3),
                fire_mode=dict(dist='discrete_', values='0,1,2', weights='0.07,0.59,0.34', lbound=None, ubound=None),
            )

    input_data_dict_flatten = {k_: InputParser.flatten_dict(v_) for k_, v_ in input_data_dict.items()}
    dict_to_xlsx(input_data_dict_flatten, fp_xlsx.as_posix())
    convert_input_file_xlsx_to_cases(fp_xlsx)

    input_data_dict_from_xlsx = xlsx_to_dict(fp_xlsx.as_posix())
    folders = [item for item in path_work.iterdir() if item.is_dir()]
    folders.sort()
    process_multiple_cases(folders, n_proc=23)
