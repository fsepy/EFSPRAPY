import pathlib
import time
from os import chdir, getcwd

from efsprapy.goal_seek import sep_parallel_any_br187
from efsprapy.mcs.parser import InputParser
from efsprapy.mcs.xlsx import dict_to_xlsx, xlsx_to_dict
from efsprapy.mcs1 import EXAMPLE_INPUT
from efsprapy.prepare_inputs import convert_input_file_xlsx_to_cases
from efsprapy.run_analysis import process_multiple_cases

if __name__ == '__main__':
    chdir(r'C:\Users\IanFu\Desktop\~1_CURRENT\efsprapy\01-analysis\ftp_selection')
    name = 'safir_40mm'

    path_work = pathlib.Path(getcwd()) / name
    path_work.mkdir(parents=True, exist_ok=True)
    fp_xlsx = path_work / f'{name}.xlsx'

    with open(pathlib.Path(__file__).parents[0] / 'hf_ft_40mm.in', 'r') as f:
        therm1d_model = f.read()

    input_data_dict = dict()
    input_data_dict['safir_40mm'] = EXAMPLE_INPUT | dict(
        ftp_chf=13.4e3, ftp_index=2.0, ftp_target=34592,

        opening_width=12,
        opening_height=3,
        room_width=12,
        room_height=3,

        receiver_separation=sep_parallel_any_br187(12, 3, 84, 12.6),
        fire_combustion_efficiency=1, receiver_ignition_temperature=300 + 273.15, safir_input_file_s=therm1d_model,
        fire_fuel_density_MJm2=dict(dist="gumbel_r_", lbound=10, ubound=3000, mean=780, sd=234),
        fire_hrr_density_kWm2=dict(dist="uniform_", lbound=0.15, ubound=0.65),
        fire_mode=dict(dist='constant_', lbound=0, ubound=0, values=None, weights=None),
    )

    input_data_dict_flatten = {k_: InputParser.flatten_dict(v_) for k_, v_ in input_data_dict.items()}
    dict_to_xlsx(input_data_dict_flatten, fp_xlsx.as_posix())
    convert_input_file_xlsx_to_cases(fp_xlsx)
    time.sleep(0.1)

    input_data_dict_from_xlsx = xlsx_to_dict(fp_xlsx.as_posix())
    folders = [item for item in path_work.iterdir() if item.is_dir()]
    folders.sort()
    process_multiple_cases(folders, n_proc=12)
