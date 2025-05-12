import pathlib
import shutil

from efsprapy.mcs.parser import InputParser
from efsprapy.mcs.xlsx import dict_to_xlsx
from efsprapy.mcs1 import EXAMPLE_INPUT
from efsprapy.prepare_inputs import convert_input_project_data_to_cases
from efsprapy.run_analysis import process_multiple_cases

if __name__ == '__main__':
    dir_project = pathlib.Path(r'C:\Users\IanFu\Desktop\~1_CURRENT\efsprapy\01-analysis\sensitivity_sprinkler_type')
    case_name = 'rti'

    dir_case = dir_project / case_name
    shutil.rmtree(dir_case, ignore_errors=True)
    dir_case.mkdir(parents=True, exist_ok=True)
    fp_xlsx = dir_case / f'{case_name}.xlsx'

    input_data_dict = dict()
    for rti in (30, 50, 80, 100, 150, 200, 250, 300, 350):
        input_data_dict[f'{rti:05d}'] = EXAMPLE_INPUT | dict(
            n_simulations=100,

            opening_width=10,
            opening_height=2.5,
            room_width=10,
            room_height=2.5,
            detector_to_fire_vertical_distance=2.5 - 0.02,

            fire_fuel_density_MJm2=dict(dist="gumbel_r_", lbound=10, ubound=3000, mean=780, sd=234),
            fire_hrr_density_kWm2=dict(dist="uniform_", lbound=0.15e3, ubound=0.65e3),

            fire_combustion_efficiency=1, receiver_ignition_temperature=-1, safir_input_file_s=None,
            fire_mode=dict(dist='discrete_', values='0,1,2', weights='0.07,0.59,0.34', lbound=None, ubound=None),
            receiver_separation=.4,

            detector_act_temp=68 + 273.15,
            detector_response_time_index=rti,
            detector_conduction_factor=0.65,
        )

    dict_to_xlsx({k_: InputParser.flatten_dict(v_) for k_, v_ in input_data_dict.items()}, fp_xlsx.as_posix())
    convert_input_project_data_to_cases(dir_case, input_data_dict, n_proc=30)

    folders = [item for item in dir_case.iterdir() if item.is_dir()]
    folders.sort()
    process_multiple_cases(folders, n_proc=20)
