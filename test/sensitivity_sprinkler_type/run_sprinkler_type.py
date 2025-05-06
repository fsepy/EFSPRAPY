import pathlib
from os import getcwd, chdir

import numpy as np

from efsprapy.mcs.parser import InputParser
from efsprapy.mcs.xlsx import dict_to_xlsx, xlsx_to_dict
from efsprapy.mcs1 import EXAMPLE_INPUT
from efsprapy.prepare_inputs import convert_input_file_xlsx_to_cases
from efsprapy.run_analysis import process_multiple_cases

if __name__ == '__main__':
    chdir()
    name = 'sprinkler_type'

    path_work = pathlib.Path(getcwd()) / name
    path_work.mkdir(parents=True, exist_ok=True)
    fp_xlsx = path_work / f'{name}.xlsx'

    sprinkler_data = dict(
        residential=dict(
            detector_act_temp=68 + 273.15,
            detector_response_time_index=50,
            detector_conduction_factor=0.02,
        ),
        office=dict(
            detector_act_temp=93 + 273.15,
            detector_response_time_index=150,
            detector_conduction_factor=0.04,
        ),
        retail=dict(
            detector_act_temp=141 + 273.15,
            detector_response_time_index=200,
            detector_conduction_factor=0.05,
        ),
        basecase=dict(
            detector_act_temp=93 + 273.15,
            detector_response_time_index=250,
            detector_conduction_factor=0.65,
        )
    )

    input_data_dict = dict()
    for k, v in sprinkler_data.items():
        for sep_dist in np.linspace(4, 10, 15):
            input_data_dict[f'{k}-{sep_dist:.5f}'] = EXAMPLE_INPUT['CASE_1'] | dict(
                n_simulations=5000,

                fire_fuel_density_MJm2=dict(dist="gumbel_r_", lbound=10, ubound=3000, mean=780, sd=234),
                fire_hrr_density_kWm2=dict(dist="uniform_", lbound=0.15, ubound=0.65),

                fire_combustion_efficiency=1, receiver_ignition_temperature=-1, safir_input_file_s=None,
                fire_mode=dict(dist='discrete_', values='0,1,2', weights='0.07,0.59,0.34', lbound=None, ubound=None),
                receiver_separation=sep_dist,

                **v,
            )

    input_data_dict_flatten = {k_: InputParser.flatten_dict(v_) for k_, v_ in input_data_dict.items()}
    dict_to_xlsx(input_data_dict_flatten, fp_xlsx.as_posix())
    convert_input_file_xlsx_to_cases(fp_xlsx)

    input_data_dict_from_xlsx = xlsx_to_dict(fp_xlsx.as_posix())
    folders = [item for item in path_work.iterdir() if item.is_dir()]
    folders.sort()
    process_multiple_cases(folders, n_proc=24)
