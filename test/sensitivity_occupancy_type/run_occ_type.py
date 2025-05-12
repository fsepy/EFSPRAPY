import pathlib
import shutil

import numpy

from efsprapy.mcs.parser import InputParser
from efsprapy.mcs.xlsx import dict_to_xlsx
from efsprapy.mcs1 import EXAMPLE_INPUT
from efsprapy.prepare_inputs import convert_input_file_xlsx_to_cases
from efsprapy.run_analysis import process_multiple_cases_2

if __name__ == '__main__':
    dir_project = pathlib.Path(r'C:\Users\IanFu\Desktop\~1_CURRENT\efsprapy\01-analysis\sensitivity_occupancy_type')
    case_name = 'occ_type'

    dir_case = dir_project / case_name
    shutil.rmtree(dir_case, ignore_errors=True)
    dir_case.mkdir(parents=True, exist_ok=True)
    fp_xlsx = dir_case / f'{case_name}.xlsx'

    fuel_load_data = dict(
        dwelling=dict(dist="gumbel_r_", lbound=10, ubound=2775, mean=780, sd=234),
        office=dict(dist="gumbel_r_", lbound=10, ubound=1500, mean=420, sd=126),
        retail=dict(dist="gumbel_r_", lbound=10, ubound=2134, mean=600, sd=180),
        storage=dict(dist="lognorm_", lbound=10, ubound=3641, mean=615, sd=355),
        cinema=dict(dist="gumbel_r_", lbound=10, ubound=940, mean=300, sd=75),
        plant_room=dict(dist="gumbel_r_", lbound=10, ubound=1216, mean=235, sd=115),
        cycle_store=dict(dist="gumbel_r_", lbound=10, ubound=1500, mean=190, sd=105),
        school=dict(dist="gumbel_r_", lbound=10, ubound=1500, mean=285, sd=85.5),
        kitchen=dict(dist="lognorm_", lbound=10, ubound=1500, mean=314, sd=160.14),
        restaurant=dict(dist="lognorm_", lbound=10, ubound=1112, mean=298, sd=190.72),
        car_park=dict(dist="norm_", lbound=222, ubound=282, mean=252, sd=7.02),
        library=dict(dist="gumbel_r_", lbound=10, ubound=3000, mean=1500, sd=450),
    )
    hrrpua_data = dict(
        dwelling=dict(dist="uniform_", lbound=0.32e3, ubound=0.57e3),
        office=dict(dist="uniform_", lbound=0.15e3, ubound=0.65e3),
        retail=dict(dist="uniform_", lbound=0.27e3, ubound=1e3),
        storage=dict(dist="uniform_", lbound=0.27e3, ubound=1e3),
        cinema=dict(dist="uniform_", lbound=0.499e3, ubound=0.501e3),
        plant_room=dict(dist="uniform_", lbound=0.09e3, ubound=0.62e3),
        cycle_store=dict(dist="uniform_", lbound=0.4e3, ubound=1e3),
        school=dict(dist="uniform_", lbound=0.249e3, ubound=0.251e3),
        kitchen=dict(dist="uniform_", lbound=0.27e3, ubound=1e3),
        restaurant=dict(dist="uniform_", lbound=0.15e3, ubound=0.65e3),
        car_park=dict(dist="uniform_", lbound=0.09e3, ubound=0.62e3),
        library=dict(dist="uniform_", lbound=0.4e3, ubound=2e3),
    )
    growth_factor_data = dict(
        dwelling=0.012,
        office=0.012,
        retail=0.047,
        storage=0.047,
        cinema=0.012,
        plant_room=0.012,
        cycle_store=0.012,
        school=0.012,
        kitchen=0.047,
        restaurant=0.047,
        car_park=0.035,
        library=0.047,
    )

    input_data_dict = dict()
    for k, v in fuel_load_data.items():
        for sep_dist in numpy.arange(2, 15.001, 0.5):
            input_data_dict[f'{k}-{sep_dist:.5f}'] = EXAMPLE_INPUT | dict(
                fire_fuel_density_MJm2=v,
                fire_hrr_density_kWm2=hrrpua_data[k],
                fire_growth_factor=growth_factor_data[k],

                n_simulations=10_000,
                fire_combustion_efficiency=1, receiver_ignition_temperature=-1, safir_input_file_s=None,
                fire_mode=dict(dist='constant_', lbound=0, ubound=0, values=None, weights=None),
                receiver_separation=sep_dist,
            )

    dict_to_xlsx({k_: InputParser.flatten_dict(v_) for k_, v_ in input_data_dict.items()}, fp_xlsx.as_posix())
    convert_input_file_xlsx_to_cases(fp_xlsx, 20)

    folders = [item for item in dir_case.iterdir() if item.is_dir()]
    folders.sort()
    process_multiple_cases_2(folders, n_proc=20)
