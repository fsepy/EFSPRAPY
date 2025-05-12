import itertools
import pathlib
import shutil

import numpy as np

from efsprapy.goal_seek import sep_parallel_any_br187
from efsprapy.mcs1 import EXAMPLE_INPUT
from efsprapy.prepare_inputs import convert_input_project_data_to_cases
from efsprapy.run_analysis import process_multiple_cases_2

if __name__ == '__main__':
    dir_project = pathlib.Path(r'C:\Users\IanFu\Desktop\~1_CURRENT\efsprapy\01-analysis\failure_probability_contour')
    case_name = 'contour_with_2m_limit'

    dir_case = dir_project / case_name
    shutil.rmtree(dir_case, ignore_errors=True)
    dir_case.mkdir(parents=True, exist_ok=True)
    fp_xlsx = dir_case / f'{case_name}.xlsx'

    kwargs_ = dict()
    for w, h in itertools.product((np.arange(0.5, 21.001, 0.5)), np.arange(6, 9.001, 0.5)):
        kwargs_[f'{fp_xlsx.stem}-{w:06.3f}-{h:06.3f}'] = EXAMPLE_INPUT | dict(
            receiver_separation=max(2, sep_parallel_any_br187(w, h, 84, 12.6)),
            room_width=w, room_height=h, opening_width=w, opening_height=h,

            fire_combustion_efficiency=1, receiver_ignition_temperature=-1, safir_input_file_s=None,

            fire_fuel_density_MJm2=dict(dist="gumbel_r_", lbound=10, ubound=3000, mean=780, sd=234),
            fire_hrr_density_kWm2=dict(dist="uniform_", lbound=0.15e3, ubound=0.65e3),
            fire_mode=0,
        )

    # dict_to_xlsx({k_: InputParser.flatten_dict(v_) for k_, v_ in kwargs_.items()}, fp_xlsx.as_posix())
    convert_input_project_data_to_cases(dir_case, kwargs_, n_proc=30)

    folders = [item for item in dir_case.iterdir() if item.is_dir()]
    folders.sort()
    process_multiple_cases_2(folders, n_proc=20)
