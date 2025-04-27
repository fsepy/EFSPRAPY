import pathlib

import numpy as np

from efsprapy.mcs1 import EXAMPLE_INPUT
from efsprapy.prepare_inputs import prepare_inputs_with_var_sep_dist
from efsprapy.run_analysis import process_multiple_cases

if __name__ == '__main__':
    name = 'cfast_sprinkler'
    path_work = pathlib.Path(__file__).parents[0] / name
    path_work.mkdir(parents=True, exist_ok=True)

    prepare_inputs_with_var_sep_dist(
        n_simulations=10_000,
        kwargs=EXAMPLE_INPUT['CASE_1'] | dict(
            ftp_chf=13.3e3, ftp_index=1.7, ftp_target=9094,

            fire_combustion_efficiency=1, receiver_ignition_temperature=-1, safir_input_file_s=None,

            fire_fuel_density_MJm2=dict(dist="gumbel_r_", lbound=123.24, ubound=1494.3, mean=420, sd=126),
            fire_hrr_density_kWm2=dict(dist="uniform_", lbound=0.15, ubound=0.65),
            fire_mode=dict(dist='discrete_', values='0,1,2', weights='0.07,0.59,0.34', lbound=None, ubound=None),
        ),
        receiver_separations=np.linspace(2, 10, 100).tolist(),
        fp_xlsx=path_work / f'{name}.xlsx',
    )

    folders = [item for item in path_work.iterdir() if item.is_dir()]
    folders.sort()
    process_multiple_cases(folders)
