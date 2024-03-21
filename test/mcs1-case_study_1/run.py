from copy import deepcopy
from os import path

import numpy as np
from sfeprapy.func.xlsx import dict_to_xlsx
from sfeprapy.mcs import InputParser
from tqdm import tqdm

from efsprapy.mcs1 import MCS1, EXAMPLE_INPUT

if __name__ == '__main__':
    n_simulations = 1000
    base_case_kwargs = deepcopy(EXAMPLE_INPUT['CASE_1'])
    base_case_kwargs.update(dict(
        n_simulations=n_simulations,

        ftp_chf=13.3e3,
        ftp_index=1.7,
        ftp_target=9094,

        fire_fuel_density_MJm2=dict(dist="br187_fuel_load_density_", lbound=10, ubound=5000, mean=None, sd=None, ),
        fire_combustion_efficiency=1,
        fire_hrr_density_kWm2=dict(dist='br187_hrr_density_', lbound=None, ubound=None),
        fire_mode=dict(dist='constant_', lbound=0, ubound=0, values=None, weights=None),
        receiver_ignition_temperature=-1,
        safir_input_file_s=None,
    ))

    kwargs = {'br187_084': deepcopy(base_case_kwargs)}

    for receiver_separation in np.linspace(2, 10, 20):
        kwargs[f'br187_084_{receiver_separation:06.3f}'] = deepcopy(base_case_kwargs)
        kwargs[f'br187_084_{receiver_separation:06.3f}'].update(
            receiver_separation=receiver_separation,
            fire_mode=dict(dist='constant_', lbound=0, ubound=0, values=None, weights=None),
            n_simulations=n_simulations,
        )

    for receiver_separation in np.linspace(2, 10, 20):
        kwargs[f'office_{receiver_separation:06.3f}'] = deepcopy(base_case_kwargs)
        kwargs[f'office_{receiver_separation:06.3f}'].update(
            receiver_separation=receiver_separation,
            fire_fuel_density_MJm2=dict(dist="gumbel_r_", lbound=10, ubound=5000, mean=420, sd=126),
            fire_hrr_density_kWm2=dict(dist="uniform_", lbound=0.15, ubound=0.65),
            fire_mode=dict(dist='constant_', lbound=0, ubound=0, values=None, weights=None),
            n_simulations=n_simulations,
        )

    for receiver_separation in np.linspace(2, 10, 20):
        kwargs[f'office_sprinklered_{receiver_separation:06.3f}'] = deepcopy(base_case_kwargs)
        kwargs[f'office_sprinklered_{receiver_separation:06.3f}'].update(
            receiver_separation=receiver_separation,
            # fire_fuel_density_MJm2=dict(dist="gumbel_r_", lbound=10, ubound=5000, mean=420, sd=126),
            # fire_hrr_density_kWm2=dict(dist="uniform_", lbound=0.15, ubound=0.65),
            fire_mode=dict(dist='discrete_', values='0,1,2', weights='0.07,0.59,0.34', lbound=None, ubound=None),
            n_simulations=n_simulations,
        )

    fp_input = path.join(path.dirname(__file__), "case_study_1.xlsx")
    kwargs_input = {k_: InputParser.flatten_dict(v_) for k_, v_ in kwargs.items()}
    dict_to_xlsx(kwargs_input, fp_input)
    mcs = MCS1()
    mcs.set_inputs_file_path(fp_input)
    pbar = tqdm()
    mcs.run(
        10,
        lambda _: pbar.update(1),
        lambda _: setattr(pbar, "total", _),
        save=True,
        save_archive=False,
        concurrency_strategy=1,
    )
    pbar.close()
