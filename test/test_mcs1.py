

def test_mcs0_deterministic():
    from efsprapy.mcs1 import hf_ft_40mm_ft
    from efsprapy.mcs1.calcs import main

    res = main(
        t_end=120. * 60.,
        t_step=5.,
        t_end=180. * 60.,
        t_step=10.,

        opening_width=5,
        opening_height=3,
        opening_ventilation_fraction=0.8,

        room_width=12,
        room_height=3,
        room_width_depth_ratio=0.5,

        # fire_mode=dict(dist='discrete_', values='0,1,2', weights='0.07,0.59,0.34'),
        fire_mode=0,
        fire_fuel_density_MJm2=600,
        fire_combustion_efficiency=1,
        fire_hrr_density_kWm2=0.5,
        fire_growth_factor=0.0117,
        fire_t_lim=20. * 60.,

        detector_to_fire_vertical_distance=3 - 0.02,
        detector_to_fire_horizontal_distance=2.83,
        detector_act_temp=93 + 273.15,
        detector_response_time_index=250,

        lining_rho=2000,
        lining_c=1000,
        lining_k=1.13,

        receiver_ignition_temperature=273.15 + 368,
        receiver_separation=7.320449646472931,

        ftp_chf=13.4e3,
        ftp_index=2.0,
        ftp_target=34592,

        safir_input_file_s=hf_ft_40mm_ft,
    )
    res_ = (24261.282744670953, 375.0, 4806344.018119446, 240.0, 2125.0, 1047.76, 0, np.inf, 800000000.0)
    print(f'{res}\n{res_}')
    assert len(res) == len(res_), 'mismatch length'
    for i in range(len(res)):
        assert abs(res[0] - res_[0]) <= 1e-2, f'{res[0]}!={res_[0]}'


def test_mcs1():
    from efsprapy.mcs1 import EXAMPLE_INPUT
    from efsprapy.mcs1 import MCS1
    from tqdm import tqdm

    kwargs = EXAMPLE_INPUT.copy()
    kwargs['CASE_1']['n_simulations'] = 100

    mcs = MCS1()
    mcs.set_inputs_dict(kwargs)
    pbar = tqdm()
    mcs.run(
        10, lambda _: pbar.update(1),
        lambda _: setattr(pbar, 'total', _),
        save=True,
        save_archive=False,
        concurrency_strategy=2
    )
    pbar.close()


def test_radiation_separation_solver():
    from fsetools.lib.fse_thermal_radiation_v2 import phi_solver
    res = phi_solver(W=10, H=2, w=5, h=1, theta=0, Q=84, Q_a=12.6, UA=1)
    res_ = (0.150000048313314, 12.600004058318376, 5.244168619479054, None)
    print(f'{res}\n{res_}')
    assert len(res) == len(res_), 'mismatch length'
    for i in range(len(res)):
        assert abs(res[0] - res_[0]) <= 1e-2, f'{res[0]}!={res_[0]}'


def test_br187_cases_batch():
    from os import path, mkdir

    import numpy as np
    from tqdm import tqdm
    from itertools import product

    from efsprapy.mcs1 import MCS1
    from sfeprapy.mcs import InputParser
    from fsetools.lib.fse_thermal_radiation_v2 import phi_solver

    dir_work = path.join(path.dirname(__file__), 'br187_base_case_analysis')
    try:
        mkdir(dir_work)
    except:
        pass

    Wv, Hv = np.meshgrid(
        np.arange(3., 130 + 3., 3.),
        np.arange(3., 30. + 3., 3.)
    )
    Sv = np.zeros_like(Wv)

    for i, j in product(range(Wv.shape[0]), range(Wv.shape[1])):
        _, _, Sv[i, j], _ = phi_solver(
            W=Wv[i, j], H=Hv[i, j], w=Wv[i, j] / 2, h=Hv[i, j] / 2, theta=0, Q=84, Q_a=12.6, UA=1)

    np.savetxt(
        path.join(dir_work, 'cases_with_s.csv'),
        np.column_stack((Wv.flatten(), Hv.flatten(), Sv.flatten())),
        delimiter=',', header='W,H,S', comments='', fmt='%.3f'
    )

    def make_mcs_kwargs(W: float, H: float, S: float):
        return dict(CASE_1=dict(
            n_simulations=100000,
            t_end=180. * 60.,
            t_step=10.,

            opening_width=W,
            opening_height=H,
            opening_ventilation_fraction=dict(dist='lognorm_mod_', lbound=1e-4, ubound=1 - 1e-4, mean=0.2, sd=0.2),

            room_width=W,
            room_width_depth_ratio=dict(dist='uniform_', lbound=0.4, ubound=0.6),
            room_height=H,

            # fire_mode=dict(dist='discrete_', values='0,1,2', weights='0.07,0.59,0.34'),
            fire_mode=0,
            fire_fuel_density_MJm2=dict(dist="br187_fuel_load_density_", lbound=10, ubound=5000),
            fire_combustion_efficiency=dict(dist='uniform_', lbound=0.8, ubound=1.0),
            fire_hrr_density_kWm2=dict(dist='br187_hrr_density_'),
            fire_growth_factor=0.0117,
            fire_t_lim=20. * 60.,
            fire_convection_factor=0.7,

            detector_to_fire_vertical_distance=H - 0.02,
            detector_to_fire_horizontal_distance=2.83,
            detector_act_temp=93 + 273.15,
            detector_response_time_index=250,
            detector_conduction_factor=0.65,

            lining_rho=2000,
            lining_c=1000,
            lining_k=1.13,

            receiver_ignition_temperature=-(300 + 273.15),
            receiver_emissivity=1,
            receiver_separation=S,

            ftp_chf=13.4e3,
            ftp_index=2.0,
            ftp_target=34592,
        ))

    kwargs = dict()
    Iv = np.zeros_like(Wv, dtype=int)
    for i in range(Wv.shape[0]):
        for j in range(Wv.shape[1]):
            W, H, S = Wv[i, j], Hv[i, j], Sv[i, j]
            case_name = f'{W:03g}-{H:03g}-{S:06.3f}'
            kwargs[case_name] = make_mcs_kwargs(W=W, H=H, S=S)['CASE_1'].copy()

        fp_input = path.join(dir_work, f'br187_base_cases_w{i:03f}.xlsx')
        dict_to_xlsx({i: InputParser.flatten_dict(v) for i, v in kwargs.items()}, fp_input)

        pbar = tqdm(total=sum([v['n_simulations'] for k, v in kwargs.items()]), )
        mcs = MCS1()
        mcs.set_inputs_file_path(fp_input)
        mcs.run(n_proc=2, set_progress=lambda x: pbar.update(1), save=True, save_archive=True)
        pbar.close()

        for j in range(Wv.shape[1]):
            W, H, S = Wv[i, j], Hv[i, j], Sv[i, j]
            case_name = f'{W:03g}-{H:03g}-{S:06.3f}'
            t_ig_ftp = mcs[case_name].output[:, 1]
            Iv[i, j] = int(sum(np.bitwise_and(t_ig_ftp != np.inf, t_ig_ftp != np.nan)))

    _ = np.column_stack((Wv.flatten(), Hv.flatten(), Sv.flatten(), Iv.flatten()))
    np.savetxt(
        path.join(dir_work, 's_and_i.csv'), _, delimiter=',', header='W,H,S,I', comments='', fmt='%.3f'
    )


def test_ig_probability_and_ftp_calibration():
    from efsprapy.mcs1 import EXAMPLE_INPUT
    EXAMPLE_INPUT_ = EXAMPLE_INPUT['CASE_1'].copy()
    EXAMPLE_INPUT_['n_simulations'] = 10_000
    EXAMPLE_INPUT_.pop('ftp_chf')
    EXAMPLE_INPUT_.pop('ftp_index')
    EXAMPLE_INPUT_.pop('ftp_target')

    kwargs = {
        "40 mm macrocarpa": dict(
            ftp_chf=13.4e3,
            ftp_index=2.0,
            ftp_target=34592,
        ),
        "40 mm beech": dict(
            ftp_chf=10.9e3,
            ftp_index=2.0,
            ftp_target=42750,
        ),
        "40 mm mdf": dict(
            ftp_chf=9.3e3,
            ftp_index=1.5,
            ftp_target=10599,
        ),
        "40 mm radiata pine": dict(
            ftp_chf=8.6351e3,
            ftp_index=2.0,
            ftp_target=41209,
        ),
        "40 mm rimu": dict(
            ftp_chf=8.3e3,
            ftp_index=2.0,
            ftp_target=61474,
        ),
        "40 mm plywood": dict(
            ftp_chf=7.3e3,
            ftp_index=1.9,
            ftp_target=39663,
        )
    }

    for v in kwargs.values():
        v.update(**EXAMPLE_INPUT_)

    from os import path
    dir_work = ''
    fp_input = path.join(dir_work, f'test-ftp_calibration.xlsx')
    from sfeprapy.mcs import InputParser
    from sfeprapy.func.xlsx import dict_to_xlsx
    dict_to_xlsx(
        {k_: InputParser.flatten_dict(v_) for k_, v_ in kwargs.items()},
        fp_input
    )
    from efsprapy.mcs1 import MCS1
    mcs = MCS1()
    mcs.set_inputs_file_path(fp_input)
    from tqdm import tqdm
    pbar = tqdm(total=sum([i.n_sim for i in mcs.mcs_cases.values()]))
    mcs.run(n_proc=4, set_progress=lambda x: pbar.update(1), save=True, save_archive=False)
    pbar.close()
    for mcs_case in mcs.mcs_cases.values():
        t_ig = mcs_case.output[:, 1]
        print(f'{mcs_case.name:<20.20}: {sum(np.logical_and(t_ig > 0, t_ig < np.inf)):<10}')


if __name__ == "__main__":
    # test_calculate_incident_heat_flux_from_controlled_fire_cfast()
    # test_mcs0_deterministic()
    test_mcs1()
