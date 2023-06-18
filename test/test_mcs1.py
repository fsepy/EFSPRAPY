

def test_mcs0_deterministic():
    from efsprapy.mcs1 import hf_ft_40mm_ft
    from efsprapy.mcs1.calcs import main

    res = main(
        t_end=180. * 60.,
        t_step=10.,

        opening_width=12,
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
    res_ = [510.4296006211892, 1.0, 1.1550856590671998, 0.3899564855008652, 0]
    print(res, '\n', res_)
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


if __name__ == "__main__":
    # test_calculate_incident_heat_flux_from_controlled_fire_cfast()
    # test_mcs0_deterministic()
    test_mcs1()
