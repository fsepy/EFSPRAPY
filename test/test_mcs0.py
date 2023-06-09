def test_mcs0_deterministic():
    from efsprapy.mcs0.calcs import main
    from efsprapy.mcs0.inputs import EXAMPLE_INPUT_DETERMINISTIC
    phi_1, phi_2, ftp, t_ig_ftp, t_ig_safir, t_max_safir, T_max_safir, fire_mode, t_d = main(
        **EXAMPLE_INPUT_DETERMINISTIC
    )
    print(phi_1, phi_2, ftp, t_ig_ftp, t_ig_safir, t_max_safir, T_max_safir, fire_mode)
    assert abs(phi_1 - 0.291670) <= 1e-2, f'{phi_1}!=0.291670'
    assert abs(phi_2 - 0.443859) <= 1e-2, f'{phi_2}!=0.443859'
    assert abs(ftp - 56531.682) <= 1e-2, f'{ftp}!=56531.682'
    assert abs(t_ig_ftp - 2235.) <= 1e-2, f'{t_ig_ftp}!=2250.'
    assert abs(t_ig_safir - 450.) <= 1e-2, f'{t_ig_safir}!=450.'


def test_mcs0():
    from efsprapy.mcs0 import MCS0, EXAMPLE_INPUT
    mcs = MCS0()
    mcs.set_inputs_dict(EXAMPLE_INPUT.copy())
    mcs.run(n_proc=4, save=True, save_archive=False)


def test_mcs0_research():
    from efsprapy.mcs0 import MCS0, EXAMPLE_INPUT

    custom_input = EXAMPLE_INPUT.copy()
    custom_input['CASE_1']['n_simulations'] = 10
    custom_input['CASE_1']['chf'] = 13_400
    custom_input['CASE_1']['ftp_index'] = 2
    custom_input['CASE_1']['ftp_target'] = 34_592
    custom_input['CASE_1']['emissivity'] = dict(dist='uniform_', lbound=0.5, ubound=1.0)
    mcs = MCS0()
    mcs.set_inputs_dict(custom_input)
    mcs.run(n_proc=1, save=True, save_archive=False)


def test_fire_mode_0():
    import numpy as np

    from fsetools.lib.fse_bs_en_1991_1_2_parametric_fire import temperature as param_temperature
    from efsprapy.mcs0.calcs import calculate_ignition_time_temperature

    _t_ = np.arange(0, 10800 + 0.5, 10)
    _T_ = param_temperature(
        A_t=360, A_f=100, A_v=36.1, h_eq=1., q_fd=600e6, lbd=1., rho=1., c=2250000, t_lim=20 * 60, t=_t_, T_0=293.15
    )
    _q_ = 5.67e-8 * 1.0 * _T_ ** 4

    t_ig, *_ = calculate_ignition_time_temperature(t=_t_, q_inc=_q_, T_ig=273.15 + 150)

    assert abs(t_ig - 170.) < 5., f'{t_ig}!=170.'


def test_fire_mode_1():
    import numpy as np
    from efsprapy.mcs0.calcs import calculate_ignition_time_temperature
    from efsprapy.mcs0.calcs import calculate_incident_heat_flux_from_controlled_fire

    _t_ = np.arange(0, 10800 + 0.5, 1)
    q_f, phi_f, q_s, phi_s, t_d = calculate_incident_heat_flux_from_controlled_fire(
        fire_mode=1,
        t_arr=_t_,
        hrr_density_kWm2=510,
        alpha=0.012,
        H_d=2.4,
        R_d=2.75,
        RTI=115,
        C=0.4,
        T_act=150 + 273.15,
        H=2.75,
        S=4,
        W_v=36.1,
        H_v=2.,
        A_t=360,
        A_f=100,
        t_alpha=300,
        b=1500,
        q_x_d=700,
    )

    t_ig, t_max, T_max = calculate_ignition_time_temperature(t=_t_, T_ig=273.15 + 50, q_inc=q_f + q_s, )

    print(t_ig, t_max, T_max)
    assert abs(t_ig - 202.) <= 1, f'{t_ig}!=202.'
    assert abs(t_max - 10800.) <= 1, f'{t_max}!=10800.'
    assert abs(T_max - 1005.54) <= 1, f'{T_max}!=1005.54'


def test_fire_mode_2():
    import numpy as np
    from efsprapy.mcs0.calcs import calculate_ignition_time_temperature
    from efsprapy.mcs0.calcs import calculate_incident_heat_flux_from_controlled_fire

    _t_ = np.arange(0, 10800 + 0.5, 1)
    q_f, phi_f, q_s, phi_s = calculate_incident_heat_flux_from_controlled_fire(
        fire_mode=2,
        t_arr=_t_,
        hrr_density_kWm2=510,
        alpha=0.012,
        H_d=2.4,
        R_d=2.75,
        RTI=115,
        C=0.4,
        T_act=150 + 273.15,
        H=2.75,
        S=4,
        W_v=6. / 2.5,
        H_v=2.5,
        A_t=80,
        A_f=16,
        t_alpha=300,
        b=1500,
        q_x_d=511e6,
    )

    t_ig, t_max, T_max = calculate_ignition_time_temperature(t=_t_, T_ig=273.15 + 50, q_inc=q_f + q_s, )

    print(t_ig, t_max, T_max)
    assert abs(t_ig - 202.) <= 1, f'{t_ig}!=202.'
    assert abs(t_max - 510.) <= 1, f'{t_max}!=510.'
    assert abs(T_max - 998.43) <= 1, f'{T_max}!=998.43'


def test_mcs0_research_v2():
    EXAMPLE_INPUT = dict(
        CASE_1=dict(
            n_simulations=100,
            t_end=180. * 60.,
            t_step=5.,

            opening_width=12,
            opening_height=3,

            room_height=3,
            room_floor_area=96,
            room_total_surface_area=312,

            fire_mode=dict(dist='discrete_', values='0,1,2', weights='0.07,0.59,0.34'),
            fire_fuel_density_MJm2=dict(dist="gumbel_r_", lbound=10, ubound=1200, mean=420, sd=126),  # 420
            fire_hrr_density_kWm2=250,
            fire_growth_factor=0.0117,
            fire_t_lim=20. * 60.,
            fire_din_growth_factor=300,  # todo
            fire_convection_factor=0.7,

            detector_to_fire_vertical_distance=3,
            detector_to_fire_horizontal_distance=5.6,  # todo, use sprinkler effectiveness?
            detector_act_temp=93 + 273.15,
            detector_response_time_index=135,
            detector_conduction_factor=0.65,

            lining_rho=2000,
            lining_c=1000,
            lining_k=1.13,
            lining_thermal_effusivity=1500,

            receiver_ignition_temperature=300 + 273.15,
            receiver_emissivity=1,
            receiver_separation=2,

            ftp_chf=12.6e3,
            ftp_index=dict(dist='uniform_', lbound=1, ubound=2),
            ftp_target=34592,
        )
    )

    from efsprapy.mcs0 import MCS0

    mcs = MCS0()
    mcs.set_inputs_dict(EXAMPLE_INPUT.copy())
    mcs.run(n_proc=2, save=True, save_archive=False)


if __name__ == "__main__":
    test_mcs0_deterministic()
    test_mcs0()
    test_fire_mode_0()
    test_fire_mode_1()
    test_fire_mode_2()
    test_mcs0_research_v2()
