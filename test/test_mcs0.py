def test_mcs0_deterministic():
    from efsprapy.mcs0.calcs import main
    from efsprapy.mcs0.inputs import EXAMPLE_INPUT_DETERMINISTIC
    *_, ftp = main(**EXAMPLE_INPUT_DETERMINISTIC)
    assert abs(ftp - 63004611.94969697) <= 1e-5


def test_mcs0():
    from efsprapy.mcs0 import MCS0, EXAMPLE_INPUT
    mcs = MCS0()
    mcs.set_inputs_dict(EXAMPLE_INPUT.copy())
    mcs.run(n_proc=4, save=True, save_archive=False)


def test_mcs0_research():
    from efsprapy.mcs0 import MCS0, EXAMPLE_INPUT

    custom_input = EXAMPLE_INPUT.copy()
    custom_input['CASE_1']['chf'] = 13_400
    custom_input['CASE_1']['ftp_index'] = 2
    custom_input['CASE_1']['ftp_target'] = 34_592
    custom_input['CASE_1']['emissivity'] = dict(dist='uniform_', lbound=0.5, ubound=1.0)
    mcs = MCS0()
    mcs.set_inputs_dict(custom_input)
    mcs.run(n_proc=1, save=True, save_archive=False)


def test_parametric_fire():
    import numpy as np
    from efsprapy.mcs0.calcs import IgnitionTimeSolverIncidentHeatFluxSafir

    from fsetools.lib.fse_bs_en_1991_1_2_parametric_fire import temperature as param_temperature

    _t_ = np.arange(0, 10800 + 0.5, 10)
    _T_ = param_temperature(
        A_t=360, A_f=100, A_v=36.1, h_eq=1., q_fd=600e6, lbd=1., rho=1., c=2250000, t_lim=20 * 60, t=_t_, T_0=293.15
    )
    _q_ = 5.67e-8 * 1.0 * _T_ ** 4

    run = IgnitionTimeSolverIncidentHeatFluxSafir(t=_t_, q_1=_q_, phi_1=0.1, epsilon_1=1.0)
    t_ig, n_iter = run.solve(solver_t_ig_tol=10, T_ig=273.15 + 300)
    print(t_ig, n_iter)
    assert abs(t_ig - 760) < 10.


def test_controlled_fire():
    import numpy as np
    from efsprapy.func.controlled_fire import heat_flux_from_controlled_fire
    from efsprapy.mcs0.calcs import IgnitionTimeSolverIncidentHeatFluxSafir

    _t_ = np.arange(0, 10800 + 0.5, 10)
    q_f, q_s = heat_flux_from_controlled_fire(
        t=np.arange(0, 10800 + 0.5, 10),  # Fire duration [mins]
        fire_hrr_density_kWm2=1000,  # HRR per unit area [kw/m²]
        fire_alpha=0.047,  # Fire growth parameter [kW/s²]
        fire_conv_frac=0.7,  # Convective fraction
        H=4.73,  # Height of ceiling [m]
        S=3,
        W_o=4,
        detector_to_fire_vertical_distance=4.73,
        detector_act_temp=93 + 273.15,
        detector_to_fire_horizontal_distance=5.66,  # Radial distance from the ceiling impingement point [m]
        detector_response_time_index=135,  # Sprinkler response time index [(m.s)^1/2]
        detector_conduction_factor=0.65,  # Conduction factor [(m.s)^1/2]
    )

    _t_ = np.arange(0, 10800 + 0.5, 10)

    run = IgnitionTimeSolverIncidentHeatFluxSafir(
        t=_t_,
        q_1=q_f, phi_1=0.5, epsilon_1=1.0,
        q_2=q_s, phi_2=0.8, epsilon_2=1.0,
    )
    t_ig, n_iter = run.solve(solver_t_ig_tol=5, T_ig=273.15 + 200)
    print(t_ig, n_iter)


if __name__ == "__main__":
    pass
