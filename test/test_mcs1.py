def test_calculate_incident_heat_flux_from_controlled_fire_cfast():
    import numpy as np
    t = np.arange(0., 7200. + 10. / 2., 10)
    from efsprapy.mcs1.calcs import calculate_incident_heat_flux_from_controlled_fire_cfast
    res = calculate_incident_heat_flux_from_controlled_fire_cfast(
        fire_mode=0, t_arr=t, W=5, D=3, H=2, W_v=5, H_v=2, q_x_d=800e6, hrr_density_kWm2=250, alpha=0.0117, H_d=2,
        R_d=3, RTI=250, C=0.5, T_act=68 + 273.15, S=2
    )
    print([i if isinstance(i, (int, float)) else np.average(i) for i in res])


def test_mcs0_deterministic():
    from efsprapy.mcs1.calcs import main
    res = main(
        t_end=180. * 60.,
        t_step=5.,

        opening_width=12,
        opening_height=3,

        room_height=3,
        room_width=5,
        room_depth=96 / 5,

        fire_mode=1,
        fire_fuel_density_MJm2=800,  # 420
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

        receiver_ignition_temperature=-(300 + 273.15),
        receiver_separation=2,

        ftp_chf=12.6e3,
        ftp_index=2,
        ftp_target=34592,
    )
    # print(phi_1, phi_2, ftp, t_ig_ftp, t_ig_safir, t_max_safir, T_max_safir, fire_mode)
    # assert abs(phi_1 - 0.291670) <= 1e-2, f'{phi_1}!=0.291670'
    # assert abs(phi_2 - 0.443859) <= 1e-2, f'{phi_2}!=0.443859'
    # assert abs(ftp - 56531.682) <= 1e-2, f'{ftp}!=56531.682'
    # assert abs(t_ig_ftp - 2235.) <= 1e-2, f'{t_ig_ftp}!=2250.'
    # assert abs(t_ig_safir - 450.) <= 1e-2, f'{t_ig_safir}!=450.'
    print(res)


def test_mcs0_research_v2():
    input_kwargs = dict(
        CASE_1=dict(
            n_simulations=100,
            t_end=180. * 60.,
            t_step=5.,

            opening_width=12,
            opening_height=3,

            room_height=3,
            room_width=5,
            room_depth=96 / 5,

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

    from efsprapy.mcs1 import MCS1
    from tqdm import tqdm

    pbar = tqdm(total=sum([v['n_simulations'] for k, v in input_kwargs.items()]))
    mcs = MCS1()
    mcs.set_inputs_dict(input_kwargs)
    mcs.run(n_proc=2, set_progress=lambda x: pbar.update(1), save=True, save_archive=False)


if __name__ == "__main__":
    # test_calculate_incident_heat_flux_from_controlled_fire_cfast()
    test_mcs0_research_v2()
    # test_mcs0_deterministic()
