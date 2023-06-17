EXAMPLE_INPUT_DETERMINISTIC = dict(
    t_end=180. * 60.,
    t_step=10.,

    opening_width=9,
    opening_height=3,
    opening_ventilation_fraction=dict(dist='lognorm_mod_', lbound=1e-4, ubound=1 - 1e-4, mean=0.2, sd=0.2),

    room_width=9,
    room_width_depth_ratio=dict(dist='uniform_', lbound=0.4, ubound=0.6),
    room_height=3,

    # fire_mode=dict(dist='discrete_', values='0,1,2', weights='0.07,0.59,0.34'),
    fire_mode=0,
    fire_fuel_density_MJm2=dict(dist="br187_fuel_load_density_", lbound=10, ubound=5000),
    fire_combustion_efficiency=dict(dist='uniform_', lbound=0.8, ubound=1.0),
    fire_hrr_density_kWm2=dict(dist='br187_hrr_density_'),
    fire_growth_factor=0.0117,
    fire_t_lim=20. * 60.,
    fire_convection_factor=0.7,

    detector_to_fire_vertical_distance=3 - 0.02,
    detector_to_fire_horizontal_distance=2.83,
    detector_act_temp=93 + 273.15,
    detector_response_time_index=250,
    detector_conduction_factor=0.65,

    lining_rho=2000,
    lining_c=1000,
    lining_k=1.13,

    receiver_ignition_temperature=-1,
    receiver_emissivity=1,
    receiver_separation=6.586120496273041,

    ftp_chf=13.4e3,
    ftp_index=2.0,
    ftp_target=34592,
)

EXAMPLE_INPUT = dict(
    CASE_1=dict(
        n_simulations=1000,
        t_end=180. * 60.,
        t_step=10.,

        opening_width=9,
        opening_height=3,
        opening_ventilation_fraction=dict(dist='lognorm_mod_', lbound=1e-4, ubound=1 - 1e-4, mean=0.2, sd=0.2),

        room_width=9,
        room_width_depth_ratio=dict(dist='uniform_', lbound=0.4, ubound=0.6),
        room_height=3,

        # fire_mode=dict(dist='discrete_', values='0,1,2', weights='0.07,0.59,0.34'),
        fire_mode=0,
        fire_fuel_density_MJm2=dict(dist="br187_fuel_load_density_", lbound=10, ubound=5000),
        fire_combustion_efficiency=dict(dist='uniform_', lbound=0.8, ubound=1.0),
        fire_hrr_density_kWm2=dict(dist='br187_hrr_density_'),
        fire_growth_factor=0.0117,
        fire_t_lim=20. * 60.,
        fire_convection_factor=0.7,

        detector_to_fire_vertical_distance=3 - 0.02,
        detector_to_fire_horizontal_distance=2.83,
        detector_act_temp=93 + 273.15,
        detector_response_time_index=250,
        detector_conduction_factor=0.65,

        lining_rho=2000,
        lining_c=1000,
        lining_k=1.13,

        receiver_ignition_temperature=-1,
        receiver_emissivity=1,
        receiver_separation=6.586120496273041,

        ftp_chf=13.4e3,
        ftp_index=2.0,
        ftp_target=34592,
    )
)

if __name__ == "__main__":
    print(EXAMPLE_INPUT_DETERMINISTIC, "\n")
