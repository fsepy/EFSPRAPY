EXAMPLE_INPUT_DETERMINISTIC = dict(
    t_end=180. * 60.,
    t_step=5.,
    opening_width=31 / 2.,  # 31/1.565
    opening_height=2.,  # 1.565
    room_height=3,
    room_floor_area=150,
    room_total_surface_area=500,
    fire_mode=1,
    fire_fuel_density_MJm2=420,  # 420
    fire_hrr_density_kWm2=250,
    fire_growth_factor=0.0117,
    fire_t_lim=20. * 60.,
    fire_din_growth_factor=300,
    fire_convection_factor=0.7,
    detector_to_fire_vertical_distance=3,
    detector_to_fire_horizontal_distance=2.1,
    detector_act_temp=93 + 273.15,
    detector_response_time_index=250,
    detector_conduction_factor=0.4,
    lining_rho=2000,
    lining_c=1000,
    lining_k=1.13,
    lining_thermal_effusivity=1500,
    receiver_ignition_temperature=300 + 273.15,
    receiver_emissivity=1,
    receiver_separation=2,
    ftp_chf=12.6e3,
    ftp_index=1,
    ftp_target=10000,
)

EXAMPLE_INPUT = dict(
    CASE_1=dict(
        n_simulations=100,
        t_end=180. * 60.,
        t_step=5.,

        opening_width=31 / 1.565,
        opening_height=1.565,

        room_height=3,
        room_floor_area=150,
        room_total_surface_area=500,

        fire_mode=dict(dist='discrete_', values='0,1,2', weights='0.07,0.59,0.34'),
        fire_fuel_density_MJm2=dict(dist="gumbel_r_", lbound=200, ubound=1200, mean=780, sd=234),  # 420
        fire_hrr_density_kWm2=250,
        fire_growth_factor=0.0117,
        fire_t_lim=20. * 60.,
        fire_din_growth_factor=300,
        fire_convection_factor=0.7,

        detector_to_fire_vertical_distance=3,
        detector_to_fire_horizontal_distance=2.1,
        detector_act_temp=93 + 273.15,
        detector_response_time_index=250,
        detector_conduction_factor=0.4,

        lining_rho=2000,
        lining_c=1000,
        lining_k=1.13,
        lining_thermal_effusivity=1500,

        receiver_ignition_temperature=300 + 273.15,
        receiver_emissivity=1,
        receiver_separation=2,

        ftp_chf=12.6e3,
        ftp_index=1,
        ftp_target=10000,
    )
)

if __name__ == "__main__":
    print(EXAMPLE_INPUT_DETERMINISTIC, "\n")
