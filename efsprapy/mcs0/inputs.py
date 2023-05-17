EXAMPLE_INPUT_DETERMINISTIC = dict(
    t_end=120. * 60.,
    t_step=1.,
    vent_width=31 / 1.565,
    vent_height=1.565,
    room_floor_area=150,
    room_total_surface_area=500,
    fuel_density=420,
    rho=2000,
    c=1000,
    k=1.13,
    t_lim=15,
    emissivity=1,
    emitter_width=20,
    emitter_height=30,
    emitter_receiver_separation=10,
    chf=12.6e3,
    ftp_index=1,
)

EXAMPLE_INPUT = dict(
    CASE_1=dict(
        n_simulations=1000,
        t_end=120. * 60.,
        t_step=1.,
        vent_width=31 / 1.565,
        vent_height=1.565,
        room_floor_area=150,
        room_total_surface_area=500,
        fuel_density=dict(dist="gumbel_r_", lbound=10, ubound=1200, mean=780, sd=234),
        rho=2000,
        c=1000,
        k=1.13,
        t_lim=15,
        emissivity=1,
        emitter_width=20,
        emitter_height=30,
        emitter_receiver_separation=10,
        chf=12.6e3,
        ftp_index=1,
    )
)

if __name__ == "__main__":
    print(EXAMPLE_INPUT_DETERMINISTIC, "\n")
