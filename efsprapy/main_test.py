#  Test file
from efsprapy.main import main

if __name__ == "__main__":
    ftp = main(
        N=1000,  # Number of simulations
        fire_duration=120,  # fire duration [min]
        time_step=1,  # time step [sec]
        h_opening=2.5,  # Height of opening [m]
        area_opening=25,  # Area of opening [m²]
        area_floor=100,  # Floor area [m²]
        area_total=300,  # Total compartment area [m²]
        fire_load_density_mean=420,  # Mean fire load density [MJ/m²]
        comb_eff=0.8,  # Combustion efficiency
        lining_rho=850,  # Density of compartment lining [kg/m³]
        lining_cp=1050,  # Specific heat of compartment [J/kgK]
        lining_k=0.2,  # Thermal conductivity of compartment [W/mK]
        limiting_time=20,  # Time to reach maximum temperature for fuel controlled fire [min]
        panel_width=20,  # Width of radiating panel [m]
        panel_height=15,  # Height of radiating panel [m]
        boundary_distance=10,  # Emitter to boundary [m]):
        heat_flux_critical=12.6,  # Critical heat flux of receiver material [kW/m²]
        n_factor=1,  # FTP index [n_factor = 1 (for thermally thin); n_factor = 2 (for thermally thick); n_factor = 1.5 (for intermediate)]
    )
    print(ftp)
