# EFS functions
import matplotlib.pyplot as plt
import numpy as np
from fsetools.lib.fse_bs_en_1991_1_2_parametric_fire import temperature as param_fire


# ISO fire curve
def iso_fire(
        fire_duration: float,  # Fire duration in minutes
        t_step: float  # Time step in sec
):
    time = np.arange(0, fire_duration * 60 + t_step, t_step)
    fire_temp = 20 + 345 * np.log10(8 * (time / 60) + 1)  # [deg C]
    return [time, fire_temp]


def en_parametric_fire(
        fire_duration: float,  # fire duration [min]
        time_step: float,  # time step [sec]
        h_opening: float,  # Height of opening [m]
        area_opening: float,  # Area of opening [m²]
        area_floor: float,  # Floor area [m²]
        area_total: float,  # Total compartment area [m²]
        fire_load_density: float,  # Fire load density [MJ/m²]
        comb_eff: float,  # Combustion efficiency
        lining_rho: float,  # Density of compartment lining [kg/m³]
        lining_cp: float,  # Specific heat of compartment [J/kgK]
        lining_k: float,  # Thermal conductivity of compartment [W/mK]
        limiting_time: float):  # Time to reach maximum temperature for fuel controlled fire [min]

    time = np.arange(0, fire_duration * 60 + time_step, time_step, dtype=float)  # Time in [sec]
    time_min = time / 60  # Time in [min]
    time_hr = time / 3600  # Time in [hr]
    ambient_temperature = 20  # Ambient temperature [deg C]

    # Calculation of compartment boundary thermal inertia (b-factor)
    b_factor = np.sqrt(lining_rho * lining_cp * lining_k)  # Thermal inertia of compartment walls [J/m^(2)s^(0.5)K]
    if b_factor < 100:
        b_factor = 100
    elif b_factor > 2200:
        b_factor = 2200
    # b_factor = 1849.1
    #  Calculation of opening factor
    O_factor = area_opening * np.sqrt(h_opening) / area_total  # Opening factor [m^(0.5)]
    if O_factor < 0.02:
        O_factor = 0.02
    elif O_factor > 0.2:
        O_factor = 0.2

    #  Calculation of design fire load density per total compartment area
    q_fd = comb_eff * fire_load_density
    q_td = q_fd * area_floor / area_total
    if q_td < 50:
        q_td = 50
    elif q_td > 1000:
        q_td = 1000

    #  Calculation of intermediate parameters

    Tau = ((O_factor / b_factor) / (0.04 / 1160)) ** 2
    t_lim = limiting_time / 60  # Time at maximum temperature for fuel controlled fires [hours]
    t_max = np.amax([0.2 * (10 ** -3) * q_td / O_factor, t_lim])

    O_factor_lim = 0.1 * (10 ** -3) * q_td / t_lim
    Tau_lim = ((O_factor_lim / b_factor) / (0.04 / 1160)) ** 2

    if O_factor > 0.04 and q_td < 75 and b_factor < 1160:
        Tau_lim = Tau_lim * (1 + ((O_factor - 0.04) / 0.04) * ((q_td - 75) / 75) * ((1160 - b_factor) / 1160))

    if t_max >= t_lim:  # Ventilation controlled fire
        t_x = time_hr * Tau
        t_x_max = t_max * Tau
        x_factor = 1.0
        fire_type = "Ventilation controlled fire"
    else:  # Fuel controlled fire
        t_x = time_hr * Tau_lim
        t_x_max = t_lim * Tau
        x_factor = t_lim / t_max
        fire_type = "Fuel controlled fire"

    T_max = 20 + 1325 * (1 - 0.324 * np.exp(-0.2 * t_x_max) - 0.204 * np.exp(-1.7 * t_x_max)
                         - 0.472 * np.exp(-19 * t_x_max))

    # Calculation of fire temperature

    fire_temp = np.zeros_like(time)  # Initialization of variables

    # Cooling phase parameters
    t_x_cooling = Tau * time_hr
    t_x_max_cooling = Tau * t_max

    for t in range(0, len(time), 1):
        if time_hr[t] <= t_max:  # Heating phase
            fire_temp[t] = 20 + 1325 * (1 - 0.324 * np.exp(-0.2 * t_x[t]) - 0.204 * np.exp(-1.7 * t_x[t])
                                        - 0.472 * np.exp(-19 * t_x[t]))
        else:
            if t_x_max <= 0.5:
                fire_temp[t] = T_max - 625 * (t_x_cooling[t] - t_x_max_cooling * x_factor)
            elif 0.5 < t_x_max < 2:
                fire_temp[t] = T_max - 250 * (3 - t_x_max_cooling) * (t_x_cooling[t] - t_x_max_cooling * x_factor)
            elif t_x_max >= 2:
                fire_temp[t] = T_max - 250 * (t_x_cooling[t] - t_x_max_cooling * x_factor)

        if fire_temp[t] < ambient_temperature:
            fire_temp[t] = ambient_temperature

    return [fire_temp]


def din_parametric_fire(
        fire_duration: float,  # fire duration [min]
        time_step: float,  # time step [sec]
        h_opening: float,  # Height of opening [m]
        area_opening: float,  # Area of opening [m²]
        area_floor: float,  # Floor area [m²]
        area_total: float,  # Total compartment area [m²]
        fire_load_density: float,  # Design fire load density [MJ/m²]
        comb_eff: float,  # Combustion efficiency
        lining_rho: float,  # Density of compartment lining [kg/m³]
        lining_cp: float,  # Specific heat of compartment [J/kgK]
        lining_k: float):  # Thermal conductivity of compartment [W/mK]

    time = np.arange(0, fire_duration * 60 + time_step, time_step, dtype=float)
    time_min = time / 60  # Time [min]
    ambient_temperature = 20  # Ambient temperature [deg C]

    b_factor = np.sqrt(lining_rho * lining_cp * lining_k)  # Thermal inertia of compartment walls [J/m^(2)s^(0.5)K]
    if b_factor < 100:
        b_factor = 100
    elif b_factor > 2200:
        b_factor = 2200

    A_total_less_opening = area_total - area_opening  # Total area of compartment WITHOUT opening [m^(2)]
    O_factor = area_opening * np.sqrt(h_opening) / area_total  # Opening factor [m^(0.5)]
    if O_factor < 0.02:
        O_factor = 0.02
    elif O_factor > 0.2:
        O_factor = 0.2

    # flf = 1  # Fire load factor
    # q_fk = 511  # Design value of the fire load density related to the surface area of the floor [MJ/m^(2)]
    q_fd = comb_eff * fire_load_density  # Design value of fire load density related to floor area [MJ/m^(2)]

    # Definition of reference fire load case parameters

    q_f_ref = 1300  # Reference fire load density related to the surface area of the floor [MJ/m^(2)]
    Q_f_ref = q_f_ref * area_floor  # Reference total fire load [MJ]
    h_u = 17.5  # Net Calorific value if wood [MJ/kg]
    t_g = 300  # Time constant for fire growth rate [sec]
    rhr_f = 250  # Maximum fire load density [kW/m^(2)]
    rhr_o = 1.0  # RHR at time Tg [MW]
    # rhr_flashover = 0.0078 * A_total_less_opening + 0.378 * area_opening * np.sqrt(h_opening)  # RHR at flashover [MW]
    # t_flashover = np.sqrt((t_g ** 2) * rhr_flashover)  # Time to flashover [sec]
    rhr_max_fuel = (rhr_f / 1000) * area_floor  # Maximum RHR for fuel controlled fire [MW]
    rhr_max_vent = 0.1 * comb_eff * h_u * area_opening * np.sqrt(h_opening)  # Maximum RHR for fuel controlled fire [MW]
    rhr_max = np.minimum(rhr_max_fuel, rhr_max_vent)  # Maximum RHR [MW]
    # rhr_growth = rhr_o * ((time / t_g) ** 2)  # RHR for fire growth stage

    t1_ref = t_g * np.sqrt(rhr_max)
    Q1_ref = rhr_o * (t1_ref ** 3) / (3 * (t_g ** 2))
    Q2_ref = 0.7 * Q_f_ref - Q1_ref
    t2_ref = t1_ref + Q2_ref / rhr_max
    Q3_ref = Q_f_ref - (Q1_ref + Q2_ref)
    t3_ref = t2_ref + 2 * Q3_ref / rhr_max

    k = ((rhr_max ** 2) / (area_opening * np.sqrt(h_opening) * A_total_less_opening * b_factor)) ** (1 / 3)
    if rhr_max == rhr_max_fuel:
        fire_type = 'FUEL CONTROLLED'
        if k <= 0.04:
            T1_ref = 24000 * k + 20
            T2_ref = 33000 * k + 20
            T3_ref = 16000 * k + 20
        else:
            T1_ref = 980
            T2_ref = 1340
            T3_ref = 660
    else:
        fire_type = 'VENTILATION CONTROLLED'
        T1_ref = -8.75 / O_factor - 0.1 * b_factor + 1175
        T2_ref = (0.04 * b_factor - 17) / O_factor - 0.4 * b_factor + 2175
        if T2_ref > 1340:
            T2_ref = 1340
        T3_ref = -5 / O_factor - 0.16 * b_factor + 1060

    # Definition of particular fire load case parameters
    Q_f = q_fd * area_floor
    Q1 = Q1_ref
    t1 = t1_ref
    T1 = T1_ref
    if Q1_ref < 0.7 * Q_f:
        Q2 = 0.7 * Q_f - Q1
        t2 = t1_ref + (0.7 * Q_f - rhr_o * ((t1_ref ** 3) / (3 * (t_g ** 2)))) / rhr_max
        T2 = (T2_ref - T1_ref) * np.sqrt((t2 - t1_ref) / (t2_ref - t1_ref)) + T1_ref
    else:
        Q1 = 0.7 * Q_f
        Q2 = 0
        t1 = (0.7 * Q_f * 3 * (t_g ** 2)) ** 3
        t2 = t1
        T1 = (T1_ref - ambient_temperature) * ((t1 / t1_ref) ** 2) + ambient_temperature
        T2 = rhr_o * (t2 ** 2) / (t_g ** 2)
    Q3 = Q_f - (Q1 + Q2)
    t3 = t2 + 0.6 * Q_f / rhr_max
    T3 = (T3_ref / np.log10(t3_ref / 60 + 1)) * np.log10(t3 / 60 + 1)

    # HRR Curve
    rhr = np.zeros_like(time)  # Initialization of variables
    for t in range(0, len(time), 1):
        if time[t] <= t1:
            rhr[t] = rhr_o * ((time[t] / t_g) ** 2)
        elif t1 < time[t] <= t2:
            rhr[t] = rhr_max
        elif t2 < time[t] <= t3:
            rhr[t] = rhr_max * (t3 - time[t]) / (t3 - t2)
        else:
            rhr[t] = 0

    # Fire temperature time curve [℃]
    fire_temp = np.zeros_like(time)  # Initialization of variables
    for t in range(0, len(time), 1):
        if time[t] <= t1:
            fire_temp[t] = (T1 - ambient_temperature) * ((time[t] / t1) ** 2) + ambient_temperature
        elif t1 < time[t] <= t2:
            fire_temp[t] = (T2 - T1) * np.sqrt((time[t] - t1) / (t2 - t1)) + T1
        elif time[t] > t2 and T2 > T3:
            fire_temp[t] = (T3 - T2) * np.sqrt((time[t] - t2) / (t3 - t2)) + T2
        else:
            fire_temp[t] = (T2 - T3) * np.sqrt((time[t] - t2) / (t3 - t2)) + T2

        if fire_temp[t] < ambient_temperature:
            fire_temp[t] = ambient_temperature
    """
        # RHR and fire temperature curves
        fig, rhr_curve = plt.subplots()
        fig.suptitle('Plot of RHR and fire temperature against time', fontsize=18)
    
        rhr_curve.set_xlabel('time [sec]', fontsize=14)
        rhr_curve.set_ylabel('Rate of Heat Release [MW]', color="red", fontsize=14)
        rhr_curve.plot(time / 60, rhr, color="red")
        rhr_curve.tick_params(axis='y', labelcolor="red")
    
        fire_temp_curve = rhr_curve.twinx()
        fire_temp_curve.set_ylabel('Temperature [℃]', color="blue", fontsize=14)
        fire_temp_curve.plot(time / 60, fire_temp, color="blue")
        fire_temp_curve.tick_params(axis='y', labelcolor="blue")
    
        fig.tight_layout()
        plt.show()
    """
    return [fire_temp]


def travelling_fire(
        fire_duration: float,  # fire duration [min]
        time_step: float,  # time step [sec]
        fire_load_density: float,  # Design fire load density [MJ/m²]
        comb_eff: float,  # Combustion efficiency
        HRRPUA: float,  # HRR per unit area [kW/m²]
        alpha: float,  # Fire growth rate [kW/s²]
        L: float,  # Length of compartment [m]
        W: float,  # Width of compartment [m]
        spread_rate: float,  # Spread rate [m/s]
        H: float,  # Height of compartment [m]
        h_opening: float,  # Height of opening [m]
        A_opening: float,  # Area of opening [m²]
        Hu: float,  # Heat of combustion of wood [MJ/kg]
        x: float):  # Point in space of interest (beam length) [m]

    time = np.arange(0, fire_duration * 60 + time_step, time_step, dtype=float)  # Time [sec]
    q_fd = comb_eff * fire_load_density
    A_floor = L * W  # compartment floor area [m²]
    A_total = 2 * (A_floor + L * H + W * H)
    O_factor = A_opening * np.sqrt(h_opening) / A_total
    #  O_factor = 0.01

    #  Fire spread rate
    s = np.sqrt(2 * alpha / (HRRPUA * np.pi))
    s = spread_rate

    #  Calculation of burning duration/stages
    t_burn = max(q_fd * 1000 / HRRPUA, 900.0)
    t_decay = max(t_burn, L / s)
    t_lim = min(t_burn, L / s)

    t_steady = q_fd * 1000 / HRRPUA
    t_max = L / s

    #  Calculation of HRR
    Q_growth = HRRPUA * W * s * time
    Q_max = min(HRRPUA * W * s * t_burn, HRRPUA * W * L)
    #  Q_decay = max((Q_max - (time - t_decay) * W * s * HRRPUA), 0)

    Q = np.zeros_like(time)  # Initialization of variables

    for t in range(0, len(time), 1):
        if time[t] < t_steady:
            Q[t] = min(alpha * (time[t] ** 2), HRRPUA * W * s * time[t], 0.1 * comb_eff * Hu * A_opening * (np.sqrt(h_opening)))
        elif t_steady <= time[t] < t_max:
            Qs = 0.5 * HRRPUA * ((np.pi * ((s * time[t]) ** 2)) - (np.pi * ((s * time[t] - s * t_steady) ** 2)))
            Q[t] = min(Qs, HRRPUA * W * s * t_steady, 0.1 * comb_eff * Hu * A_opening * (np.sqrt(h_opening)))
        else:
            Q[t] = HRRPUA * W * (s * t_steady - s * (time[t] - t_max))
        if Q[t] < 0:
            Q[t] = 0

    for t in range(0, len(time), 1):
        if time[t] < t_lim:
            Q[t] = HRRPUA * W * s * time[t]
        elif t_lim <= time[t] < t_decay:
            Q[t] = min(HRRPUA * W * s * t_burn, HRRPUA * W * L)
        else:
            Q[t] = max(HRRPUA * W * (s * t_steady - s * (time[t] - t_max)), 0)
        if Q[t] < 0:
            Q[t] = 0

    # Calculation of temperature [℃]
    fire_temp = np.zeros_like(time)  # Initialization of variables
    T_ambient = 20  # Ambient temperature [℃]

    if O_factor <= 0.02:
        T_near_field = 1200
    elif 0.02 < O_factor < 0.2:
        T_near_field = 1200 - (O_factor - 0.02) * (1200 - 900) / (0.2 - 0.02)
    else:
        T_near_field = 900
    # T_near_field = 1200

    for t in range(0, len(time), 1):
        if time[t] < t_max:
            r = np.abs(x - s * time[t])
        else:
            r = np.abs(x - s * t_max)

        if r / H > 0.18:
            fire_temp[t] = 5.38 * np.power(Q[t] / r, 2 / 3) / H + T_ambient
        else:
            fire_temp[t] = 16.9 * np.power(Q[t], 2 / 3) / np.power(H, 5 / 3) + T_ambient

        if fire_temp[t] > T_near_field:
            fire_temp[t] = T_near_field
        if fire_temp[t] < T_ambient:
            fire_temp[t] = T_ambient

    return [fire_temp]


def travelling_fire_2(
        fire_duration: float,  # fire duration [min]
        time_step: float,  # time step [sec]
        fire_load_density: float,  # Design fire load density [MJ/m²]
        comb_eff: float,  # Combustion efficiency
        HRRPUA: float,  # HRR per unit area [kW/m²]
        L: float,  # Length of compartment [m]
        W: float,  # Width of compartment [m]
        H: float,  # Height of compartment [m]
        spread_rate: float,  # Spread rate [m/s]
        T_near_field: float,  # Near field temperature [℃]
        x: float):  # Point in space of interest (beam length) [m]

    time = np.arange(0, fire_duration * 60 + time_step, time_step, dtype=float)  # Time [sec]
    m = 0.8  # Combustion efficiency
    q_fd = comb_eff * fire_load_density
    A_floor = L * W  # compartment floor area [m²]

    #  Fire spread rate
    s = spread_rate

    #  Calculation of burning duration/stages
    t_burn = max(q_fd * 1000 / HRRPUA, 900.0)
    t_decay = max(t_burn, L / s)
    t_lim = min(t_burn, L / s)

    #  Calculation of HRR
    Q_growth = HRRPUA * W * s * time
    Q_max = min(HRRPUA * W * s * t_burn, HRRPUA * W * L)
    Q_decay = Q_max - (time - t_decay) * W * s * HRRPUA

    Q = np.zeros_like(time)  # Initialization of variables

    for t in range(0, len(time), 1):
        if time[t] < t_burn:
            Q[t] = Q_growth[t]
        elif t_lim <= time[t] < t_decay:
            Q[t] = Q_max
        elif time[t] > t_burn:
            Q[t] = Q_decay[t]
        if Q[t] < 0:
            Q[t] = 0

    # Calculation of temperature
    fire_temp = np.zeros_like(time)  # Initialization of variables
    T_ambient = 20  # Ambient temperature [℃]

    for t in range(0, len(time), 1):

        l_fire_front = s * time[t]
        if l_fire_front < 0:
            l_fire_front = 0.0
        if l_fire_front > L:
            l_fire_front = L
        l_fire_end = s * (time[t] - t_lim)
        if l_fire_end < 0:
            l_fire_end = 0.0
        if l_fire_end > L:
            l_fire_end = L
        l_fire_median = (l_fire_front + l_fire_end) / 2.0
        r = np.abs(x - l_fire_median)

        if r / H > 0.18:
            fire_temp[t] = 5.38 * np.power(Q[t] / r, 2 / 3) / H + T_ambient
        else:
            fire_temp[t] = 16.9 * np.power(Q[t], 2 / 3) / np.power(H, 5 / 3) + T_ambient

        if fire_temp[t] > T_near_field:
            fire_temp[t] = T_near_field
        if fire_temp[t] < T_ambient:
            fire_temp[t] = T_ambient

    return [fire_temp]


# Time equivalence
def time_eq(
        q_fk: float,  # Characteristic fire load density [MJ/m²]
        m: float,  # Combustion factor
        delta_1: float,  # Sprinkler multiplication factor
        k_b: float,  # Conversion factor in sec
        m_risk: float,  # Risk multiplication factor
        A_f: float,  # Floor area [m²]
        A_v: float,  # Area of vertical opening [m²]
        A_h: float,  # Area of horizontal opening [m²]
        h_comp: float  # Height of compartment [m]
):
    q_fd = q_fk * m * delta_1
    alpha_v = A_v / A_f
    if alpha_v < 0.0025:
        alpha_v = 0.0025
    if alpha_v > 0.25:
        alpha_v = 0.25

    alpha_h = A_h / A_f

    b_v = 12.5 * (1 + 10 * alpha_v - (alpha_v ** 2))
    if b_v < 10:
        b_v = 10

    w_f = ((6 / h_comp) ** 0.3) * (0.62 + 90 * ((0.4 - alpha_v) ** 4)) / (1 + b_v * alpha_h)
    if w_f < 0.5:
        w_f = 0.5

    t_ed = (q_fd * k_b * w_f) * m_risk
    return [alpha_v, alpha_h, b_v, w_f, t_ed]


def heat_flux_source(
        temp: float,  # Temperature of radiating object [℃]
        emissivity: float):  # Emissivity of radiating object
    boltz = 5.67 * (10 ** -11)
    heat_flux = boltz * emissivity * ((temp + 273.15) ** 4)
    return heat_flux


def view_factor_parallel(
        width: float,  # Width of panel [m]
        height: float,  # Height of panel [m]
        boundary_distance: float):  # Emitter to boundary [m]):

    # Renaming of variables for ease of use in calculations
    W = width
    H = height
    s_half = boundary_distance

    # Calculation of view_factor
    s = 2 * s_half  # Separation distance [m]
    X = W / (2 * s)
    Y = H / (2 * s)
    View_factor = (2 / np.pi) * ((X / np.sqrt(1 + (X ** 2))) * np.arctan((Y / np.sqrt(1 + (X ** 2))))
                                 + (Y / np.sqrt(1 + (Y ** 2))) * np.arctan((X / np.sqrt(1 + (Y ** 2)))))
    return View_factor


# 1-D heat transfer solver
def calculate_1d_temperature(num_points, time_step, num_steps, heat_flux, density, specific_heat, thermal_conductivity):
    # Constants
    ambient_temperature = 20  # °C

    # Initialize the temperature array
    temperatures = np.zeros(num_points)

    # Set the initial condition
    temperatures[0] = ambient_temperature

    # Calculate the temperature at each time step
    for step in range(num_steps):
        # Calculate the temperature change due to conduction
        conduction_term = (
                thermal_conductivity * np.diff(temperatures, 2) / (density * specific_heat)
        )

        # Calculate the temperature change at the boundary due to the heat flux
        boundary_term = heat_flux / (density * specific_heat)

        # Update the temperature array
        temperatures[1] += (conduction_term[0] + boundary_term) * time_step

        # Update the temperature at the boundary
        temperatures[0] += boundary_term * time_step

    return temperatures


# Calculate critical flux time product for a give incidence heat flux
def flux_time_product_critical(
        heat_flux: float,  # Heat flux received [kw/m²]
        heat_flux_crit: float,  # Critical heat flux received [kw/m²]
        ignition_time: float,  # time to ignition of receive material [sec]
        n_factor: float,  # FTP index [n_factor = 1 (for thermally thin); n_factor = 2 (for thermally thick); n_factor = 1.5 (for intermediate)]
):
    ftp_critical = ((heat_flux - heat_flux_crit) ** n_factor) * ignition_time
    return ftp_critical


# Calculate flux time product for a variable heat flux received
def flux_time_product(
        heat_flux_variable: list,  # Heat flux received [kw/m²]
        heat_flux_crit: float,  # Critical heat flux received [kw/m²]
        time: list,  # time to ignition of receive material [sec]
        time_step,  # time to ignition of receive material [sec]
        n_factor: float,  # FTP index [n_factor = 1 (for thermally thin); n_factor = 2 (for thermally thick); n_factor = 1.5 (for intermediate)]
):
    ftp = []
    ftp[0] = 0
    ftp[1] = ((heat_flux_variable[1] - heat_flux_crit) ** n_factor) * time_step
    if heat_flux_variable[1] <= heat_flux_crit:
        ftp[1] = 0
    for i in range(2, time_step, time[-1]):
        ftp[i] = ((heat_flux_variable[i] - heat_flux_crit) ** n_factor) * time_step + ftp[i-1]
        if heat_flux_variable[i] <= heat_flux_crit:
            ftp[i] = 0
    return ftp


# Probability distributions

def gumbel_dist(
        mean: float,  # Steel temperature [℃]
        std: float,  # Random factor
        N: int,  # Number of simulations
):
    beta = std * np.sqrt(6) / np.pi
    mu = mean - 0.5772 * beta
    ran_var = np.random.gumbel(mu, beta, N)  # stats.gumbel_r.rvs(420, 126, N)

    return ran_var


def lognormal_dist(
        mean: float,  # Steel temperature [℃]
        std: float,  # Random factor
        N: int,  # Number of simulations
):
    V = std ** 2  # Variance
    mu = np.log((mean ** 2) / np.sqrt(V + (mean ** 2)))  # Location parameter
    sigma = np.sqrt(np.log(V / (mean ** 2) + 1))  # Scale parameter
    ran_var = np.random.lognormal(mu, sigma, N)  # Generate random variable

    return ran_var


def trunc_lognormal_dist(
        mean: float,  # Steel temperature [℃]
        std: float,  # Random factor
        lim_value: float,  # Limiting value for truncation
        N: int,  # Number of simulations
):
    V = std ** 2  # Variance
    mu = np.log((mean ** 2) / np.sqrt(V + (mean ** 2)))  # Location parameter
    sigma = np.sqrt(np.log(V / (mean ** 2) + 1))  # Scale parameter
    ran_var = np.zeros_like(range(0, N), dtype=float)
    i = 0
    while i < N:
        ran_var[i] = np.random.lognormal(mu, sigma, 1)  # Generate random variable
        if ran_var[i] <= lim_value:
            i = i + 1
        else:
            i = i

    return ran_var


def gamma_dist(
        mean: float,  # Steel temperature [℃]
        std: float,  # Random factor
        N: int,  # Number of simulations
):
    alpha = (mean / std) ** 2  # Shape parameter
    beta = (std ** 2) / mean  # Scale parameter
    ran_var = np.random.gamma(alpha, beta, N)

    return ran_var
