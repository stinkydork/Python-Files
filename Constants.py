import numpy as np
# This is the file for storing some data and constants

G = 6.67430e-20 # km^3 / kg^1.s^2

Earth = {
    'name' : 'Earth',
    'mass' : 5.9721e24, # kg
    'mu' : 398600.4418, # km^3 / s^2
    'radius' : 6378.137, # km
    'J2' : 1.082635854e-3,
    'J3' : -2.5326613168e-6,
    'at_rot_vec' : np.array([0,0,7.29211e-5]) # rad / s
}

# 1U CubeSat (Made-up)
CubeSat_1U = {
    'name' : 'CubeSat_1U',
    'mass' : 1.2, # kg
    'shape' : 'cube', 
    'sigma_N' : 1,
    'sigma_T' : 1,
    'T_wall' : 300, # k
    'max_slew_rate' : (3 * (np.pi/180)), # rad/s
    'lx' : 0.1, # m
    'ly' : 0.1, # m
    'lz' : 0.1, # m
    'covarience_matrix' : np.diag([0.1, 0.1, 0.1]), # km^2
}
# Made-up Data
Debri = {
    'name' : 'debri',
    'mass' : 25, # kg
    'shape' : 'sphere',
    'radius' : 1, # m 
    'covarience_matrix' : np.diag([0.3, 0.3, 0.3]) # km^2
}

# TERP Raptor
CubeSat_TerpRaptor = {
    'name' : 'CubeSat_TerpRaptor',
    'mass' : 30.4, # kg
    'shape' : 'cube', 
    'sigma_N' : 1,
    'sigma_T' : 1,
    'T_wall' : 300, # k
    'max_slew_rate' : (3 * (np.pi/180)), # rad/s
    'drag_ref_area' : 0.12, # m^2
    'SRP_area' : 0.12, # m^2
    'HBR_radius' : 0.3, # m
    'covarience_matrix' : np.diag([0.1, 0.1, 0.1]), # km^2
}

# Made-up Data (Based on Starlink v2 mini)
Satellite = {
    'name' : 'Satellite',
    'mass' : 800, # kg
    'Cd' : 2.2, # (Not Sure)
    'Drag_Area' : 1.08e-5, # km^2
    'SRP_Area' : 1.08e-5, # km^2
    'HBR_Radius' : 0.004, # km
    'Covarience_Matrix' : np.diag([0.1, 0.1, 0.1]), # km^2
    # Source - https://www.space-propulsion.com/spacecraft-propulsion/hydrazine-thrusters/20n-hydrazine-thruster.html
    'Impluse' : 230, # seconds
    'Thrust' : 0.02 # kN
}

# Starlink v2 mini (Guesstimated)
# Source 1 - https://dishycentral.com/how-big-are-starlink-satellites
# Source 2 - https://www.aviationtoday.com/2023/02/28/spacex-shares-details-higher-bandwidth-v2-mini-starlink-satellites/
Starlink_v2_mini = {
    'name' : 'Starlink_v2_mini',
    'mass' : 800, # kg
    'Cd' : 2.2, # (Not Sure)
    'Drag_Area' : 1.08e-5, # km^2 
    'HBR_Radius' : 0.004, # km
    'Covarience_Matrix' : np.diag([0.1, 0.1, 0.1]), # km^2
    'Low_Impluse' : 2500, # seconds
    'Low_Thrust' : 0.00017 # kN
}