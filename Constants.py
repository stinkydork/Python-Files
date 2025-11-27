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

# Made-up CubeSat
CubeSat = {
    'name' : 'CubeSat',
    'mass' : 1.5, # kg
    'Cd' : 2.2, # (Not Sure)
    'Drag_Area' : 9e-8, # km^2
    'SRP_Area' : 9e-8, # km^2
    'HBR_Radius' : 0.0003, # km
    'Covarience_Matrix' : np.diag([0.1, 0.1, 0.1]), # km^2
    # Source - https://www.space-propulsion.com/spacecraft-propulsion/hydrazine-thrusters/20n-hydrazine-thruster.html
    'Impluse' : 230, # seconds
    'Thrust' : 0.02 # kN
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

# Made-up Data
Debri = {
    'name' : 'Satellite',
    'mass' : 25, # kg
    'Cd' : 2.2,
    'Drag_Area' : 1e-7, # km^2
    'HBR_Radius' : 0.1, # km
    'Covarience_Matrix' : np.diag([3, 3, 3]) # km^2
}