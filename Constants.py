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
    'mass' : 29.6, # kg
    'shape' : 'cube', 
    'sigma_N' : 1,
    'sigma_T' : 1,
    'T_wall' : 300, # k
    'max_slew_rate' : (3 * (np.pi/180)), # rad/s
    'lx' : 0.4, # m
    'ly' : 0.2, # m
    'lz' : 0.2, # m
    'covarience_matrix' : np.diag([0.1, 0.1, 0.1]), # km^2
}