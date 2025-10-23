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

Satellite = {
    'name' : 'Satellite',
    'mass' : 25, # kg
    'Cd' : 2.2,
    'Drag_Area' : 1e-7, # km^2
    'HBR_Radius' : 0.01, # km
    'Covarience_Matrix' : np.diag([1, 1, 1]) # km^2
}

Debri = {
    'name' : 'Satellite',
    'mass' : 25, # kg
    'Cd' : 2.2,
    'Drag_Area' : 1e-7, # km^2
    'HBR_Radius' : 0.1, # km
    'Covarience_Matrix' : np.diag([3, 3, 3]) # km^2
}