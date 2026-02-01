import numpy as np
from OrbitPropagator import OrbitPropagator, null_perturbations
import Constants as ct
import time

# Use this file for copying the syntax.
if __name__ == '__main__':
    # Initial conditions of satellite
    r0_S, v0_S = OrbitPropagator.Kep2State(7000, 0, 51.6, 0, 0, 0, cb=ct.Earth)

    print(r0_S)
    print(v0_S)
    r0_D = np.array([7000+(1e-6), 1e-6, 1e-6])
    v0_D = np.array([0, 4.687, 5.913])


    collision_Prob = OrbitPropagator.Pc2D(r0_S, v0_S, ct.CubeSat['Covarience_Matrix'], r0_D, v0_D, ct.Debri['Covarience_Matrix'], (ct.CubeSat['HBR_Radius']+ct.Debri['HBR_Radius']))
    print(collision_Prob)
    print(1e-6*1000000)