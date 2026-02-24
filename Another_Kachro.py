import numpy as np
from OrbitPropagator import OrbitPropagator, null_perturbations
import Constants as ct
import time

# Use this file for copying the syntax.
if __name__ == '__main__':
    start_time = time.time()
    print("Orbit propagator is running", flush=True)

    # Initial conditions of satellite
    r0_S = OrbitPropagator.State2Kep(np.array([[3634.1,5926,1206.6]]), np.array([[-6.9049,4.3136,2.6163]]), cb=ct.Earth)

    print(r0_S)