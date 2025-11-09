import numpy as np
from OrbitPropagator import OrbitPropagator, null_perturbations
import Constants as ct
import time

if __name__ == '__main__':
    print(OrbitPropagator.State2Kep(np.array([[-5107.28759593, -2725.03856697, -3697.97987789]]), np.array([[-5.08860732, 3.65387192, 4.33905895]]), cb=ct.Earth))