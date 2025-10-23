import numpy as np
from OrbitPropagator import OrbitPropagator, null_perturbations
import Constants as ct
import time

# This scenario is from a research paper.
# Source - https://digitalcommons.usu.edu/cgi/viewcontent.cgi?article=6317&context=smallsat

if __name__ == '__main__':
    start_time = time.time()
    print("Orbit propagator is running", flush=True)

    # Initial conditions of satellite
    r0_S, v0_S = OrbitPropagator.Kep2State(6880, 0, 51.6, 0, 0, 0, cb=ct.Earth)

    # Initial conditions of debri
    r0_D, v0_D = OrbitPropagator.Kep2State(6880, 0, 97.8, 350, 187, 180.00415, cb=ct.Earth)

    # Time span and time step
    tspan = (100) * (60*60) # hours
    dt = 10 # seconds

    # Propagate orbit
    perturbations = null_perturbations()
    perturbations['J2'] = False
    perturbations['J3'] = False
    perturbations['Drag'] = False
    ob_S = OrbitPropagator(r0_S,v0_S,tspan,dt,cb=ct.Earth,ob=ct.Satellite,perturbations=perturbations)
    ob_D = OrbitPropagator(r0_D,v0_D,tspan,dt,cb=ct.Earth,ob=ct.Debri,perturbations=perturbations)
    end_time = time.time()
    print('Time Taken to execute propagator code: '+ str(end_time-start_time) +' seconds')

    # Collision Probability 
    Pc_array = np.zeros(ob_S.rs.shape[0])
    for k in range(ob_S.rs.shape[0]):
        Pc_array[k] = OrbitPropagator.Pc2D(ob_S.rs[k], ob_S.vs[k], ct.Satellite['Covarience_Matrix'], 
            ob_D.rs[k], ob_D.vs[k], ct.Debri['Covarience_Matrix'], 
            (ct.Satellite['HBR_Radius'] + ct.Debri['HBR_Radius']))
        
    r_array = np.zeros(((ob_S.rs.shape[0]), 3))
    for k in range(ob_S.rs.shape[0]):
        r_array[k] = (ob_S.rs[k] - ob_D.rs[k])

    distances = np.linalg.norm(r_array, axis=1)
    
    print(r_array)
    print(distances)
    print(np.max(Pc_array))
    print(np.min(distances))
    # Plots
    # OrbitPropagator.PlotOrbits([ob_S.rs, ob_D.rs], labels=['Satellite', 'Debri'], cb=ct.Earth, title='Synthetic Collision Situation')


