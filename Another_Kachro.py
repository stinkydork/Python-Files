import numpy as np
from OrbitPropagator import OrbitPropagator, null_perturbations
import Constants as ct
import time

if __name__ == '__main__':
    start_time = time.time()
    print("Orbit propagator is running", flush=True)

    # Initial conditions of satellite
    r0_S, v0_S = OrbitPropagator.Kep2State(6880, 0, 51.6, 0, 0, 0, cb=ct.Earth)

    # Initial conditions of debri
    r0_D, v0_D = OrbitPropagator.Kep2State(6880, 0, 97.8, 350, 187, 180.00415, cb=ct.Earth)

    # Time span and time step
    tspan = (11) * (60*60) # hours
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

    # OrbitPropagator.PlotOrbits([ob_S.rs], labels=['Satellite'], cb=ct.Earth, title='Synthetic Collision Situation')
    OrbitPropagator.PlotOrbits([ob_S.rs, ob_D.rs], labels=['Satellite', 'Debri'], cb=ct.Earth, title='Synthetic Collision Situation')
    # OrbitPropagator.PlotKepOrbits(ob_S.rs,ob_S.vs,ob_S.ts,cb=ct.Earth,title="Evolution of Satellite's Keplerian Elements")
    
    print(ob_S.rs[-1])
    print(ob_D.rs[-1])
    Pc_array = OrbitPropagator.Pc2D(ob_S.rs[-1], ob_S.vs[-1], ct.Satellite['Covarience_Matrix'], 
            ob_D.rs[-1], ob_D.vs[-1], ct.Debri['Covarience_Matrix'], 
            (ct.Satellite['HBR_Radius'] + ct.Debri['HBR_Radius']))
    print(Pc_array)