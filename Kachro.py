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
    r0_D = [-6653.02579823, 1071.92988144, 1329.48660569]
    v0_D = [1.89493254, 4.58390171, 5.78944253]

    # Time span and time step
    tspan = (12) * (60*60) # hours
    dt = 10 # seconds

    # Propagate orbit
    perturbations = null_perturbations()
    perturbations['J2'] = True
    perturbations['J3'] = False
    perturbations['Drag'] = False
    
    ob_S = OrbitPropagator(r0_S,v0_S,tspan,dt,cb=ct.Earth,ob=ct.Satellite,perturbations=perturbations)
    ob_D = OrbitPropagator(r0_D,v0_D,tspan,dt,cb=ct.Earth,ob=ct.Debri,perturbations=perturbations)
    end_time = time.time()
    print('Time Taken to execute propagator code: '+ str(end_time-start_time) +' seconds')

    print(r0_S)
    print(ob_S.rs[-1])
    print(ob_S.vs[-1])
    print(ob_D.rs[-1])
    print(ob_D.rs[-1] - r0_S)

    # OrbitPropagator.PlotOrbits([ob_S.rs], labels=['Satellite'], cb=ct.Earth, title='Synthetic Collision Situation')
    # OrbitPropagator.PlotOrbits([ob_S.rs, ob_D.rs], labels=['Satellite', 'Debri'], cb=ct.Earth, title='Synthetic Collision Situation')
    # OrbitPropagator.PlotKepOrbits(ob_S.rs,ob_S.vs,ob_S.ts,cb=ct.Earth,title="Evolution of Satellite's Keplerian Elements")

