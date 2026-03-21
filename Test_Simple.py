import numpy as np
from OrbitPropagator import OrbitPropagator, all_perturbations
import Constants as ct
import time

# Use this file for copying the syntax.
if __name__ == '__main__':
    start_time = time.time()
    print("Orbit propagator is running", flush=True)

    # Initial conditions of satellite
    r0_S, v0_S = OrbitPropagator.Kep2State(7000, 0, 0, 0, 0, 0, cb=ct.Earth)

    # Initial conditions of debri
    r0_D, v0_D = OrbitPropagator.Kep2State(6880, 0, 97.8, 350, 187, 180.00415, cb=ct.Earth)

    # Time span and time step
    tspan = (12) * (60*60) # hours
    dt = 10 # seconds
    perturbations = all_perturbations()
    perturbations['SRP'] = False

    # Propagate orbit
    ob_S = OrbitPropagator(r0_S,v0_S,tspan,dt,cenb=ct.Earth,ob=ct.CubeSat_1U,pointing_mode='nadir',perturbations=perturbations)
    # ob_D = OrbitPropagator(r0_D,v0_D,tspan,dt,cenb=ct.Earth,ob=ct.Debri,perturbations=perturbations)
    end_time = time.time()
    print('Time Taken to execute propagator code: '+ str(end_time-start_time) +' seconds')

    # OrbitPropagator.PlotOrbits([ob_S.rs], labels=['Satellite'], cb=ct.Earth, title='3D plot of a propagated LEO trajectory around Earth')
    # OrbitPropagator.PlotOrbits([ob_S.rs, ob_D.rs], labels=['Satellite', 'Debri'], cb=ct.Earth, title='Synthetic Collision Situation')
    OrbitPropagator.PlotKepOrbits(ob_S.rs,ob_S.vs,ob_S.ts,cb=ct.Earth,title="Evolution of Orbiting Body's Keplerian Elements")
    # OrbitPropagator.PlotEnergy(ob_S.rs,ob_S.vs,ob_S.ts,cb=ct.Earth,title="Specific Energy vs. Time")