import numpy as np
from OrbitPropagator import OrbitPropagator, null_perturbations
import Constants as ct
import time

# Use this file for copying the syntax.
if __name__ == '__main__':
    start_time = time.time()
    print("Orbit propagator is running", flush=True)

    # Initial conditions of satellite
    r0_S, v0_S = OrbitPropagator.Kep2State(11438, 0.3, 0, 0, 0, 300, cb=ct.Earth)

    # Initial conditions of debri
    r0_D, v0_D = OrbitPropagator.Kep2State(6880, 0, 97.8, 350, 187, 180.00415, cb=ct.Earth)

    # Time span and time step
    tspan = 1484.06388559# (12) * (60*60) # hours
    dt = 1 # seconds

    # Propagate orbit
    perturbations = null_perturbations()
    perturbations['J2'] = False
    perturbations['J3'] = False
    perturbations['Drag'] = False
    perturbations['SRP'] = False
    perturbations['External Thrust'] = False
    
    ob_S = OrbitPropagator(r0_S,v0_S,0,tspan,dt,cb=ct.Earth,ob=ct.CubeSat_TerpRaptor,perturbations=perturbations)
    # ob_D = OrbitPropagator(r0_D,v0_D,tspan,dt,cb=ct.Earth,ob=ct.Debri,perturbations=perturbations)
    end_time = time.time()
    print('Time Taken to execute propagator code: '+ str(end_time-start_time) +' seconds')

    # OrbitPropagator.PlotOrbits([ob_S.rs], labels=['Satellite'], cb=ct.Earth, title='3D plot of a propagated LEO trajectory around Earth')
    # OrbitPropagator.PlotOrbits([ob_S.rs, ob_D.rs], labels=['Satellite', 'Debri'], cb=ct.Earth, title='Synthetic Collision Situation')
    OrbitPropagator.PlotKepOrbits(ob_S.rs,ob_S.vs,ob_S.ts,cb=ct.Earth,title="Evolution of Orbiting Body's Keplerian Elements")
    # OrbitPropagator.PlotEnergy(ob_S.rs,ob_S.vs,ob_S.ts,cb=ct.Earth,title="Specific Energy vs. Time")