import numpy as np
from OrbitPropagator import OrbitPropagator, null_perturbations
import Constants as ct
import time

# This scenario is from a research paper.
# Source - https://digitalcommons.usu.edu/cgi/viewcontent.cgi?article=6317&context=smallsat

if __name__ == '__main__':
    
    Pc_array = OrbitPropagator.Pc2D(np.array([7500,0,0]), np.array([0,10,0]), ct.Satellite['Covarience_Matrix'], 
            np.array([7500.5,0,0]), np.array([0,10.5,0]), ct.Debri['Covarience_Matrix'], 
            ct.Satellite['HBR_Radius'] + ct.Debri['HBR_Radius'])


    print(Pc_array)