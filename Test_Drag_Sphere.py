import numpy as np
from scipy.special import erf
import ussa1976
import Constants as ct

# Source - https://arc.aiaa.org/doi/suppl/10.2514/1.A33606
# Sample Inputs
h = 500   # altitude in km
v = np.sqrt(ct.Earth["mu"]/(h+ct.Earth["radius"]))    # orbital velocity in km/s
T_wall = 300.0
sigma_N = 1.0 # Good asumption
sigma_T = 1.0 # Good asumption

# Atmosphere at given altitude
ds = ussa1976.compute(z=np.array([h*1000]), variables=["p", "rho", "t"])
rho = ds["rho"].values # kg/m^3
T_inf = ds["t"].values # K
p = ds["p"].values # N/m^2

# Molecular Speed Ratio
s = (v*1000) / np.sqrt(2 * (p/rho))

# C_D formula
term1 = (2 - sigma_N + sigma_T) * (0.5 +
    0.5 * erf(s) * (1 + 1/s**2 - 1/(4*s**4)) +
    (1 + 2*s**2) * np.exp(-s**2) / (4*s**3*np.sqrt(np.pi)))
term2 = (sigma_N / (3*s)) * np.sqrt(np.pi * (T_wall/T_inf)) * (1 + erf(s))
term3 = (2 - sigma_N)/(2*s**2) 
term4 = (sigma_N/(6*s**4)) * (1 + (2*s**2 - 1)*np.exp(-s**2)) * np.sqrt(T_wall/T_inf)
CD = term1 + term2 + term3 + term4

# Drag force
A_projected = np.pi * (1)**2 # Unit circle
F = 0.5 * rho * (v*1000)**2 * CD * A_projected

print(f"rho     = {rho[0]:.8e} kg/m3")
print(f"T_inf   = {T_inf[0]:.8f} K")
print(f"s       = {s[0]:.8f}")
print(f"CD      = {CD[0]:.8f}")
print(f"F_drag  = {F[0]:.8e} N")