import numpy as np
from scipy.special import erf
import ussa1976
import Constants as ct
import matplotlib.pyplot as plt
from scipy.spatial.transform import Rotation as R

# Source - https://arc.aiaa.org/doi/suppl/10.2514/1.A33606
# Sample Inputs
h = 500   # altitude in km
v_cir = np.sqrt(ct.Earth["mu"]/(h+ct.Earth["radius"]))    # orbital velocity in km/s
T_wall = 300.0
sigma_N = 1.0 # Good asumption
sigma_T = 1.0 # Good asumption
L = 1 # Cube Side Length in m
lx, ly, lz = L, L, L
A_ref = ly * lz   # m^2

# Orientation angles
alpha = 0 * (np.pi/180) # radians
beta  = 0 * (np.pi/180) # radians

# Atmosphere at given altitude
ds = ussa1976.compute(z=np.array([h*1000]), variables=["p", "rho", "t"])
rho = ds["rho"].values # kg/m^3
T_inf = ds["t"].values # K
p = ds["p"].values # N/m^2

# Molecular Speed Ratio
s = (v_cir*1000) / np.sqrt(2 * (p/rho))

# Trig shortcuts
ca, sa = np.cos(alpha), np.sin(alpha)
cb, sb = np.cos(beta),  np.sin(beta)
Tr = T_wall / T_inf # Temperature ratio

# Force Coefficients Formulas
# Note - # Signum Function is np.sign function
# C_A
u = ca * cb
term1 = ( ((2 - sigma_N)/(s * np.sqrt(np.pi)) * u) + (np.sign(u) * (sigma_N/(2*s**2)) * np.sqrt(Tr)) ) * (np.exp(-s**2 * u**2))
term2 = (2 - sigma_N) * (u**2 + (1/(2*s**2))) * (np.sign(u) + erf(s*u))
term3 = (sigma_N/(2*s) * u * np.sqrt(np.pi*Tr)) * (1 + np.sign(u)*erf(s*u))
term4 = (sigma_T*u*(lx/ly)) * ( (1/(s*np.sqrt(np.pi)) * np.exp(-s**2 * sb**2)) + (sb*(np.sign(sb)+erf(s*sb))) )
term5 = (sigma_T*u*(lx/lz)) * ( (1/(s*np.sqrt(np.pi)) * np.exp(-s**2 * sa**2 * cb**2)) + (sa*cb*(erf(s*sa*cb)+np.sign(sa*cb))) )
C_A = term1 + term2 + term3 + term4 + term5

# C_S
term6 = (lx/ly) * ( ((2 - sigma_N)/(s * np.sqrt(np.pi)) * sb) + (np.sign(sb) * (sigma_N/(2*s**2)) * np.sqrt(Tr)) ) * (np.exp(-s**2 * sb**2))
term7 = (lx/ly) * (2 - sigma_N) * (sb**2 + (1/(2*s**2))) * (np.sign(sb) + erf(s*sb))
term8 = (lx/ly) * (sigma_N/(2*s) * sb * np.sqrt(np.pi*Tr)) * (1 + np.sign(sb)*erf(s*sb))
term9 = (sigma_T*sb) * ( (1/(s*np.sqrt(np.pi)) * np.exp(-s**2 * u**2)) + (u*(erf(s*u)+np.sign(u))) )
term10 = (sigma_T*sb*(lx/lz)) * ( (1/(s*np.sqrt(np.pi)) * np.exp(-s**2 * sa**2 * cb**2)) + (sa*cb*(erf(s*sa*cb)+np.sign(sa*cb))) )
C_S = term6 + term7 + term8 + term9 + term10

# C_N
y = sa * cb
term11 = (lx/lz) * ( ((2 - sigma_N)/(s * np.sqrt(np.pi)) * y) + (np.sign(y) * (sigma_N/(2*s**2)) * np.sqrt(Tr)) ) * (np.exp(-s**2 * y**2))
term12 = (lx/lz) * (2 - sigma_N) * (y**2 + (1/(2*s**2))) * (np.sign(y) + erf(s*y))
term13 = (lx/lz) * (sigma_N/(2*s) * y * np.sqrt(np.pi*Tr)) * (1 + np.sign(y)*erf(s*y))
term14 = (sigma_T*y) * ( (1/(s*np.sqrt(np.pi)) * np.exp(-s**2 * u**2)) + (u*(erf(s*u)+np.sign(u))) )
term15 = (sigma_T*y*(lx/ly)) * ( (1/(s*np.sqrt(np.pi)) * np.exp(-s**2 * sb**2)) + (sb*(erf(s*sb)+np.sign(sb))) )
C_N = term11 + term12 + term13 + term14 + term15

# Inertial frame
r = np.array([h*1000,0,0]) # m
v = np.array([0,v_cir*1000,0]) # m/s
# VNB frame
v_hat = v/np.linalg.norm(v)
n_hat = np.cross(r,v)/np.linalg.norm(np.cross(r,v))
b_hat = np.cross(v_hat,n_hat)

# Force vector in body frame
q = 0.5 * rho * (v_cir*1000)**2   # dynamic pressure (N/m²)
F_body_axial = -q * A_ref * C_A
F_body_side = -q * A_ref * C_S
F_body_normal = -q * A_ref * C_N
F_body = np.array([F_body_axial, F_body_normal, F_body_side]).flatten()

# Rotation matrix (VNB → Body)
r_vnb2body = R.from_euler('zy', [alpha, beta])
R_vnb2body = r_vnb2body.as_matrix()
# Rotation matrix (Body → VNB)
R_body2vnb = R_vnb2body.T

# Force vector in VNB frame
F_vnb = R_body2vnb @ F_body
print(F_vnb)
print(F_body)
print(C_A, C_N, C_S)
print(np.linalg.norm(F_body))
print(np.linalg.norm(F_vnb))