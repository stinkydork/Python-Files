import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erf
import ussa1976
import Constants as ct

# Inputs
ob = ct.CubeSat_TerpRaptor
r_cir = 7000  # km
v_cir = np.sqrt(ct.Earth["mu"]/r_cir)    # orbital velocity in km/s
h = r_cir - ct.Earth["radius"] # Altitude in km
T_wall = ob["T_wall"]
sigma_N = ob["sigma_N"]
sigma_T = ob["sigma_T"]
lx, ly, lz = ob["lx"], ob["ly"], ob["lz"]
A_ref = ly * lz   # m^2
n = np.sqrt(ct.Earth["mu"] / (r_cir)**3) # Mean motion

# Atmosphere at given altitude
ds = ussa1976.compute(z=np.array([h*1000]), variables=["p", "rho", "t"])
rho = ds["rho"].values # kg/m^3
T_inf = ds["t"].values # K
p = ds["p"].values # N/m^2
Tr = T_wall / T_inf # Temperature ratio

# Molecular Speed Ratio
s = (v_cir*1000) / np.sqrt(2 * (p/rho))

# Orientation angle values
alpha = 0
beta_vals = np.linspace(0, 80*(np.pi/180), 500)

# Force calculations
ca, cb, sa, sb = np.cos(alpha), np.cos(beta_vals), np.sin(alpha), np.sin(beta_vals)
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
# Forces in body frame
q = 0.5 * rho * (v_cir*1000)**2
F_body_axial  = -q * A_ref * C_A
F_body_side   = -q * A_ref * C_S
F_body_normal = -q * A_ref * C_N
# Forces in VNB frame
F_v = (ca*cb)*F_body_axial + (-sb)*F_body_side + (cb*sa)*F_body_normal
F_n = (ca*sb)*F_body_axial + (cb)*F_body_side + (sa*sb)*F_body_normal
F_b = (-sa)*F_body_axial + (0)*F_body_side + (ca)*F_body_normal

# Plot
plt.figure(figsize=(8, 5))
for day in range(1, 8):
    t = day * 86400 # seconds
    # Note - In below, subtrated term is of force when s/c is at default orientation.
    y = ((4*(1 - np.cos(n*t))/n**2 - 3*t**2/2) * ((F_v-F_v[0]) / ob["mass"]))
    plt.plot(beta_vals*(180/np.pi), y, label=f"t = {day} days")
plt.axhline(y=-1000, color='black', linestyle='--', linewidth=1.5, label='|Δy| = 1000 m')

plt.xlabel("β (deg)")
plt.ylabel("Achievable Δy (m)")
plt.title("Required β for given Δy | Terp Raptor")
plt.legend()
plt.grid(True)
plt.show()