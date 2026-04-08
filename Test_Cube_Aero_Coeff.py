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


# Atmosphere at given altitude
ds = ussa1976.compute(z=np.array([h*1000]), variables=["p", "rho", "t"])
rho = ds["rho"].values # kg/m^3
T_inf = ds["t"].values # K
p = ds["p"].values # N/m^2
Tr = T_wall / T_inf # Temperature ratio

# Molecular Speed Ratio
s = (v_cir*1000) / np.sqrt(2 * (p/rho))

# Meshgrid over orientation angles
alpha_vals = np.linspace(-20 * np.pi/180, 20 * np.pi/180, 200)
beta_vals  = np.linspace(-20 * np.pi/180, 20 * np.pi/180, 200)
ALPHA, BETA = np.meshgrid(alpha_vals, beta_vals)
ca, sa = np.cos(ALPHA), np.sin(ALPHA)
cb, sb = np.cos(BETA),  np.sin(BETA)

# Force Coefficients Formulas
# Note - # Signum Function is np.sign function
u = ca * cb
y = sa * cb
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

# Plot
fig, axes = plt.subplots(1, 3, figsize=(18, 5))
alpha_deg = np.degrees(alpha_vals)
beta_deg  = np.degrees(beta_vals)

for ax, (Z, title, label) in zip(axes, [(C_A, "C$_A$", ""), (C_S, "C$_S$", ""), (C_N, "C$_N$", "")]):
    cp = ax.contourf(alpha_deg, beta_deg, Z, levels=40, cmap="YlOrRd")
    plt.colorbar(cp, ax=ax, label=label)
    ax.set_title(f"{title}")
    ax.set_xlabel("α (deg)")

axes[0].set_ylabel("β (deg)")
plt.suptitle(f"Aerodynamic Coefficients — h={h}km, L={L}m, v={v_cir*1000:.0f}m/s", fontsize=13)
plt.tight_layout()
plt.show()