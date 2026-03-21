import numpy as np
from scipy.special import erf
import ussa1976
import Constants as ct
import matplotlib.pyplot as plt
from scipy.spatial.transform import Rotation as R

# Sample Inputs
h = 500   # altitude in km
v_cir = np.sqrt(ct.Earth["mu"]/(h+ct.Earth["radius"]))    # orbital velocity in km/s
T_wall = 300.0
sigma_N = 1.0
sigma_T = 1.0
L = 1
lx, ly, lz = L, L, L
A_ref = ly * lz

# Atmosphere
ds = ussa1976.compute(z=np.array([h*1000]), variables=["p", "rho", "t"])
rho   = ds["rho"].values
T_inf = ds["t"].values
p     = ds["p"].values
s     = (v_cir*1000) / np.sqrt(2 * (p/rho))
Tr    = T_wall / T_inf
q     = 0.5 * rho * (v_cir*1000)**2

# Meshgrid
alpha_vals = np.linspace(-np.pi/2, np.pi/2, 500)
beta_vals  = np.linspace(-np.pi/2, np.pi/2, 500)
ALPHA, BETA = np.meshgrid(alpha_vals, beta_vals)

ca, sa = np.cos(ALPHA), np.sin(ALPHA)
cb, sb = np.cos(BETA),  np.sin(BETA)

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


# Force in body frame
F_body_a = -q * A_ref * C_A   # axial
F_body_s = -q * A_ref * C_S   # side
F_body_n = -q * A_ref * C_N   # normal

# Rotate body → VNB for each (alpha, beta) pair
F_vnb_v = np.zeros_like(ALPHA)
F_vnb_n = np.zeros_like(ALPHA)
F_vnb_b = np.zeros_like(ALPHA)

for i in range(ALPHA.shape[0]):
    for j in range(ALPHA.shape[1]):
        a = ALPHA[i, j]
        b = BETA[i, j]

        r_vnb2body = R.from_euler('zy', [a, b])
        R_body2vnb = r_vnb2body.as_matrix().T

        F_body_vec = np.array([F_body_a[i, j], F_body_n[i, j], F_body_s[i, j]])
        F_vnb_vec = R_body2vnb @ F_body_vec

        F_vnb_v[i, j] = F_vnb_vec[0]
        F_vnb_n[i, j] = F_vnb_vec[1]
        F_vnb_b[i, j] = F_vnb_vec[2]

# Plotting
fig, axes = plt.subplots(1, 3, figsize=(18, 6))

data = [
    (F_vnb_v, "F_V (along velocity)",  "F$_V$ (N)"),
    (F_vnb_n, "F_N (orbit normal)",     "F$_N$ (N)"),
    (F_vnb_b, "F_B (binormal)",         "F$_B$ (N)"),
]

for ax, (Z, title, label) in zip(axes, data):
    cp = ax.contourf(np.degrees(ALPHA), np.degrees(BETA), Z, levels=40, cmap="plasma")
    plt.colorbar(cp, ax=ax, label=label)
    ax.set_xlabel("α (deg)")
    ax.set_ylabel("β (deg)")
    ax.set_title(title)

plt.suptitle(f"Forces in VNB Frame — h={h}km, L={L}m", fontsize=13)
plt.tight_layout()
plt.show()