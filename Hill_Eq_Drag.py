import numpy as np
from scipy.special import erf
import ussa1976
import Constants as ct
import matplotlib.pyplot as plt

# This file is for initial analysis on how Drag can be used to increase the relative position of satellite.
# This model is based on Hill's Equation with constant drag force in hill frame. 
# Source - https://ensatellite.com/hills-equations/

def hills_drag_rel_state(t, n, f_B, f_V, f_N):
    """
    Compute relative state vector in hill frame under constant acceleration.
    Note - Reference orbit is assumed circular in derivation of Hill's Equation
    Note - The initial conditions are assumed to be zero, that means initially reference and actual state vector are same.

    Parameters:
        t   : time (scalar or array) (seconds)
        n   : mean motion of reference orbit (Circular Orbit)
        f_B : radial acceleration in hill frame (m/s^2)
        f_V : along-track acceleration in hill frame (m/s^2)
        f_N : cross-track acceleration in hill frame (m/s^2)

    Returns:
        relative state vector: [x, y, z, xdot, ydot, zdot] in hill frame (in meters)
    """
    # Matrix
    Phi = np.array([[(1 - np.cos(n*t))/n**2,           (2*t)/n - (2*np.sin(n*t))/n**2,                                 0],
                    [(2*np.sin(n*t))/n**2 - (2*t)/n,   (4*(1 - np.cos(n*t)))/n**2 - (3*t**2)/2,                        0],
                    [0,                                0,                                         (1 - np.cos(n*t))/n**2],
                    [np.sin(n*t)/n,                    (2*(1 - np.cos(n*t)))/n,                                        0],
                    [(2*(np.cos(n*t) - 1))/n,          (4*np.sin(n*t))/n - 3*t,                                        0],
                    [0,                                0,                                                  np.sin(n*t)/n]])

    # Force vector
    f_vec = np.array([f_B, f_V, f_N])

    # Returns state vector
    return Phi @ f_vec


def drag_acceleration_vnb(r_eci, v_eci, alpha, beta, ob, cenb):
    """
    Compute aerodynamic drag force in VNB frame

    Inputs:
        r_eci  : position vector (km) in Inertial frame
        v_eci  : velocity vector (km/s) in Inertial frame
        alpha  : angle of attack (rad)
        beta   : sideslip angle (rad)
        ob     : orbiting body
        cenb   : central body
    
    Returns:
        a_vnb  : acceleration vector (km/s^2) in VNB frame
    """
    # Geometry
    lx = ob["lx"]
    ly = ob["ly"]
    lz = ob["lz"]
    A_ref = ly * lz
    sigma_N = ob["sigma_N"]
    sigma_T = ob["sigma_T"]

    # Altitude
    r_norm = np.linalg.norm(r_eci)
    h = r_norm - cenb["radius"]  # km
    # Atmosphere
    ds = ussa1976.compute(z=np.array([h*1000]), variables=["p", "rho", "t"])
    rho = ds["rho"].values[0]
    T_inf = ds["t"].values[0]
    p = ds["p"].values[0]

    # Molecular speed ratio
    v_mag = np.linalg.norm(v_eci) * 1000  # m/s
    s = v_mag / np.sqrt(2 * (p / rho))

    # Coefficient
    ca, sa = np.cos(alpha), np.sin(alpha)
    cb, sb = np.cos(beta),  np.sin(beta)
    Tr = ob["T_wall"] / T_inf
    u = ca * cb
    y = sa * cb
    # C_A
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

    # Force vector in body frame
    q = 0.5 * rho * v_mag**2
    F_body = np.array([-q*A_ref*C_A, -q*A_ref*C_S, -q*A_ref*C_N])
    # Rotation matrix (Body → VNB)
    R_body2vnb = np.array([[ca*cb, -sb, cb*sa],
                           [ca*sb,  cb, sa*sb],
                           [-sa,     0,    ca]])

    return (R_body2vnb @ F_body)/(1000*ob["mass"])


# Simulation Parameters
t_final = 2*(24)*3600  # seconds
dt = 10 # seconds
t_array = np.arange(0, t_final+1, dt)

# Reference orbit (circular)
r0 = np.array([7000, 0, 0])  # km
r_cir = np.linalg.norm(r0) # magnitude in km
v_cir = np.sqrt(ct.Earth["mu"]/(np.linalg.norm(r0))) # magnitude of orbital velocity in km/s
v0 = np.array([0, v_cir, 0])   # km/s
n = np.sqrt(ct.Earth["mu"] / (np.linalg.norm(r0))**3) # Mean motion

# Orientation angles
betas = [0, 5, 10, 15, 20, 25, 30, 35]  # degrees
alpha = 0

all_states = {}
for beta_deg in betas:
    beta = np.deg2rad(beta_deg)
    states = []
    for t in t_array:
        a_drag = drag_acceleration_vnb(r0, v0, alpha, beta, ct.CubeSat_TerpRaptor, ct.Earth)
        rel_state = hills_drag_rel_state(t, n, a_drag[2]*1000, a_drag[0]*1000, a_drag[1]*1000)
        x_rel, y_rel, z_rel = rel_state[:3]/1000

        theta = n * t
        x_ref = r_cir * np.cos(theta)
        y_ref = r_cir * np.sin(theta)
        z_ref = 0
        states.append([x_ref + x_rel, y_ref + y_rel, z_ref + z_rel])
    all_states[beta_deg] = np.array(states)

# Plot
fig, axes = plt.subplots(2, 2, figsize=(12, 8))
axes = axes.flatten()

labels = ['Δx (m)', 'Δy (m)', 'Δz (m)', '|Δr| (m)']
time_hr = t_array / 3600

for beta in betas:
    diff = all_states[beta] - all_states[0]
    mag = np.linalg.norm(diff, axis=1)
    data = [diff[:, 0], diff[:, 1], diff[:, 2], mag]

    for ax, d, label in zip(axes, data, labels):
        ax.plot(time_hr, d*1000, label=f'β = {beta}°')
        ax.set_ylabel(label)
        ax.grid(True)

for ax in axes:
    ax.set_xlabel('Time (hours)')
    ax.legend()

fig.suptitle('Relative Position Differences for Varying β', fontsize=14)
plt.tight_layout()
plt.show()