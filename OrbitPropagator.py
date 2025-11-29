# Main Orbit Propagator File
import numpy as np
import matplotlib.pyplot as plt 
from scipy.integrate import ode
import Constants as ct
import ussa1976
import spiceypy as spice
import os
plt.style.use('dark_background')

# For selecting which type of perturbation is needed.
def null_perturbations():
    return {
        'J2' : False,
        'J3' : False,
        'Drag' : False,
        'SRP' : False,
        'External Thrust' : False
    }


# Kernels for NASA Spice files.
KERNEL_DIR = r"C:\Engineering Program Files\spice_kernels" # This is where the kernels are located.
LSK  = fr"{KERNEL_DIR}\lsk\naif0012.tls"
SPK1 = fr"{KERNEL_DIR}\spk\de441_part-1.bsp"
SPK2 = fr"{KERNEL_DIR}\spk\de441_part-2.bsp"
def load_kernels():
    # clear any existing kernels (optional)
    try:
        spice.kclear()
    except Exception:
        pass

    # Furnsh (load) kernels
    for k in (LSK, SPK1, SPK2):
        if not os.path.exists(k):
            raise FileNotFoundError(f"Kernel not found: {k}")
        spice.furnsh(k)


# Standard Unit of Distance is Kilometers.
# Starting time reference point is (2000-01-01 12:00:00 UTC).
class OrbitPropagator:
    # Integrater Function
    # Using DOP853 method.
    # Source - http://www.youtube.com/@alfonsogonzalez-astrodynam2207
    def __init__(self,r0,v0,tspan,dt,cb,ob,perturbations=null_perturbations()):
        load_kernels()
        self.r0 = r0
        self.v0 = v0
        self.tspan = tspan
        self.dt = dt
        self.cb = cb
        self.ob = ob
        self.altitude = np.linspace(0,1000000,1001)
        self.rho_ds = ussa1976.compute(z=self.altitude, variables=["rho"])
        self.rho_table = self.rho_ds["rho"].values.squeeze()

        # Total number of steps
        self.n_steps=int(np.ceil(self.tspan/self.dt))

        # Initializing arrays for state vector and time
        self.ys=np.zeros((self.n_steps,6))
        self.ts=np.zeros((self.n_steps,1))

        # Initial conditions
        self.y0 = np.hstack((self.r0, self.v0))
        self.ts[0] = 0
        self.ys[0] = np.array(self.y0)
        self.step=1

        # Initiate solver
        self.solver=ode(self.diffy_q)
        self.solver.set_integrator('dop853', rtol=1e-10, atol=1e-12)
        self.solver.set_initial_value(self.y0,0)
    
        # Pertubations types
        self.perturbations = perturbations

        # Propagate Orbit
        while self.solver.successful() and self.step < self.n_steps:
            self.solver.integrate(self.solver.t + self.dt)
            self.ts[self.step] = self.solver.t
            self.ys[self.step] = self.solver.y
            self.step += 1
            print(self.solver.y[:3], self.solver.t) # Use only for debugging
        self.rs = self.ys[:,:3]
        self.vs = self.ys[:,3:]


    # Force Model Function
    # Two Body acceleration is always on.
    # Source 1 - http://www.youtube.com/@alfonsogonzalez-astrodynam2207
    # Source 2 - Fundamentals of Astrodynamics (Textbook)
    # Source 3 - https://farside.ph.utexas.edu/teaching/celestial/Celestial/node94.html
    # Source 4 - https://ussa1976.readthedocs.io/en/latest/index.html
    def diffy_q(self,t,y):
        # Making proper arrays
        self.rx, self.ry, self.rz, self.vx, self.vy, self.vz = y
        self.r = np.array([self.rx, self.ry, self.rz])
        self.v = np.array([self.vx, self.vy, self.vz])

        # Norm of the position vector
        self.norm_r = np.linalg.norm(self.r)

        # Two-body acceleration
        self.a = -self.r * self.cb['mu'] / self.norm_r**3

        # J2 perturbation
        if self.perturbations['J2']:
            self.j2x = self.r[0] * ((5 * self.r[2]**2 / self.norm_r**2) - 1)
            self.j2y = self.r[1] * ((5 * self.r[2]**2 / self.norm_r**2) - 1)
            self.j2z = self.r[2] * ((5 * self.r[2]**2 / self.norm_r**2) - 3)
            self.a_j2 = (1.5 * self.cb['J2'] * self.cb['mu'] * self.cb['radius']**2 / self.norm_r**5) * np.array([self.j2x, self.j2y, self.j2z])
            self.a += self.a_j2

        # J3 perturbation
        if self.perturbations['J3']:
            self.j3x = (5 * self.r[0]) * ((7 * self.r[2]**3 / self.norm_r**3) - (3 * self.r[2] / self.norm_r))
            self.j3y = (5 * self.r[1]) * ((7 * self.r[2]**3 / self.norm_r**3) - (3 * self.r[2] / self.norm_r))
            self.j3z = ((35 * self.r[2]**4 / self.norm_r**4) - (30 * self.r[2]**2 / self.norm_r**2) + 3)
            self.a_j3 = (0.5 * self.cb['J3'] * self.cb['mu'] * self.cb['radius']**3 / self.norm_r**6) * np.array([self.j3x, self.j3y, self.j3z])
            self.a += self.a_j3

        # Drag perturbation
        if self.perturbations['Drag']:
            self.z = (self.norm_r - self.cb['radius']) * 1000.0
            if self.z < 1000000:
                self.rho = self.rho_table[round(self.z/1000)]
                self.vrel = self.v - np.cross(self.cb['at_rot_vec'], self.r)
                self.a_drag = -((0.5 * self.rho * self.ob['Cd'] * self.ob['Drag_Area']) / self.ob['mass']) * np.linalg.norm(self.vrel) * self.vrel
                self.a += self.a_drag
            else:
                self.a = self.a

        # Solar Radiation Pressure
        if self.perturbations['SRP']:
            self.sun_pos, self.sun_vel = OrbitPropagator.sun_position_icrf(self.solver.t)
            self.right_ascension_sun = np.atan2(self.sun_pos[1], self.sun_pos[0])
            if self.right_ascension_sun < 0:
                self.right_ascension_sun += 2*np.pi
            self.a_srp_x = ((-4.5e-5) * ( (self.ob['SRP_Area']*1e+10)/(self.ob['mass']*1000) )) * np.cos(self.right_ascension_sun)
            self.a_srp_y = ((-4.5e-5) * ( (self.ob['SRP_Area']*1e+10)/(self.ob['mass']*1000) )) * np.sin(self.right_ascension_sun) * np.cos(np.radians(23.4349))
            self.a_srp_z = ((-4.5e-5) * ( (self.ob['SRP_Area']*1e+10)/(self.ob['mass']*1000) )) * np.sin(self.right_ascension_sun) * np.sin(np.radians(23.4349))
            self.a_srp = np.array([self.a_srp_x, self.a_srp_y, self.a_srp_z]) / 100000
            self.a += self.a_srp

        # This is the custom code block for the external thrust
        if self.perturbations['External Thrust']:
            self.a = self.a

        return [self.vx, self.vy, self.vz, self.a[0], self.a[1], self.a[2]]
    

    # Keplerian-Orbital elements to State Vector Converter
    # Note - This function only take degrees as input.
    # Source - https://orbital-mechanics.space/classical-orbital-elements/orbital-elements-and-the-state-vector.html
    def Kep2State(a,e,i_deg,RAAN_deg,omega_deg,nu_deg, cb):
        i = np.radians(i_deg)
        RAAN = np.radians(RAAN_deg)
        omega = np.radians(omega_deg)
        nu = np.radians(nu_deg)
        mu = cb['mu']

        # Rotation matrix
        rotation_matrix = np.array([ [np.cos(RAAN)*np.cos(omega) - np.sin(RAAN)*np.sin(omega)*np.cos(i),
             -np.cos(RAAN)*np.sin(omega) - np.sin(RAAN)*np.cos(omega)*np.cos(i),
             np.sin(RAAN)*np.sin(i)],

            [np.sin(RAAN)*np.cos(omega) + np.cos(RAAN)*np.sin(omega)*np.cos(i),
             -np.sin(RAAN)*np.sin(omega) + np.cos(RAAN)*np.cos(omega)*np.cos(i),
             -np.cos(RAAN)*np.sin(i)],

            [np.sin(omega)*np.sin(i),
             np.cos(omega)*np.sin(i),
             np.cos(i)] ])

        # Calculating the position vector
        r = (a * (1 - e**2)) / (1 + e * np.cos(nu))
        x = r * np.cos(nu)
        y = r * np.sin(nu)
        z = 0.0
        cor_r = rotation_matrix @ np.array([x , y, z])

        # Calculating the velocity vector
        p = a * (1 - e**2)
        v_factor = np.sqrt(mu / p)
        v_x = -v_factor * np.sin(nu)
        v_y = v_factor * (e + np.cos(nu))
        v_z = 0.0
        cor_v = rotation_matrix @ np.array([v_x , v_y, v_z])

        return cor_r, cor_v
    

    # State Vector to Keplerian-Orbital elements Converter
    # Note - This function return angles in degrees.
    # Source - https://downloads.rene-schwarz.com/download/M002-Cartesian_State_Vectors_to_Keplerian_Orbit_Elements.pdf
    def State2Kep(R,V,cb):
        # Getting mu value of Central Body
        mu = cb['mu']
        # Setting up the function
        N = R.shape[0]
        kep_array = np.zeros((N, 6))
        for k in range(N):
            r = R[k]
            v = V[k]
            # Specific angular momentum
            h = np.cross(r, v)

            # Eccentricity vector
            e = (np.cross(v, h)/mu) - (r/np.linalg.norm(r))

            # True anomaly (in radians)
            if np.dot(r, v) >= 0:
                nu = np.arccos(np.dot(e, r) / (np.linalg.norm(e)*np.linalg.norm(r)))
            else:
                nu = (2*np.pi) - np.arccos(np.dot(e, r) / (np.linalg.norm(e)*np.linalg.norm(r)))

            # n vector
            n = np.array([-h[1], h[0], 0])

            # Inclination
            i = np.arccos(h[2] / np.linalg.norm(h))

            # Eccentricity
            e_mag = np.linalg.norm(e)

            # RAAN (in radians)
            if n[1] >= 0:
                RAAN = np.arccos(n[0] / np.linalg.norm(n))
            else:
                RAAN = (2*np.pi) - np.arccos(n[0] / np.linalg.norm(n))

            # Omega (in radians)
            if e[2] >= 0:
                omega = np.arccos(np.dot(n, e) / (np.linalg.norm(n)*np.linalg.norm(e)))
            else:
                omega = (2*np.pi) - np.arccos(np.dot(n, e) / (np.linalg.norm(n)*np.linalg.norm(e)))

            # Semi-major axis
            a = 1 / ( (2/np.linalg.norm(r)) - (np.linalg.norm(v)**2/mu) )

            # For each iteration
            kep_array[k] = [a, e_mag, np.degrees(i), np.degrees(RAAN), np.degrees(omega), np.degrees(nu)]

        return kep_array


    # Collision Probabilty calculator Function
    # Source - https://amostech.com/TechnicalPapers/2014/Conjunction_Assessment/DUNCAN.pdf
    def Pc2D(r1,v1,cov1,r2,v2,cov2,HBR):
        r_rel = r1 - r2
        v_rel = v1 - v2
        r_hat = r_rel / np.linalg.norm(r_rel)
        v_hat = v_rel / np.linalg.norm(v_rel)
        z_hat = np.cross(r_hat, v_hat)

        # Transformation Matrix
        M = np.vstack([r_hat, v_hat, z_hat])

        # Transformations
        p_conj = M @ r_rel
        C_comb = M @ (cov1 + cov2) @ M.T

        # Reducing to their projections in the conjunction plane
        p_conj = np.delete(p_conj, 1, axis=0)
        C_comb = np.delete(C_comb, 1, axis=0)
        C_comb = np.delete(C_comb, 1, axis=1)

        # Final integral
        # Calculated constants
        Pc_scaler = 1 / (2 * np.pi * np.sqrt(np.linalg.det(C_comb)))
        Pc_e_scaler = -0.5 * (p_conj @ np.linalg.inv(C_comb) @ p_conj)
        # After solving integral on paper
        Pc = Pc_scaler * np.pi * (HBR**2) * (np.e ** Pc_e_scaler)

        return Pc


    # 3D-PLotter
    # Source - http://www.youtube.com/@alfonsogonzalez-astrodynam2207
    def PlotOrbits(rs,labels,cb,title):
        # Defines a 3D PLot
        fig = plt.figure(figsize=(16,8))
        ax = fig.add_subplot(111,projection='3d')

        # Plot trajectory and starting point
        n = 0
        for r in rs:
            ax.plot(r[:,0], r[:,1], r[:,2], label=labels[n])
            ax.scatter(r[0,0], r[0,1], r[0,2], s=20, marker='o',label='Initial Position of '+labels[n])
            n += 1

        # Plot earth
        _u,_v = np.mgrid [0:2*np.pi:200j, 0:np.pi:100j]
        _x = cb['radius'] * np.cos(_u) * np.sin(_v)
        _y = cb['radius'] * np.sin(_u) * np.sin(_v)
        _z = cb['radius'] * np.cos(_v)
        ax.plot_surface(_x, _y, _z, rstride=4, cstride=4, color='b', alpha=0.3, linewidth=0)

        # Check for custom axes Limits
        max_val=np.max(np.abs(rs))
        # Set Labels and title
        ax.set_xlim([-max_val, max_val])
        ax.set_ylim([-max_val, max_val])
        ax.set_zlim([-max_val, max_val])
        ax.set_xlabel('X (km)') 
        ax.set_ylabel('Y (km)')
        ax.set_zlabel('Z (km)')
        ax.set_aspect('equal')
        ax.set_title(title)
        plt.legend()
        plt.show()
        pass


    # Keplerian-Orbital elements Plotter
    # Standard Unit of Time is Hours.
    # Source - http://www.youtube.com/@alfonsogonzalez-astrodynam2207
    def PlotKepOrbits(rs,vs,ts,cb,title):
        # Defining a plot
        fig, axs = plt.subplots(nrows=2, ncols=3, figsize=(16, 8))
        fig.suptitle(title,fontsize=25)
        xlabel = 'Time Elapsed (hours)'

        # Calculating Keplarian Orbital elements
        kep_elements = OrbitPropagator.State2Kep(rs, vs, cb)

        # Making time array
        ts = ts / 3600

        # Plotting Semi-major axis
        axs[0,0].plot(ts,kep_elements[:,0])
        axs[0,0].set_title('Semi-Major Axis vs. Time')
        axs[0,0].grid(False)
        axs[0,0].set_ylabel('a (km)')
        axs[0,0].set_xlabel(xlabel)

        axs[0,1].plot(ts,kep_elements[:,1])
        axs[0,1].set_title('Eccentricity vs. Time')
        axs[0,1].grid(False)
        axs[0,1].set_ylabel('e')
        axs[0,1].set_xlabel(xlabel)

        axs[0,2].plot(ts,kep_elements[:,2])
        axs[0,2].set_title('Inclination vs. Time')
        axs[0,2].grid(False)
        axs[0,2].set_ylabel('i (degrees)')
        axs[0,2].set_xlabel(xlabel)

        axs[1,0].plot(ts,kep_elements[:,3])
        axs[1,0].set_title('Right Ascension of the Ascending Node vs. Time')
        axs[1,0].grid(False)
        axs[1,0].set_ylabel('RAAN (degrees)')
        axs[1,0].set_xlabel(xlabel)

        axs[1,1].plot(ts,kep_elements[:,4])
        axs[1,1].set_title('Argument of Periapsis vs. Time')
        axs[1,1].grid(False)
        axs[1,1].set_ylabel('ω (degrees)')
        axs[1,1].set_xlabel(xlabel)

        axs[1,2].plot(ts,kep_elements[:,5])
        axs[1,2].set_title('True Anomaly vs. Time')
        axs[1,2].grid(False)
        axs[1,2].set_ylabel('v (degrees)')
        axs[1,2].set_xlabel(xlabel)

        fig.tight_layout(rect=[0, 0, 1, 1])
        plt.show()
        pass

   
    # Sun Position Calculator
    # Refrence fram is ICRF.
    # Note - The starting reference point is (2000-01-01 12:00:00 UTC)
    # Source - https://naif.jpl.nasa.gov/naif/data.html
    def sun_position_icrf(et):
        state, _ = spice.spkez(10, (et-64.18392728473108), "J2000", "NONE", 399) 
        # Note - 64ish seconds have been subtracted so that the initial time point is at (2000-01-01 12:00:00 UTC).
        pos = np.array(state[:3])   # km
        vel = np.array(state[3:6])  # km/s
        return pos, vel