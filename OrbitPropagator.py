# Main Orbit Propagator File
import numpy as np
import matplotlib.pyplot as plt 
from scipy.integrate import ode
from scipy.special import erf
import Constants as ct
import ussa1976
import spiceypy as spice
import os
import sys


# For selecting which type of perturbation is needed.
def all_perturbations():
    return {
        'J2' : True,
        'J3' : True,
        'Drag' : True,
        'SRP' : False,
        'External Thrust' : False
    }


# Kernels for NASA Spice files.
# Do not touch this portion of code
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
# Reference Frame is ICRF.
class OrbitPropagator:
    # Integrater Function
    # Using DOP853 method.
    # Source - http://www.youtube.com/@alfonsogonzalez-astrodynam2207
    def __init__(self,r0,v0,tspan,dt,cenb,ob,pointing_mode,perturbations=all_perturbations()):
        load_kernels()
        self.r0 = r0 # Position Vector [km]
        self.v0 = v0 # Velocity Vector [km/s]
        self.tspan = tspan # Total time of propagation [s]
        self.dt = dt # Time step [s]
        self.cenb = cenb # Central Body (This extracts data of given object from Constants.py file)
        self.ob = ob # Orbiting Body (This extracts data of given object from Constants.py file)

        # Pointing Mode
        # As of now this only supports nadir pointing.
        self.pointing_mode = pointing_mode
        if self.pointing_mode=="nadir":
            if self.ob['shape'] == 'sphere_debri':
                print("Invalid Combination of pointing mode and shape of orbiting body")
                sys.exit()
        if self.pointing_mode=="N/A":
            if self.ob['shape'] == 'cube':
                print("Invalid Combination of pointing mode and shape of orbiting body")
                sys.exit()

        # Altitude data
        self.altitude = np.linspace(0,1000000,1000001)
        self.altitude_ds = ussa1976.compute(z=self.altitude, variables=["p", "rho", "t"])
        self.rho_table = self.altitude_ds["rho"].values.squeeze() # Atmospheric Density table [kg/m^3]

        # Total number of steps
        self.n_steps=int(np.ceil(self.tspan/self.dt))+1

        # Initializing arrays for state vector and time
        self.ys=np.zeros((self.n_steps,6)) # State Vector
        self.ts=np.zeros((self.n_steps,1)) # Time Vector
        self.os = np.zeros((self.n_steps, 3, 3)) # Orientation

        # Initial orientation
        if self.pointing_mode=="nadir":
            self.C0 = OrbitPropagator.compute_lvlh_frame(self.r0, self.v0)
        if self.pointing_mode=="N/A":
            self.C0 = np.zeros((3, 3))# This is just a placeholder.

        # Initial conditions
        self.y0 = np.hstack((self.r0, self.v0))
        self.ts[0] = 0
        self.ys[0] = np.array(self.y0)
        self.os[0] = self.C0
        self.C_cur = self.C0.copy()

        # Initiate solver
        self.step=1
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
            self.os[self.step] = self.C_cur.copy()
            self.step += 1
            # Use only for debugging
            print(f"{'─'*60}")
            print(f"  t   = {self.solver.t:.6f} s")
            print(f"  r   = [{self.solver.y[0]:.6f}, {self.solver.y[1]:.6f}, {self.solver.y[2]:.6f}] km")
            print(f"  v   = [{self.solver.y[3]:.6f}, {self.solver.y[4]:.6f}, {self.solver.y[5]:.6f}] km/s")
            print(f"  C (in inertial)  =")
            print(f"       [{self.C_cur[0,0]:.4f}, {self.C_cur[0,1]:.4f}, {self.C_cur[0,2]:.4f}]")
            print(f"       [{self.C_cur[1,0]:.4f}, {self.C_cur[1,1]:.4f}, {self.C_cur[1,2]:.4f}]")
            print(f"       [{self.C_cur[2,0]:.4f}, {self.C_cur[2,1]:.4f}, {self.C_cur[2,2]:.4f}]")
            print(f"{'─'*60}")
        self.rs = self.ys[:, 0:3]
        self.vs = self.ys[:, 3:6]

    # Simple Helper function to compute LVLH frame
    # Source - https://ai-solutions.com/_freeflyeruniversityguide/attitude_reference_frames.htm 
    def compute_lvlh_frame(r, v):
        r_hat = r / np.linalg.norm(r)
        h_hat = np.cross(r, v)/np.linalg.norm(np.cross(r, v))
        z = -r_hat
        y = -h_hat
        x = np.cross(y, z)
        x = x/np.linalg.norm(x) # To avoid small numerical errors
        return np.vstack((x, y, z)).T


    # Simple Helper function to compute VNB frame
    # Source - https://ai-solutions.com/_freeflyeruniversityguide/attitude_reference_frames.htm 
    def compute_vnb_frame(r, v):
        v_hat = v / np.linalg.norm(v)
        n_hat = np.cross(r, v) / np.linalg.norm(np.cross(r, v))
        b_hat = np.cross(v_hat, n_hat)
        b_hat = b_hat / np.linalg.norm(b_hat)
        return np.vstack((v_hat, n_hat, b_hat)).T


    # Force Model Function
    # Two Body acceleration is always on.
    # Source - http://www.youtube.com/@alfonsogonzalez-astrodynam2207
    def diffy_q(self,t,y):
        # Making proper arrays
        self.rx, self.ry, self.rz = y[0:3]
        self.vx, self.vy, self.vz = y[3:6]
        self.r = np.array([self.rx, self.ry, self.rz])
        self.v = np.array([self.vx, self.vy, self.vz])

        # Norm of the position vector
        self.norm_r = np.linalg.norm(self.r)

        # Two-body acceleration
        self.a = -self.r * self.cenb['mu'] / self.norm_r**3

        # Orientation matrix
        # Nadir pointing - Corresponds to LVLH frmae as, z-vector is pointed towards negative of position vector.
        if self.pointing_mode == 'nadir':
            # Getting desired Nadir (LVLH) orientation
            self.C_des = OrbitPropagator.compute_lvlh_frame(self.r, self.v)
            # Getting current achievable orientation via GetOrien function
            self.C_cur = self.GetOrien(self.C_cur, self.C_des, self.ob)   # C_cur is current orientation.

        # J2 perturbation
        if self.perturbations['J2']:
            self.j2x = self.r[0] * ((5 * self.r[2]**2 / self.norm_r**2) - 1)
            self.j2y = self.r[1] * ((5 * self.r[2]**2 / self.norm_r**2) - 1)
            self.j2z = self.r[2] * ((5 * self.r[2]**2 / self.norm_r**2) - 3)
            self.a_j2 = (1.5 * self.cenb['J2'] * self.cenb['mu'] * self.cenb['radius']**2 / self.norm_r**5) * np.array([self.j2x, self.j2y, self.j2z])
            self.a += self.a_j2

        # J3 perturbation
        if self.perturbations['J3']:
            self.j3x = (5 * self.r[0]) * ((7 * self.r[2]**3 / self.norm_r**3) - (3 * self.r[2] / self.norm_r))
            self.j3y = (5 * self.r[1]) * ((7 * self.r[2]**3 / self.norm_r**3) - (3 * self.r[2] / self.norm_r))
            self.j3z = ((35 * self.r[2]**4 / self.norm_r**4) - (30 * self.r[2]**2 / self.norm_r**2) + 3)
            self.a_j3 = (0.5 * self.cenb['J3'] * self.cenb['mu'] * self.cenb['radius']**3 / self.norm_r**6) * np.array([self.j3x, self.j3y, self.j3z])
            self.a += self.a_j3

        # Drag perturbation
        # Source - https://arc.aiaa.org/doi/suppl/10.2514/1.A33606
        if self.perturbations['Drag']:
            self.z = (self.norm_r - self.cenb['radius']) * 1000.0
            if self.z < 1000000:
                self.vrel = (self.v - np.cross(self.cenb['at_rot_vec'], self.r))*1000 # m/s
                self.rho_inf = self.rho_table[round(self.z)]

                # Sphere
                if self.ob['shape'] == 'sphere_debri':
                    self.C_d = OrbitPropagator.get_DragCoefficient_Sphere(self.norm_r, np.linalg.norm(self.vrel), self.ob, self.cenb)
                    self.a_drag = (-((0.5 * self.rho_inf * self.C_d * (np.pi*self.ob['radius']**2)) / self.ob['mass']) * np.linalg.norm(self.vrel) * (self.vrel))/1000
                
                # Cube
                if self.ob['shape'] == 'cube':
                    # Nadir Pointing
                    if self.pointing_mode == 'nadir':
                            self.A_ref_drag  = self.ob["ly"] * self.ob["lz"]
                            self.v_body = self.C_cur.T @ self.vrel # inertial to body (m/s)

                            # Getting angle-of-attack and side-slip angles
                            # Source - https://www.youtube.com/watch?v=4kaK569ug9Q
                            self.alpha = np.arctan2(self.v_body[2], self.v_body[0])
                            self.beta = np.arcsin(self.v_body[1] / np.linalg.norm(self.v_body))

                            # Force Coefficients
                            self.C_A, self.C_S, self.C_N = OrbitPropagator.get_AeroCoefficient_RectangularPrism(np.rad2deg(self.alpha), np.rad2deg(self.beta), self.norm_r, np.linalg.norm(self.vrel), self.ob, self.cenb)
                            # Force calculation
                            self.q_inf = 0.5 * self.rho_inf * np.linalg.norm(self.v_body)**2 # N/m²
                            self.F_body_axial = -self.q_inf * self.A_ref_drag * self.C_A # N
                            self.F_body_side = -self.q_inf * self.A_ref_drag * self.C_S # N
                            self.F_body_normal = -self.q_inf * self.A_ref_drag * self.C_N # N

                            # Rotate body → inertial
                            # C_cur: body → inertial (columns = body axes in inertial)
                            F_body_vec  = np.array([self.F_body_axial, self.F_body_side, self.F_body_normal])
                            F_inertial  = self.C_cur @ F_body_vec # N, inertial frame
                            self.a_drag = (F_inertial / self.ob['mass']) / 1000 # km/s^2

                self.a += self.a_drag
            else:
                self.a = self.a

        # Solar Radiation Pressure
        if self.perturbations['SRP']:
            self.a = self.a # Work in progress

        # This is the custom code block for the external thrust
        if self.perturbations['External Thrust']:
            self.a = self.a

        # State Derivative
        self.dydt = np.zeros(6)
        self.dydt[0:3] = [self.vx, self.vy, self.vz]
        self.dydt[3:6] = self.a
        return self.dydt
 

    # Source - https://academicflight.com/articles/kinematics/rotation-formalisms/principal-rotation-vector/
    def GetOrien(self, Orien_i, Orien_des, ob):
        self.Orien_i = Orien_i
        self.Orien_des = Orien_des
        self.ob = ob
        
        # Skipping this function for sphere, as it represents debri.
        if self.ob['shape'] == 'sphere_debri':
            return self.Orien_des.copy()
        else:
            # Compute rotation error
            self.R_err = self.Orien_i.T @ self.Orien_des
            # Principal rotation angle
            self.rot_theta = np.arccos( max(-1.0, min(1.0, (0.5 * (np.trace(self.R_err)-1)) )))
            # If angle is very small, it will return the initial orientation.
            if np.abs(self.rot_theta) < 1e-12:
                return self.Orien_des.copy()
            # Principal rotation axis
            self.rot_axis = np.array([ self.R_err[1,2] - self.R_err[2,1],
                        self.R_err[2,0] - self.R_err[0,2],
                        self.R_err[0,1] - self.R_err[1,0]]) / (2*np.sin(self.rot_theta))
            
            # Calculating the rotation rate about the principal rotation axis
            self.rot_theta_dot = self.rot_theta/self.dt
            # Calculating the angular rate component-wise
            self.bodyrate = self.rot_theta_dot * self.rot_axis

            # If Maneuver is possible.
            if np.all(np.abs(self.bodyrate) <= self.ob['max_slew_rate']):
                return self.Orien_des.copy()
            
            # If Maneuver is not possible.
            else:
                self.N_matrix = np.array([[0, -self.rot_axis[2], self.rot_axis[1]],
                                        [self.rot_axis[2], 0, -self.rot_axis[0]],
                                        [-self.rot_axis[1], self.rot_axis[0], 0]])
                self.ndiyadn = np.array([[self.rot_axis[0]*self.rot_axis[0], self.rot_axis[0]*self.rot_axis[1], self.rot_axis[0]*self.rot_axis[2]],
                                        [self.rot_axis[1]*self.rot_axis[0], self.rot_axis[1]*self.rot_axis[1], self.rot_axis[1]*self.rot_axis[2]],
                                        [self.rot_axis[2]*self.rot_axis[0], self.rot_axis[2]*self.rot_axis[1], self.rot_axis[2]*self.rot_axis[2]]])
                self.theta_dot_max = np.min(np.abs(self.ob['Max_Slew_Rate'] / self.rot_axis))
                self.R_max = ( (np.cos(self.theta_dot_max*self.dt)*np.identity(3)) 
                            + (np.sin(self.theta_dot_max*self.dt)*self.N_matrix)
                            + ((1-np.cos(self.theta_dot_max*self.dt))*self.ndiyadn) )
                return self.Orien_i @ self.R_max
        

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
            
            # Invalid Cases
            if a*(1-e_mag) <= cb['radius']: # If perigee is smaller than central body's radius
                kep_array[k, :] = np.nan

        return kep_array


    # Collision Probabilty calculator Function
    # Source - https://amostech.com/TechnicalPapers/2014/Conjunction_Assessment/DUNCAN.pdf
    def Pc2D(r1,v1,cov1,r2,v2,cov2,HBR):
        r_rel = r1 - r2
        v_rel = v1 - v2
        r_hat = r_rel / np.linalg.norm(r_rel)
        v_hat = v_rel / np.linalg.norm(v_rel)
        z_hat = np.cross(r_hat, v_hat)
        if np.linalg.norm(z_hat) < 1e-10:
            raise ValueError("Degenerate conjunction geometry: r_rel and v_rel are nearly parallel")
        z_hat = z_hat / np.linalg.norm(z_hat)

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
        Pc = Pc_scaler * np.pi * (HBR**2) * np.exp(Pc_e_scaler)

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
            # ax.scatter(r[0,0], r[0,1], r[0,2], s=20, marker='o',label='Initial Position of '+labels[n])
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
        ax.grid(False)
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
    

    # Energy vs. Time Plotter
    def PlotEnergy(rs,vs,ts,cb,title):
        # Defining a plot
        fig, ax = plt.subplots(figsize=(16, 8))
        fig.suptitle(title,fontsize=25)
        xlabel = 'Time Elapsed (hours)'

        # Calculating Energy
        kinetic = 0.5 * np.linalg.norm(vs, axis=1)**2
        potential = -cb['mu'] / np.linalg.norm(rs, axis=1)
        totel = (kinetic + potential) * 1000

        # Making time array
        ts = ts / 3600

        # Plotting Semi-major axis
        ax.plot(ts,totel)
        ax.grid(False)
        ax.set_ylabel('Specific Energy (kJ/kg)')
        ax.set_xlabel(xlabel)

        fig.tight_layout(rect=[0, 0, 1, 1])
        plt.show()
        pass


    # This function compute aerodynamic force coefficient in Body frame for a spacecraft at given orientation and at given state.
    # This function assumme that spacecraft is "Rectangular Prism".
    def get_AeroCoefficient_RectangularPrism(alpha, beta, r_mag, v_mag, ob, cenb=ct.Earth):
        """
        This function compute aerodynamic force coefficient in Body frame for a spacecraft at given orientation and at given state.
        This function assumme that spacecraft is "Rectangular Prism".
            
        Inputs
        ------
        alpha = Angle of Attack / Negative of Pitch angle (degrees)
        beta = Sideslip Angle / Negative of Yaw angle (degrees)
        r_mag = Position magnitude (km)
        v_mag = Velocity magnitude (km/s)
        ob = Orbiting Body
        cb = Central Body (Earth is default choice)
        
        Outputs
        -------
        C_A, C_S, C_N = aerodynamic force coefficient (Axial, Side, Normal)

        """
        # Object properties
        T_wall = ob["T_wall"]
        sigma_N = ob["sigma_N"]
        sigma_T = ob["sigma_T"]
        lx, ly, lz = ob["lx"], ob["ly"], ob["lz"]

        # Atmosphere
        ds    = ussa1976.compute(z=np.array([(r_mag-cenb["radius"]) * 1000]), variables=["p", "rho", "t"])
        rho   = ds["rho"].values[0]
        T_inf = ds["t"].values[0]
        p     = ds["p"].values[0]
        s     = (v_mag) / np.sqrt(2 * (p / rho))
        Tr    = T_wall / T_inf

        # Trig shorthands
        ca, sa = np.cos(np.radians(alpha)), np.sin(np.radians(alpha))
        cb, sb = np.cos(np.radians(beta)),  np.sin(np.radians(beta))

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

        return C_A, C_S, C_N


    # This function compute drag coefficient for a spacecraft at given state.
    # This function assumme that spacecraft is "Sphere".
    def get_DragCoefficient_Sphere(r_mag, v_mag, ob, cenb=ct.Earth):
        """
        This function compute drag coefficient for a spacecraft at given state.
        This function assumme that spacecraft is "Sphere".
            
        Inputs
        ------
        r_mag = Position magnitude (km)
        v_mag = Velocity magnitude (km/s)
        ob = Orbiting Body
        cb = Central Body (Earth is default choice)
        
        Outputs
        -------
        C_D = Drag coefficient

        """
        # Object properties
        T_wall = ob["T_wall"]
        sigma_N = ob["sigma_N"]
        sigma_T = ob["sigma_T"]

        # Atmosphere at given altitude
        ds = ussa1976.compute(z=np.array([(r_mag-cenb["radius"]) *1000]), variables=["p", "rho", "t"])
        rho   = ds["rho"].values[0]
        T_inf = ds["t"].values[0]
        p     = ds["p"].values[0]
        s     = (v_mag) / np.sqrt(2 * (p / rho))

        # C_D formula
        term1 = (2 - sigma_N + sigma_T) * (0.5 +
            0.5 * erf(s) * (1 + 1/s**2 - 1/(4*s**4)) +
            (1 + 2*s**2) * np.exp(-s**2) / (4*s**3*np.sqrt(np.pi)))
        term2 = (sigma_N / (3*s)) * np.sqrt(np.pi * (T_wall/T_inf)) * (1 + erf(s))
        term3 = (2 - sigma_N)/(2*s**2) 
        term4 = (sigma_N/(6*s**4)) * (1 + (2*s**2 - 1)*np.exp(-s**2)) * np.sqrt(T_wall/T_inf)
        C_D = term1 + term2 + term3 + term4

        return C_D
