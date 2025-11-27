# sun_spice.py
import os
import numpy as np
import spiceypy as spice

KERNEL_DIR = r"C:\Engineering Program Files\spice_kernels"

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

def sun_position_icrf(utc_string):
    """
    Return Sun position vector [km] in ICRF/J2000 relative to Earth's center.
    utc_string example: "2000-01-01T00:00:00"
    """
    # Convert UTC string to ephemeris seconds past J2000 (ET)
    et = spice.str2et(utc_string)      # str2et handles common ISO strings
    # spkez(target, et, ref_frame, abcorr, observer)
    state, lt = spice.spkez(10, (et), "J2000", "NONE", 399)
    pos = np.array(state[:3])   # km
    vel = np.array(state[3:6])  # km/s
    return pos, vel, lt

if __name__ == "__main__":
    load_kernels()
    utc = "2000-01-01T12:00:00"
    pos, vel, light_time = sun_position_icrf(utc)
    print("Sun position (km):", pos)
    print("Sun velocity (km/s):", vel)
    print("Light-time (s):", light_time)
    print(360+np.degrees(np.atan2(pos[1],pos[0])))
    print(-4.5e-5 * 1000)
    print(np.radians(23))