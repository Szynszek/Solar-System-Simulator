import json
import spiceypy as spice
import numpy as np
import sys

celestial_bodies = {
    "Sun": {"body_id": "10"},
    "Mercury": {"body_id": "199"},
    "Venus": {"body_id": "299"},
    "Earth": {"body_id": "399"},
    "Mars_barycenter": {"body_id": "4"},
    "Jupiter_barycenter": {"body_id": "5"},
    "Saturn_barycenter": {"body_id": "6"},
    "Uranus_barycenter": {"body_id": "7"},
    "Neptune_barycenter": {"body_id": "8"},
    "Moon": {"body_id": "301"},
    "Ceres": {"body_id": "2000001"},
    "Pallas": {"body_id": "2000002"},
    "Vesta": {"body_id": "2000004"},
}

def generate_solar_system(date_str: str, output_file: str, celestial_bodies: dict):
    spice.furnsh("de440.bsp") # JPL DE440 planetary ephemerides
    spice.furnsh("gm_de440.tpc") # Gravitational parameters
    spice.furnsh("naif0012.tls") # Leapseconds kernel for time scale conversions
    spice.furnsh("pck00011.tpc") # Planetary constants kernel (orientation and shape)
    spice.furnsh("sb441-n16.bsp") # Kernel of the biggest asteroids
    spice.furnsh("spice-sim-data.tpc") # Extended planetary constants from my GitHub

    et = spice.str2et(f"{date_str} TDB")
    ref = "ECLIPJ2000"
    output_data = {}

    for name, ids in celestial_bodies.items():
        body_id = ids["body_id"]
        state = get_state(body_id, et, ref)
        GM, radii, j2 = get_gravitational_parameters(body_id)
        longitude, latitude = get_pole_angles(name, et, ref)

        output_data[name] = {
            "mu": GM,
            "reference_frame": ref,
            "longitude": longitude,
            "latitude": latitude,
            "radius": radii,
            "j2": j2,
            "state": (np.array(state)).tolist()
        }

    spice.kclear()

    with open(output_file, "w") as json_file:
        json.dump(output_data, json_file, indent=4)

def get_state(body_id: str, et: float, ref):
    # relative to the Solar System Barycenter
    try:
        state, _ = spice.spkezr(body_id, et, ref, "NONE", "0")
        state = state * 1e3
        return state
    except spice.stypes.SpiceyError:
        spice.kclear()
        sys.exit(f"FATAL ERROR: No state for {body_id} !")

def get_gravitational_parameters(body_id: str):
    GM = 0.0
    radii = 0.0
    j2 = 0.0
    error_list = []
    try:
        GM = spice.bodvrd(body_id, "GM", 1)[1][0] * 1e9
    except spice.stypes.SpiceyError:
        spice.kclear()
        sys.exit(f"FATAL ERROR: No gravitational parameters for {body_id} !")

    try:
        j2 = spice.bodvrd(body_id, "J2", 1)[1][0]
    except spice.stypes.SpiceyError:
        error_list.append("J2")

    try:
        radii = spice.bodvrd(body_id, "RADII", 3)[1][0] * 1000.0
    except spice.stypes.SpiceyError:
        error_list.append("RADII")

    if error_list:
        print(f"Info: No {error_list} for {body_id}")
    return GM, radii, j2

def get_pole_angles(body_name, et: float, ref: str) -> tuple[float, float]:

    try:
        frame_name = f"IAU_{body_name.upper()}"

        # Extract 3x3 rotation matrix from the ECLIPTIC J2000 frame
        rot_mat = spice.pxform(ref, frame_name, et)
        pole_ecliptic = rot_mat.T @ np.array([0.0, 0.0, 1.0])

        # spice.recrad returns: radius (ignored), longitude (RA equivalent), latitude (DEC equivalent)
        _, longitude, latitude = spice.recrad(pole_ecliptic)

        return np.degrees(longitude), np.degrees(latitude)

    except spice.stypes.SpiceyError:
        print(f"Info: Spatial orientation model unmodeled for {body_name.upper()} .")
        return 0.0, 90.0

if __name__ == "__main__":
    print("Generating solar system JSON file...")
    generate_solar_system("2026-03-14", "ephemeris.json", celestial_bodies)
    print("Finished without exceptions.")

