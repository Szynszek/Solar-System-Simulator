import requests
import sys


def fetch_and_format_ephemeris(date: str) -> None:
    url = "https://ssd.jpl.nasa.gov/api/horizons.api"

    bodies_id = {
        "sun": 10, "mercury": 199, "venus": 299, "earth": 399,
        "mars": 499, "jupiter": 599, "saturn": 699, "uranus": 799, "neptune": 899
    }

    base_parameters = {
        "format": "text", "MAKE_EPHEM": "'YES'", "EPHEM_TYPE": "'VECTORS'",
        "CENTER": "'500@0'", "TLIST": f"'{date}'", "VEC_TABLE": "'2'",
        "REF_SYSTEM": "'ICRF'", "REF_PLANE": "'ECLIPTIC'", "VEC_CORR": "'NONE'",
        "CAL_TYPE": "'M'", "OUT_UNITS": "'KM-S'", "VEC_LABELS": "'YES'",
        "VEC_DELTA_T": "'NO'", "CSV_FORMAT": "'YES'", "OBJ_DATA": "'YES'"
    }

    formatted_lines = []

    for name, body_id in bodies_id.items():
        parameters = base_parameters.copy()
        parameters["COMMAND"] = f"'{body_id}'"

        try:
            response = requests.get(url, params=parameters, timeout=10)
            response.raise_for_status()
            raw_data = response.text.split('$$SOE\n')[1].split('\n$$EOE')[0].strip()
            values = raw_data.split(',')

            x, y, z = float(values[2]) * 1e3, float(values[3]) * 1e3, float(values[4]) * 1e3
            vx, vy, vz = float(values[5]) * 1e3, float(values[6]) * 1e3, float(values[7]) * 1e3

            matlab_code = (
                f"r0_{name} = [{x}, {y}, {z}];\n"
                f"v0_{name} = [{vx}, {vy}, {vz}];\n"
                f"y0_{name} = [r0_{name}, v0_{name}];\n"
            )
            formatted_lines.append(matlab_code)

        except Exception as e:
            print(f"Critical exception: {name}. Description: {e}")
            sys.exit(1)


    for line in formatted_lines:
        print(line)

    y0_string = "y0 = [" + ", ".join([f"y0_{name}" for name in bodies_id.keys()]) + "];"
    print(y0_string)


if __name__ == "__main__":
    target_date = "2026-03-14"
    print(f"% Data generated for date: {target_date}\n")
    fetch_and_format_ephemeris(target_date)