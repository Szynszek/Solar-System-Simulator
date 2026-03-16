import json
import sys
import requests

def fetch_and_format_ephemeris(date: str) -> None:
    url = "https://ssd.jpl.nasa.gov/api/horizons.api"

    bodies_id = {
        "sun": "10", "mercury": "199", "venus": "299", "earth": "399",
        "mars": "499", "jupiter": "599", "saturn": "699", "uranus": "799",
        "neptune": "899", "moon": "301", "ceres": "1;", "pallas": "2;", "vesta": "4;"
    }

    mu_data = {
        "sun": 1.3271244004193938e+20, "mercury": 2.203186855000000e+13, "venus": 3.248585920000000e+14,
        "earth": 3.986004354360000e+14, "mars": 4.282837566200000e+13, "jupiter": 1.266865319000000e+17,
        "saturn": 3.793120623400000e+16, "uranus": 5.793950610300000e+15, "neptune": 6.835099970000000e+15,
        "moon": 4.902800066000000e+12, "ceres": 6.262840000000000e+10, "pallas": 1.363000000000000e+10,
        "vesta": 1.728828000000000e+10,
    }

    base_parameters = {
        "format": "text", "MAKE_EPHEM": "'YES'", "EPHEM_TYPE": "'VECTORS'",
        "CENTER": "'500@0'", "TLIST": f"'{date}'", "VEC_TABLE": "'2'",
        "REF_SYSTEM": "'ICRF'", "REF_PLANE": "'ECLIPTIC'", "VEC_CORR": "'NONE'",
        "CAL_TYPE": "'M'", "OUT_UNITS": "'KM-S'", "VEC_LABELS": "'YES'",
        "VEC_DELTA_T": "'NO'", "CSV_FORMAT": "'YES'", "OBJ_DATA": "'YES'"
    }

    output_data = {}

    session = requests.Session()

    for name, body_id in bodies_id.items():
        parameters = base_parameters.copy()
        parameters["COMMAND"] = f"'{body_id}'"

        try:
            response = session.get(url, params=parameters, timeout=10)
            response.raise_for_status()
            text_data = response.text.replace('\r\n', '\n')

            if '$$SOE' not in text_data:
                raise ValueError(f"JPL API Error for {name}.")

            raw_data = text_data.split('$$SOE\n')[1].split('\n$$EOE')[0].strip()
            values = raw_data.split(',')

            state_vector = [
                float(values[2]) * 1e3, float(values[3]) * 1e3, float(values[4]) * 1e3,
                float(values[5]) * 1e3, float(values[6]) * 1e3, float(values[7]) * 1e3
            ]
            output_data[name] = {
                "mu": mu_data[name],
                "state": state_vector
            }

        except Exception as e:
            sys.exit(f"Critical exception: {name}. Description: {e}")

    with open("ephemeris.json", "w") as json_file:
        json.dump(output_data, json_file, indent=4)

if __name__ == "__main__":
    target_date = "2026-03-14"
    print(f"% Data generated for date: {target_date}\n")
    fetch_and_format_ephemeris(target_date)
    print("Success: Ephemeris fully generated without exceptions.")