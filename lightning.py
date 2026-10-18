"""Rayos de la última hora sobre Sudamérica austral, para la app Celaje.

Corre en GitHub Actions cada ~5 minutos. Lee los archivos del detector de rayos GLM del
satélite GOES-19 (NOAA, dominio público, en la nube pública de AWS: un archivo cada 20 s),
baja solo los nuevos desde la vuelta anterior y escribe out/lightning.json:

  {"updated": unix, "last": "<inicio del último archivo leído>",
   "flashes": [[lat*100, lon*100, unix], ...]}   # solo los de la última hora
"""

import json
import os
import re
import sys
import tempfile
import time
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timedelta, timezone

import h5py
import numpy as np

NORTH, SOUTH, WEST, EAST = -5.0, -60.0, -95.0, -45.0
WINDOW = 3600
BUCKET = "https://noaa-goes19.s3.amazonaws.com"
PREVIOUS = os.environ.get("PREVIOUS_URL", "https://jorgedsgn.github.io/celaje-datos/lightning.json")
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")


def get(url, timeout=30):
    for attempt in range(3):
        try:
            with urllib.request.urlopen(url, timeout=timeout) as r:
                return r.read()
        except Exception:
            if attempt == 2:
                raise
            time.sleep(2)


def start_of(key):
    """Hora de inicio del archivo, desde su nombre (..._sAAAADDDHHMMSSs_...)."""
    m = re.search(r"_s(\d{4})(\d{3})(\d{2})(\d{2})(\d{2})(\d)_", key)
    y, doy, hh, mm, ss, tenth = (int(g) for g in m.groups())
    return datetime(y, 1, 1, tzinfo=timezone.utc) + timedelta(days=doy - 1, hours=hh, minutes=mm, seconds=ss + tenth / 10)


def list_keys(hour):
    prefix = hour.strftime("GLM-L2-LCFA/%Y/%j/%H/")
    xml = get(f"{BUCKET}/?list-type=2&prefix={prefix}").decode()
    return re.findall(r"<Key>([^<]+\.nc)</Key>", xml)


def flashes_in(key):
    data = get(f"{BUCKET}/{key}")
    with tempfile.NamedTemporaryFile(suffix=".nc") as tmp:
        tmp.write(data)
        tmp.flush()
        with h5py.File(tmp.name, "r") as f:
            if "flash_lat" not in f or f["flash_lat"].shape[0] == 0:
                return []
            lat = f["flash_lat"][:].astype(float)
            lon = f["flash_lon"][:].astype(float)
            ok = f["flash_quality_flag"][:] == 0 if "flash_quality_flag" in f else np.ones(len(lat), bool)
    t = int(start_of(key).timestamp())
    keep = ok & (lat <= NORTH) & (lat >= SOUTH) & (lon >= WEST) & (lon <= EAST)
    return [[int(round(a * 100)), int(round(o * 100)), t] for a, o in zip(lat[keep], lon[keep])]


def main():
    now = datetime.now(timezone.utc)
    cutoff = int(now.timestamp()) - WINDOW
    flashes, last = [], None
    try:
        prev_file = os.environ.get("PREVIOUS_FILE")
        prev = json.load(open(prev_file)) if prev_file else json.loads(get(PREVIOUS))
        flashes = [f for f in prev.get("flashes", []) if f[2] >= cutoff]
        last = prev.get("last")
    except Exception as e:
        print("Sin archivo anterior:", e)

    keys = []
    for h in (now - timedelta(hours=1), now):
        try:
            keys += list_keys(h)
        except Exception as e:
            print("No se pudo listar", h, e)
    keys = [k for k in keys if start_of(k).timestamp() >= cutoff and (last is None or os.path.basename(k) > last)]
    keys.sort()
    print(f"Archivos nuevos: {len(keys)}")
    with ThreadPoolExecutor(max_workers=8) as pool:
        for got in pool.map(lambda k: (k, safe(k)), keys):
            k, fl = got
            if fl is None:
                continue
            flashes += fl
            last = os.path.basename(k)

    flashes.sort(key=lambda f: f[2])
    os.makedirs(OUT, exist_ok=True)
    out = {"source": "NOAA GOES-19 GLM", "updated": int(now.timestamp()), "last": last, "flashes": flashes}
    with open(os.path.join(OUT, "lightning.json"), "w") as f:
        json.dump(out, f, separators=(",", ":"))
    print(f"Rayos en la última hora: {len(flashes)}")


def safe(k):
    try:
        return flashes_in(k)
    except Exception as e:
        print("Falló", k, e)
        return None


if __name__ == "__main__":
    sys.exit(main())
