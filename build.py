"""Genera el mapa de lluvia de Celaje a partir de los datos abiertos del ECMWF.

Corre solo en GitHub Actions cada 3 horas. Baja la lluvia acumulada (tp) y la energía
para tormentas (mucape) del último pronóstico, recorta Sudamérica austral y guarda un
archivo chico (out/rain.bin) que la app descarga.

Datos: ECMWF Open Data, licencia CC BY 4.0 (uso comercial permitido citando la fuente).

Formato de rain.bin:
  4 bytes  largo N del encabezado (uint32, little endian)
  N bytes  encabezado JSON (zona, grilla, horas)
  resto    valores comprimidos (deflate sin encabezado): por cada intervalo de 3 h y cada
           punto (fila norte→sur, columna oeste→este) un byte:
           bits 0-6 = round(sqrt(mm/h) * 20)  (0 = sin lluvia, 127 ≈ 40 mm/h)
           bit 7    = tormenta probable (mucape alto con lluvia)
"""

import json
import math
import os
import struct
import sys
import zlib
from datetime import timedelta

import eccodes
import numpy as np
from ecmwf.opendata import Client

NORTH, SOUTH, WEST, EAST = -10.0, -60.0, -90.0, -50.0
RES = 0.25
STEPS = list(range(0, 63, 3))           # 0..60 h: alcanza para 48 h hacia adelante
STORM_CAPE = 400.0                      # J/kg
STORM_RAIN = 0.4                        # mm/h

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")


def main():
    client = Client(source=os.environ.get("ECMWF_SOURCE", "ecmwf"), model="ifs", resol="0p25")
    run = client.latest(type="fc", step=STEPS[-1], param="tp")
    print("Pronóstico", run)
    os.makedirs(OUT, exist_ok=True)
    grib = os.path.join(OUT, "data.grib2")
    client.retrieve(date=run, type="fc", step=STEPS, param=["tp", "mucape"], target=grib)

    rows = int(round((NORTH - SOUTH) / RES)) + 1
    cols = int(round((EAST - WEST) / RES)) + 1
    tp = {}
    cape = {}
    with open(grib, "rb") as f:
        while True:
            h = eccodes.codes_grib_new_from_file(f)
            if h is None:
                break
            try:
                param = eccodes.codes_get(h, "shortName")
                step = int(eccodes.codes_get(h, "endStep"))
                ni = eccodes.codes_get(h, "Ni")
                nj = eccodes.codes_get(h, "Nj")
                lat0 = eccodes.codes_get(h, "latitudeOfFirstGridPointInDegrees")
                lon0 = eccodes.codes_get(h, "longitudeOfFirstGridPointInDegrees")
                dlat = eccodes.codes_get(h, "jDirectionIncrementInDegrees")
                dlon = eccodes.codes_get(h, "iDirectionIncrementInDegrees")
                vals = eccodes.codes_get_values(h).reshape(nj, ni)
            finally:
                eccodes.codes_release(h)
            # Filas: de norte a sur. Columnas: la longitud puede partir en -180 o en 0.
            r0 = int(round((lat0 - NORTH) / dlat))
            sub_rows = vals[r0:r0 + rows]
            lon_idx = [int(round(((WEST + c * RES) - lon0) / dlon)) % ni for c in range(cols)]
            sub = sub_rows[:, lon_idx]
            (tp if param == "tp" else cape)[step] = sub

    steps = [s for s in STEPS if s in tp]
    intervals = []
    for a, b in zip(steps, steps[1:]):
        mmh = np.clip((tp[b] - tp[a]) * 1000.0 / (b - a), 0, None)     # tp viene en metros
        c = np.maximum(cape.get(a, 0), cape.get(b, 0))
        code = np.clip(np.round(np.sqrt(mmh) * 20), 0, 127).astype(np.uint8)
        storm = (c >= STORM_CAPE) & (mmh >= STORM_RAIN)
        code = code | (storm.astype(np.uint8) << 7)
        intervals.append((a, code))

    header = {
        "source": "ECMWF Open Data (CC BY 4.0)",
        "run": int(run.timestamp()),
        "north": NORTH, "west": WEST, "step": RES,
        "rows": rows, "cols": cols,
        # Hora de inicio (unix) de cada intervalo de 3 h.
        "starts": [int((run + timedelta(hours=a)).timestamp()) for a, _ in intervals],
        "hours": 3,
        "encoding": "sqrt20+storm",
    }
    payload = b"".join(code.tobytes() for _, code in intervals)
    comp = zlib.compressobj(9, zlib.DEFLATED, -15)
    body = comp.compress(payload) + comp.flush()
    head = json.dumps(header, separators=(",", ":")).encode()
    with open(os.path.join(OUT, "rain.bin"), "wb") as f:
        f.write(struct.pack("<I", len(head)) + head + body)
    os.remove(grib)
    rain_pts = sum(int((code & 127 > 6).sum()) for _, code in intervals)
    print(f"Listo: {len(intervals)} intervalos, {rows}x{cols} puntos, {len(body) / 1024:.0f} KB, puntos con lluvia: {rain_pts}")


if __name__ == "__main__":
    sys.exit(main())
