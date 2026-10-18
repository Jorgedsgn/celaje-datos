"""Pronóstico por hora de Chile con 3 modelos (ECMWF, GFS e ICON), para la app Celaje.

Corre en GitHub Actions cada 6 horas. Baja los datos abiertos de cada modelo, los lleva a
una grilla común de 0,25° sobre Chile, los pasa a horas (próximos ~7 días) y escribe un
archivo por cada zona de 2° × 2° en out/fc/ (la app baja solo el de su zona).

Fuentes (todas permiten uso comercial):
  ECMWF Open Data (CC BY 4.0) · NOAA GFS (dominio público) · DWD ICON (CC BY 4.0)
  Alturas reales: AWS Terrain Tiles (SRTM/GMTED/ETOPO), ver tools/make_dem.py.

Formato de cada archivo out/fc/<norte>_<oeste>.bin:
  uint32 largo del encabezado + encabezado JSON + deflate crudo de:
    por modelo, por variable, por hora, por punto (8×8): un byte (255 = sin dato)
    luego las alturas reales de la zona (int16, cada 1/60°)
  Variables y escala (byte = valor convertido):
    temp, dew: (°C + 64) × 2      precip: sqrt(mm/h) × 20     prob: %
    code: WMO                       gust, wind: km/h            dir: grados / 2
    cloud: %                        cape: J/kg / 20             fl: metros / 50
"""

import bz2
import json
import math
import os
import struct
import sys
import tempfile
import time
import urllib.request
import zlib
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timedelta, timezone

import eccodes
import numpy as np

NORTH, SOUTH, WEST, EAST, RES = -17.0, -56.0, -76.0, -66.0, 0.25
ROWS = int(round((NORTH - SOUTH) / RES)) + 1      # 157
COLS = int(round((EAST - WEST) / RES)) + 1        # 41
LATS = NORTH - np.arange(ROWS) * RES
LONS = WEST + np.arange(COLS) * RES
TILE = 8                                          # puntos por lado (2°)
HOURS = 192                                       # desde 12 h antes de ahora hasta ~7,5 días
VARS = ["temp", "dew", "precip", "prob", "code", "gust", "wind", "dir", "cloud", "cape", "fl"]
MODELS = ["ecmwf", "gfs", "icon"]
QUICK = os.environ.get("QUICK") == "1"     # prueba rápida con pocas horas
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")


def log(*a):
    print(*a, flush=True)


def get(url, timeout=60, tries=4):
    for k in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=timeout) as r:
                return r.read()
        except Exception as e:
            if k == tries - 1:
                raise
            time.sleep(2 + 3 * k)


def exists(url):
    try:
        req = urllib.request.Request(url, method="HEAD")
        with urllib.request.urlopen(req, timeout=30) as r:
            return r.status == 200
    except Exception:
        return False


# --------------------------------------------------------------------------- GRIB

def regular_subset(h):
    """Valores de un mensaje GRIB en grilla regular, recortados a nuestra grilla (vecino más cercano)."""
    ni, nj = eccodes.codes_get(h, "Ni"), eccodes.codes_get(h, "Nj")
    lat0 = eccodes.codes_get(h, "latitudeOfFirstGridPointInDegrees")
    lon0 = eccodes.codes_get(h, "longitudeOfFirstGridPointInDegrees")
    dlat = eccodes.codes_get(h, "jDirectionIncrementInDegrees")
    dlon = eccodes.codes_get(h, "iDirectionIncrementInDegrees")
    jpos = eccodes.codes_get(h, "jScansPositively")
    vals = eccodes.codes_get_values(h).reshape(nj, ni)
    vals = np.where(vals == eccodes.codes_get(h, "missingValue"), np.nan, vals) if eccodes.codes_get(h, "bitmapPresent") else vals
    sgn = 1 if jpos else -1
    r = np.clip(np.round((LATS - lat0) / (sgn * dlat)).astype(int), 0, nj - 1)
    c = np.round(((LONS - lon0) % 360) / dlon).astype(int) % ni
    return vals[np.ix_(r, c)]


def read_grib(path, on_message):
    with open(path, "rb") as f:
        while True:
            h = eccodes.codes_grib_new_from_file(f)
            if h is None:
                break
            try:
                on_message(h)
            finally:
                eccodes.codes_release(h)


# --------------------------------------------------------------------------- Horas

def build_hours(valid, inst, accum, t0, step_fn=None):
    """Lleva una serie por tiempos válidos a horas desde t0.

    valid: lista de datetimes ordenados. inst: dict var -> lista de arrays (instantáneos,
    se interpolan). accum: lista de arrays de lluvia acumulada (mm) desde el inicio del run,
    o None. Devuelve dict var -> array [HOURS, ROWS, COLS] (NaN si no hay dato).
    """
    secs = np.array([(v - t0).total_seconds() / 3600 for v in valid])
    out = {}
    hrs = np.arange(HOURS)
    for name, series in inst.items():
        arr = np.full((HOURS, ROWS, COLS), np.nan, np.float32)
        stack = np.stack(series)
        for k, h in enumerate(hrs):
            if h < secs[0] or h > secs[-1]:
                continue
            j = min(np.searchsorted(secs, h), len(secs) - 1)
            if secs[j] == h or j == 0:
                arr[k] = stack[j]
            else:
                a, b = secs[j - 1], secs[j]
                w = (h - a) / (b - a)
                arr[k] = stack[j - 1] * (1 - w) + stack[j] * w
        out[name] = arr
    if accum is not None:
        arr = np.full((HOURS, ROWS, COLS), np.nan, np.float32)
        acc = np.stack(accum)
        for j in range(1, len(secs)):
            a, b = secs[j - 1], secs[j]
            rate = np.clip((acc[j] - acc[j - 1]) / (b - a), 0, None)
            for h in hrs:
                # La hora h cubre (h-1, h]: se le asigna la tasa del intervalo que la contiene.
                if a < h <= b:
                    arr[h] = rate
        out["precip"] = arr
    return out


def wmo_code(precip, temp, cloud, cape, snow=None):
    """Código del tiempo (como Open-Meteo/WMO) a partir de las variables de una hora."""
    code = np.where(cloud < 20, 0, np.where(cloud < 50, 1, np.where(cloud < 85, 2, 3)))
    is_snow = (temp <= 0.8) if snow is None else (snow | (temp <= 0.3))
    convective = cape >= 300
    rain = np.where(precip < 0.3, 51, np.where(precip < 2.5, 61, np.where(precip < 7.5, 63, 65)))
    shower = np.where(precip < 2.5, 80, np.where(precip < 7.5, 81, 82))
    rain = np.where(convective & (precip >= 0.3), shower, rain)
    snowc = np.where(precip < 1, 71, np.where(precip < 3, 73, 75))
    wet = precip >= 0.1
    code = np.where(wet, np.where(is_snow, snowc, rain), code)
    thunder = (cape >= 800) & (precip >= 0.5)
    code = np.where(thunder, np.where((cape >= 2000) & (precip >= 3), 96, 95), code)
    return code


def neighborhood_prob(precip):
    """Probabilidad de lluvia de un modelo: % de puntos con lluvia (≥0,1 mm/h) en 5×5 vecinos y ±1 h."""
    wet = (np.nan_to_num(precip, nan=0) >= 0.1).astype(np.float32)
    valid = (~np.isnan(precip)).astype(np.float32)
    def box(a):
        p = np.pad(a, ((1, 1), (2, 2), (2, 2)))
        s = np.zeros_like(a)
        for dt in range(3):
            for dy in range(5):
                for dx in range(5):
                    s += p[dt:dt + a.shape[0], dy:dy + a.shape[1], dx:dx + a.shape[2]]
        return s
    prob = 100 * box(wet) / np.maximum(box(valid), 1)
    return np.where(np.isnan(precip), np.nan, prob)


def wind_from(u, v):
    speed = np.hypot(u, v) * 3.6
    deg = (np.degrees(np.arctan2(-u, -v)) + 360) % 360
    return speed, deg


# --------------------------------------------------------------------------- ECMWF

def ecmwf(t0):
    from ecmwf.opendata import Client
    client = Client(source=os.environ.get("ECMWF_SOURCE", "azure"), model="ifs", resol="0p25")
    latest = client.latest(type="fc", step=90, param="2t").replace(tzinfo=timezone.utc)
    runs = [latest]
    if latest.hour in (6, 18):
        runs.append(latest - timedelta(hours=6))   # el run anterior llega a 144 h+
    sfc = ["2t", "2d", "tp", "10u", "10v", "10fg", "tcc", "mucape", "sp", "msl", "ptype"]
    levels = [1000, 925, 850, 700, 600, 500]
    data = {}   # valid -> dict
    precip_hour = np.full((HOURS, ROWS, COLS), np.nan, np.float32)
    for k, run in enumerate(runs):
        if run.hour in (0, 12):
            steps = list(range(0, 145, 3)) + list(range(150, 187, 6))
        else:
            steps = list(range(0, 91, 3))
        if k == 1:
            # Solo lo que falta después del run nuevo (más un paso antes, para la lluvia del empalme).
            last_new = runs[0] + timedelta(hours=90)
            steps = [s for s in steps if run + timedelta(hours=s) >= last_new - timedelta(hours=6)]
        if QUICK:
            steps = steps[:5]
        log("ECMWF run", run, "pasos", len(steps))
        with tempfile.TemporaryDirectory() as d:
            p1, p2 = os.path.join(d, "sfc.grib2"), os.path.join(d, "pl.grib2")
            naive = run.replace(tzinfo=None)
            client.retrieve(date=naive, type="fc", step=steps, param=sfc, target=p1)
            client.retrieve(date=naive, type="fc", step=steps, param=["t", "gh"], levtype="pl", levelist=levels, target=p2)
            per = {}
            def take(h):
                step = int(eccodes.codes_get(h, "endStep"))
                key = eccodes.codes_get(h, "shortName")
                if eccodes.codes_get(h, "typeOfLevel") == "isobaricInhPa":
                    key = f"{key}{eccodes.codes_get(h, 'level')}"
                per.setdefault(step, {})[key] = regular_subset(h)
            read_grib(p1, take)
            read_grib(p2, take)
        # Lluvia: el acumulado solo sirve dentro del mismo run. El run nuevo tiene prioridad.
        ss = sorted(per)
        for i in range(1, len(ss)):
            a = (run + timedelta(hours=ss[i - 1]) - t0).total_seconds() / 3600
            b = (run + timedelta(hours=ss[i]) - t0).total_seconds() / 3600
            rate = np.clip((per[ss[i]].get("tp", 0) - per[ss[i - 1]].get("tp", 0)) * 1000 / (b - a), 0, None)
            for h in range(HOURS):
                if a < h <= b and np.isnan(precip_hour[h]).all():
                    precip_hour[h] = rate
        for s_ in ss:
            valid = run + timedelta(hours=s_)
            if valid not in data:
                data[valid] = per[s_]
    valid = sorted(data)
    def series(key, fn=lambda x: x):
        return [fn(data[v][key]) if key in data[v] else np.full((ROWS, COLS), np.nan) for v in valid]
    # Isoterma 0 °C desde los niveles de presión.
    def freezing(f):
        t2 = f["2t"] - 273.15
        fl = np.full((ROWS, COLS), np.nan)
        prev_t, prev_z = t2, np.full((ROWS, COLS), 0.0)
        found = np.zeros((ROWS, COLS), bool)
        for lv in levels:
            t = f[f"t{lv}"] - 273.15
            z = f[f"gh{lv}"]
            cross = (~found) & (prev_t > 0) & (t <= 0) & (z > prev_z)
            w = np.where(cross, prev_t / np.maximum(prev_t - t, 1e-6), 0)
            fl = np.where(cross, prev_z + w * (z - prev_z), fl)
            found |= cross
            prev_t, prev_z = np.where(z > prev_z, t, prev_t), np.maximum(z, prev_z)
        fl = np.where(found, fl, np.where(t2 <= 0, 0, 5600))
        return fl
    def elev(f):
        t = f["2t"]
        return np.clip(t / 0.0065 * ((f["msl"] / f["sp"]) ** 0.190263 - 1), 0, 6500)
    inst = {
        "temp": series("2t", lambda x: x - 273.15),
        "dew": series("2d", lambda x: x - 273.15),
        "u": series("10u"), "v": series("10v"),
        "gust": series("10fg", lambda x: x * 3.6),
        "cloud": series("tcc", lambda x: x * 100),
        "cape": series("mucape"),
        "fl": [freezing(data[v]) if "t1000" in data[v] else np.full((ROWS, COLS), np.nan) for v in valid],
        "snow": series("ptype", lambda x: ((x >= 5) & (x <= 8)).astype(float)),
    }
    out = build_hours(valid, inst, None, t0)
    out["precip"] = precip_hour
    out["wind"], out["dir"] = wind_from(out.pop("u"), out.pop("v"))
    snow = out.pop("snow") >= 0.5
    out["code"] = wmo_code(np.nan_to_num(out["precip"]), out["temp"], out["cloud"], np.nan_to_num(out["cape"]), snow)
    out["code"] = np.where(np.isnan(out["temp"]), np.nan, out["code"])
    first = data[valid[0]]
    return out, elev(first), runs[0]


# --------------------------------------------------------------------------- GFS

def gfs(t0):
    base = "https://nomads.ncep.noaa.gov"
    now = datetime.now(timezone.utc)
    run = None
    for back in range(0, 30, 6):
        c = (now - timedelta(hours=back)).replace(minute=0, second=0, microsecond=0)
        c = c.replace(hour=c.hour - c.hour % 6)
        url = f"{base}/pub/data/nccf/com/gfs/prod/gfs.{c:%Y%m%d}/{c:%H}/atmos/gfs.t{c:%H}z.pgrb2.0p25.f180.idx"
        if exists(url):
            run = c
            break
    if run is None:
        raise RuntimeError("GFS: no hay run completo")
    steps = list(range(0, 121)) + list(range(123, 181, 3))
    if QUICK:
        steps = steps[:4]
    log("GFS run", run, "pasos", len(steps))
    q = ("&var_TMP=on&var_DPT=on&var_UGRD=on&var_VGRD=on&var_GUST=on&var_APCP=on&var_TCDC=on"
         "&var_CAPE=on&var_HGT=on&lev_2_m_above_ground=on&lev_10_m_above_ground=on&lev_surface=on"
         "&lev_entire_atmosphere=on&lev_0C_isotherm=on&subregion=&toplat=-16&leftlon=283&rightlon=295&bottomlat=-57")
    per = {}
    orog = None
    def fetch(s):
        url = (f"{base}/cgi-bin/filter_gfs_0p25.pl?dir=%2Fgfs.{run:%Y%m%d}%2F{run:%H}%2Fatmos"
               f"&file=gfs.t{run:%H}z.pgrb2.0p25.f{s:03d}{q}")
        return s, get(url, timeout=90, tries=5)
    # NOAA bloquea si se piden más de ~120 archivos por minuto: de a 2 en paralelo alcanza.
    with ThreadPoolExecutor(max_workers=2) as pool:
        for s, blob in pool.map(fetch, steps):
            with tempfile.NamedTemporaryFile(suffix=".grib2") as tmp:
                tmp.write(blob)
                tmp.flush()
                f = {}
                def take(h):
                    nonlocal orog
                    sn, lvl = eccodes.codes_get(h, "shortName"), eccodes.codes_get(h, "typeOfLevel")
                    rng = str(eccodes.codes_get(h, "stepRange"))
                    if sn == "tp":
                        if rng.startswith("0-") or (s == 0):
                            f["tp"] = regular_subset(h)
                        return
                    if "-" in rng:        # promedios de intervalo (tcc): se usa el instantáneo
                        return
                    if sn == "orog":
                        orog = regular_subset(h)
                        return
                    key = {"2t": "temp", "2d": "dew", "10u": "u", "10v": "v", "gust": "gust", "tcc": "cloud",
                           "cape": "cape", "gh": "fl"}.get(sn)
                    if sn == "gh" and lvl != "isothermZero":
                        return
                    if key:
                        f[key] = regular_subset(h)
                read_grib(tmp.name, take)
                per[s] = f
    valid = [run + timedelta(hours=s) for s in steps]
    def series(k, fn=lambda x: x):
        return [fn(per[s][k]) if k in per[s] else np.full((ROWS, COLS), np.nan) for s in steps]
    inst = {"temp": series("temp", lambda x: x - 273.15), "dew": series("dew", lambda x: x - 273.15),
            "u": series("u"), "v": series("v"), "gust": series("gust", lambda x: x * 3.6),
            "cloud": series("cloud"), "cape": series("cape"), "fl": series("fl")}
    acc = [per[s].get("tp", np.zeros((ROWS, COLS))) for s in steps]
    out = build_hours(valid, inst, acc, t0)
    out["wind"], out["dir"] = wind_from(out.pop("u"), out.pop("v"))
    out["code"] = wmo_code(np.nan_to_num(out["precip"]), out["temp"], out["cloud"], np.nan_to_num(out["cape"]))
    out["code"] = np.where(np.isnan(out["temp"]), np.nan, out["code"])
    return out, orog if orog is not None else np.zeros((ROWS, COLS)), run


# --------------------------------------------------------------------------- ICON

def icon(t0):
    base = "https://opendata.dwd.de/weather/nwp/icon/grib"
    now = datetime.now(timezone.utc)
    run = None
    for back in range(0, 30, 6):
        c = (now - timedelta(hours=back)).replace(minute=0, second=0, microsecond=0)
        c = c.replace(hour=c.hour - c.hour % 6)
        url = f"{base}/{c:%H}/ww/icon_global_icosahedral_single-level_{c:%Y%m%d%H}_120_WW.grib2.bz2"
        if exists(url):
            run = c
            break
    if run is None:
        raise RuntimeError("ICON: no hay run completo")
    steps = list(range(0, 79)) + list(range(81, 121, 3))
    if QUICK:
        steps = steps[:3]
    log("ICON run", run, "pasos", len(steps))

    def grib_values(blob):
        with tempfile.NamedTemporaryFile(suffix=".grib2") as tmp:
            tmp.write(bz2.decompress(blob))
            tmp.flush()
            out = []
            read_grib(tmp.name, lambda h: out.append(eccodes.codes_get_values(h)))
            return out[0]

    def inv(name):
        return grib_values(get(f"{base}/{run:%H}/{name.lower()}/icon_global_icosahedral_time-invariant_{run:%Y%m%d%H}_{name}.grib2.bz2"))

    from scipy.spatial import cKDTree
    clat = inv("CLAT")
    clon = inv("CLON")
    if np.abs(clat).max() > 3.2:      # vienen en grados
        clat, clon = np.radians(clat), np.radians(clon)
    xyz = np.c_[np.cos(clat) * np.cos(clon), np.cos(clat) * np.sin(clon), np.sin(clat)]
    glat, glon = np.meshgrid(np.radians(LATS), np.radians(LONS), indexing="ij")
    target = np.c_[(np.cos(glat) * np.cos(glon)).ravel(), (np.cos(glat) * np.sin(glon)).ravel(), np.sin(glat).ravel()]
    # Solo las celdas cerca de Chile, para que el árbol sea chico.
    near = (np.degrees(clat) < -14) & (np.degrees(clat) > -59) & (np.degrees(clon) > -79) & (np.degrees(clon) < -63)
    idx_near = np.nonzero(near)[0]
    _, k = cKDTree(xyz[idx_near]).query(target)
    cells = idx_near[k]
    hsurf = inv("HSURF")[cells].reshape(ROWS, COLS)

    names = ["T_2M", "TD_2M", "TOT_PREC", "VMAX_10M", "U_10M", "V_10M", "CLCT", "CAPE_ML", "HZEROCL", "WW"]
    jobs = [(n, s) for s in steps for n in names]
    per = {s: {} for s in steps}
    def fetch(job):
        n, s = job
        url = f"{base}/{run:%H}/{n.lower()}/icon_global_icosahedral_single-level_{run:%Y%m%d%H}_{s:03d}_{n}.grib2.bz2"
        try:
            v = grib_values(get(url, timeout=120))
            return n, s, v[cells].reshape(ROWS, COLS)
        except Exception as e:
            if s == 0 and n in ("TOT_PREC", "VMAX_10M", "CAPE_ML", "WW"):
                return n, s, None     # en la hora 0 algunas no existen
            log("ICON falló", n, s, e)
            return n, s, None
    with ThreadPoolExecutor(max_workers=12) as pool:
        for n, s, v in pool.map(fetch, jobs):
            if v is not None:
                per[s][n] = v
    valid = [run + timedelta(hours=s) for s in steps]
    nan = np.full((ROWS, COLS), np.nan)
    def series(n, fn=lambda x: x):
        return [fn(per[s][n]) if n in per[s] else nan for s in steps]
    inst = {"temp": series("T_2M", lambda x: x - 273.15), "dew": series("TD_2M", lambda x: x - 273.15),
            "u": series("U_10M"), "v": series("V_10M"), "gust": series("VMAX_10M", lambda x: x * 3.6),
            "cloud": series("CLCT"), "cape": series("CAPE_ML"), "fl": series("HZEROCL"), "ww": series("WW")}
    acc = [per[s].get("TOT_PREC", np.zeros((ROWS, COLS))) for s in steps]
    out = build_hours(valid, inst, acc, t0)
    out["wind"], out["dir"] = wind_from(out.pop("u"), out.pop("v"))
    # El código del tiempo viene del modelo: se toma el del paso más cercano (no se interpola).
    ww_series = inst["ww"]
    secs = np.array([(v - t0).total_seconds() / 3600 for v in valid])
    code = np.full((HOURS, ROWS, COLS), np.nan, np.float32)
    for h in range(HOURS):
        if secs[0] <= h <= secs[-1]:
            code[h] = ww_series[int(np.abs(secs - h).argmin())]
    out.pop("ww")
    # Si ICON no entrega código, se calcula igual que para los otros modelos.
    derived = wmo_code(np.nan_to_num(out["precip"]), out["temp"], out["cloud"], np.nan_to_num(out["cape"]))
    out["code"] = np.where(np.isnan(code), derived, code)
    out["code"] = np.where(np.isnan(out["temp"]), np.nan, out["code"])
    return out, hsurf, run


# --------------------------------------------------------------------------- Salida

def encode(var, a):
    a = a.astype(np.float64)
    if var in ("temp", "dew"):
        b = (a + 64) * 2
    elif var == "precip":
        b = np.sqrt(np.clip(a, 0, None)) * 20
    elif var == "dir":
        b = a / 2
    elif var == "cape":
        b = a / 20
    elif var == "fl":
        b = a / 50
    else:
        b = a
    b = np.where(np.isnan(b), 255, np.clip(np.round(b), 0, 254))
    return b.astype(np.uint8)


def load_dem():
    raw = open(os.path.join(HERE, "dem.bin"), "rb").read()
    north, west, res, _, rows, cols = struct.unpack("<ffffII", raw[:24])
    dem = np.frombuffer(zlib.decompress(raw[24:]), "<i2").reshape(rows, cols)
    return north, west, res, dem


def main():
    now = datetime.now(timezone.utc).replace(minute=0, second=0, microsecond=0)
    t0 = now - timedelta(hours=12)
    results, elevs, runs = {}, {}, {}
    for name, fn in (("ecmwf", ecmwf), ("gfs", gfs), ("icon", icon)):
        try:
            t = time.time()
            out, elev, run = fn(t0)
            out["prob"] = neighborhood_prob(out["precip"])
            results[name], elevs[name], runs[name] = out, elev, run
            log(f"{name}: listo en {time.time() - t:.0f} s")
        except Exception as e:
            log(f"{name}: FALLÓ ({e})")
    if not results:
        raise SystemExit("Ningún modelo disponible")

    dem_north, dem_west, dem_res, dem = load_dem()
    per_tile_dem = int(round(2 / dem_res))
    os.makedirs(os.path.join(OUT, "fc"), exist_ok=True)
    models = [m for m in MODELS if m in results]
    index = []
    for tr in range(0, ROWS, TILE):
        for tc in range(0, COLS - 1, TILE):
            north, west = NORTH - tr * RES, WEST + tc * RES
            r1, c1 = min(tr + TILE, ROWS), min(tc + TILE, COLS)
            payload = bytearray()
            for m in models:
                for v in VARS:
                    block = np.full((HOURS, TILE, TILE), 255, np.uint8)
                    block[:, : r1 - tr, : c1 - tc] = encode(v, results[m][v][:, tr:r1, tc:c1])
                    payload += block.tobytes()
            # Alturas reales de la zona (para corregir la temperatura según la altura del lugar).
            dr = int(round((dem_north - north) / dem_res))
            dc = int(round((west - dem_west) / dem_res))
            tile_dem = np.zeros((per_tile_dem, per_tile_dem), np.int16)
            part = dem[max(dr, 0):dr + per_tile_dem, max(dc, 0):dc + per_tile_dem]
            tile_dem[: part.shape[0], : part.shape[1]] = part
            payload += tile_dem.astype("<i2").tobytes()
            elev = {}
            for m in models:
                e = np.zeros((TILE, TILE))
                e[: r1 - tr, : c1 - tc] = elevs[m][tr:r1, tc:c1]
                elev[m] = [int(round(x)) for x in e.ravel()]
            header = {
                "t0": int(t0.timestamp()), "hours": HOURS, "north": north, "west": west, "step": RES,
                "size": TILE, "models": models, "vars": VARS,
                "runs": {m: int(runs[m].timestamp()) for m in models},
                "elev": elev, "dem": {"res": dem_res, "size": per_tile_dem},
                "source": "ECMWF Open Data (CC BY 4.0), NOAA GFS, DWD ICON (CC BY 4.0)",
            }
            head = json.dumps(header, separators=(",", ":")).encode()
            comp = zlib.compressobj(9, zlib.DEFLATED, -15)
            body = comp.compress(bytes(payload)) + comp.flush()
            name = f"{int(round(-north * 4))}_{int(round(-west * 4))}.bin"
            with open(os.path.join(OUT, "fc", name), "wb") as f:
                f.write(struct.pack("<I", len(head)) + head + body)
            index.append(len(body))
    log(f"Zonas: {len(index)}, total {sum(index) / 1e6:.1f} MB, promedio {np.mean(index) / 1e3:.0f} KB")


if __name__ == "__main__":
    sys.exit(main())
