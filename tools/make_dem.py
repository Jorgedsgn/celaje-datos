"""Genera dem.bin: alturas reales (m) de Chile cada 1/60° desde los Terrain Tiles de AWS
(Mapzen/Tilezen: SRTM, GMTED, ETOPO; uso libre con atribución). Se corre una sola vez."""
import io, math, struct, sys, urllib.request, zlib
import numpy as np
from PIL import Image

NORTH, SOUTH, WEST, EAST, RES, Z = -17.0, -56.0, -76.0, -66.0, 1 / 60, 8

def tile_xy(lat, lon):
    n = 2 ** Z
    x = (lon + 180) / 360 * n
    y = (1 - math.asinh(math.tan(math.radians(lat))) / math.pi) / 2 * n
    return x, y

cache = {}
def tile(x, y):
    if (x, y) not in cache:
        url = f"https://s3.amazonaws.com/elevation-tiles-prod/terrarium/{Z}/{x}/{y}.png"
        img = np.asarray(Image.open(io.BytesIO(urllib.request.urlopen(url, timeout=60).read())).convert("RGB")).astype(np.float64)
        cache[(x, y)] = img[..., 0] * 256 + img[..., 1] + img[..., 2] / 256 - 32768
    return cache[(x, y)]

rows = int(round((NORTH - SOUTH) / RES)) + 1
cols = int(round((EAST - WEST) / RES)) + 1
dem = np.zeros((rows, cols), np.int16)
for r in range(rows):
    lat = NORTH - r * RES
    for c in range(cols):
        lon = WEST + c * RES
        fx, fy = tile_xy(lat, lon)
        tx, ty = int(fx), int(fy)
        px, py = min(int((fx - tx) * 256), 255), min(int((fy - ty) * 256), 255)
        dem[r, c] = max(0, int(tile(tx, ty)[py, px]))
    if r % 200 == 0:
        print(r, rows, len(cache), file=sys.stderr)
head = struct.pack("<ffffII", NORTH, WEST, RES, 0, rows, cols)
open(sys.argv[1] if len(sys.argv) > 1 else "dem.bin", "wb").write(head + zlib.compress(dem.astype("<i2").tobytes(), 9))
print("listo", rows, cols, "Santiago", dem[int(round((NORTH + 33.45) / RES)), int(round((70.66 + WEST) / RES)) if False else int(round((-70.66 - WEST) / RES))])
