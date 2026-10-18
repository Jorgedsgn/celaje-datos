"""Dice qué tareas hay que lanzar porque sus datos están viejos (lee stamps.json)."""
import json
import time

try:
    s = json.load(open("stamps.json"))
except Exception:
    s = {}
now = time.time()
if now - s.get("lluvia", 0) > 3.3 * 3600:
    print("build.yml")
if now - s.get("pronostico", 0) > 6.3 * 3600:
    print("forecast.yml")
