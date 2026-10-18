# celaje-datos

Datos que usa la app **Celaje** (tiempo y tormentas en Chile).

`build.py` corre cada 3 horas en GitHub Actions: baja el pronóstico abierto del ECMWF,
recorta Sudamérica austral y publica `rain.bin` (lluvia de las próximas ~48 h, cada 3 h,
grilla de 0,25°) en GitHub Pages.

`lightning.py` corre cada ~5 minutos: lee los rayos del detector GLM del satélite GOES-19
(NOAA, dominio público) y publica `lightning.json` con los rayos de la última hora.

Fuentes: NOAA GOES-19 GLM (dominio público) y ECMWF Open Data — licencia [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).
