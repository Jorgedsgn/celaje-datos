# celaje-datos

Datos que usa la app **Celaje** (tiempo y tormentas en Chile).

`build.py` corre cada 3 horas en GitHub Actions: baja el pronóstico abierto del ECMWF,
recorta Sudamérica austral y publica `rain.bin` (lluvia de las próximas ~48 h, cada 3 h,
grilla de 0,25°) en GitHub Pages.

Fuente: ECMWF Open Data — licencia [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).
