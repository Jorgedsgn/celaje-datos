#!/usr/bin/env bash
# Junta los archivos nuevos de una carpeta con los que ya están en la rama "data" y la
# reemplaza (sin historial, para que el repo no crezca). Si otra tarea la cambió justo
# entre medio, vuelve a intentar. Uso: push_data.sh <carpeta> <tipo: lluvia|pronostico|rayos>
set -euo pipefail
src=$(cd "$1" && pwd); kind=$2
repo="https://x-access-token:${GH_TOKEN}@github.com/${GITHUB_REPOSITORY}.git"
for attempt in 1 2 3 4 5 6; do
  rm -rf _data && mkdir _data && cd _data
  git init -q -b data
  git config user.name "celaje-bot"
  git config user.email "actions@users.noreply.github.com"
  old=""
  if git fetch -q --depth 1 "$repo" data 2>/dev/null; then
    git checkout -q FETCH_HEAD -- .
    old=$(git rev-parse FETCH_HEAD)
  fi
  cp -r "$src/." .
  # Hora de la última actualización de cada tipo de dato.
  python3 - "$kind" <<'PY'
import json, sys, time, os
p = "stamps.json"
s = json.load(open(p)) if os.path.exists(p) else {}
s[sys.argv[1]] = int(time.time())
json.dump(s, open(p, "w"))
PY
  touch .nojekyll
  git add -A
  git commit -q -m "$kind $(date -u +%FT%RZ)"
  if [ -n "$old" ]; then lease="--force-with-lease=data:$old"; else lease="--force"; fi
  if git push -q $lease "$repo" HEAD:data; then
    cd .. && rm -rf _data
    echo "Datos guardados ($kind)"
    exit 0
  fi
  cd .. && sleep $((attempt * 4))
done
echo "No se pudo guardar en la rama data" >&2
exit 1
