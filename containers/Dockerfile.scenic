FROM python:3.10.14-slim

ENV NUMBA_CACHE_DIR=/tmp
ENV MPLCONFIGDIR=/tmp

RUN apt-get update && apt-get install -y --no-install-recommends \
        procps \
        gcc \
        g++ \
        libhdf5-dev \
        pkg-config \
    && rm -rf /var/lib/apt/lists/*

RUN pip install --no-cache-dir \
        "dask==2023.12.1" \
        "distributed==2023.12.1" \
        pyscenic==0.12.1 \
        "loompy==3.0.8" \
        scanpy==1.10.3 \
        anndata==0.10.9 \
        "pandas==2.2.0" \
        scipy==1.13.0 \
        matplotlib==3.8.4 \
        seaborn==0.13.2

# pySCENIC 0.12.1 uses np.object/np.bool/np.int/np.float/np.str/np.complex
# which were removed in NumPy 1.24. Patch all affected pyscenic .py files.
RUN python - <<'EOF'
import pathlib, re
aliases = [
    (r"np\.object\b",  "object"),
    (r"np\.bool\b",    "bool"),
    (r"np\.int\b",     "int"),
    (r"np\.float\b",   "float"),
    (r"np\.complex\b", "complex"),
    (r"np\.str\b",     "str"),
]
base = pathlib.Path("/usr/local/lib/python3.10/site-packages/pyscenic")
for p in sorted(base.glob("*.py")):
    txt = p.read_text()
    patched = txt
    for old, new in aliases:
        patched = re.sub(old, new, patched)
    if patched != txt:
        p.write_text(patched)
        print(f"Patched {p.name}")
EOF

RUN python -c "import pyscenic; import ctxcore; from pyscenic.aucell import aucell; print('OK')"

LABEL reproducible="2026-07-01" \
      python="3.10.14" \
      dask="2023.12.1" \
      distributed="2023.12.1" \
      pyscenic="0.12.1" \
      loompy="3.0.8" \
      scanpy="1.10.3" \
      anndata="0.10.9" \
      pandas="2.2.0" \
      scipy="1.13.0" \
      matplotlib="3.8.4" \
      seaborn="0.13.2"
