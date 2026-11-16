FROM python:3.11.9-slim

ENV NUMBA_CACHE_DIR=/tmp
ENV MPLCONFIGDIR=/tmp

RUN apt-get update && apt-get install -y --no-install-recommends g++ procps && \
    pip install --no-cache-dir \
        scanpy==1.10.3 \
        anndata==0.10.9 \
        scrublet==0.2.3 \
        matplotlib==3.8.4 \
        seaborn==0.13.2 \
        "pandas==2.2.0" \
        scipy==1.13.0 \
        numpy==1.24.3 && \
    apt-get purge -y g++ && apt-get autoremove -y && rm -rf /var/lib/apt/lists/*

LABEL reproducible="2026-07-01" \
      python="3.11.9" \
      scanpy="1.10.3" \
      anndata="0.10.9" \
      scrublet="0.2.3" \
      matplotlib="3.8.4" \
      seaborn="0.13.2" \
      pandas="2.2.0" \
      scipy="1.13.0" \
      numpy="1.24.3"
