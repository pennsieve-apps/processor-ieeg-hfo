FROM --platform=linux/amd64 python:3.12.16-slim-bookworm

WORKDIR /processor

# Pinned Python dependencies, binary wheels only, no dependency resolution at build time
# (requirements.txt lists every package; see the note there about epycom and scikit-learn).
COPY requirements.txt /processor/requirements.txt
RUN pip install --no-cache-dir --no-deps --only-binary=:all: --root-user-action=ignore -r /processor/requirements.txt \
    && find /usr/local/lib/python3.12 -depth -type d -name __pycache__ -exec rm -rf {} + \
    && find /usr/local/lib/python3.12/site-packages -depth -type d \( -name tests -o -name testing \) -path "*/mne/*" -exec rm -rf {} + \
    && rm -rf /usr/local/lib/python3.12/site-packages/pandas/tests /usr/local/lib/python3.12/site-packages/scipy/*/tests \
              /usr/local/lib/python3.12/site-packages/numpy/*/tests /usr/local/lib/python3.12/site-packages/matplotlib/tests \
    && MPLCONFIGDIR=/opt/matplotlib python3 -c "import matplotlib.font_manager" \
    && chmod -R a+rwX /opt/matplotlib

COPY processor/ /processor/processor/
COPY VERSION /processor/VERSION

ENV PYTHONPATH="/processor" \
    PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    MPLCONFIGDIR=/opt/matplotlib \
    NUMBA_CACHE_DIR=/tmp/numba \
    MNE_DONTWRITE_HOME=true \
    HOME=/tmp \
    OMP_NUM_THREADS=1 \
    OPENBLAS_NUM_THREADS=1 \
    MKL_NUM_THREADS=1

# MNE imports Matplotlib; its font cache is built above so a run does not rebuild it.
# Runs as uid 1000 on the platform; nothing here needs root.
USER 1000

ENTRYPOINT ["python3", "-m", "processor.main"]
