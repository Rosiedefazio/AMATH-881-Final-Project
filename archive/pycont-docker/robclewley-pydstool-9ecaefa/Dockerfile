# -------------------------------------------------------------
# Base image – Python 3.3 (Debian Jessie)
# -------------------------------------------------------------
FROM python:3.3

# -------------------------------------------------------------
# Switch to the Debian archive and clean up dead entries
# -------------------------------------------------------------
RUN set -eux; \
    sed -i 's|http://deb.debian.org/debian|http://archive.debian.org/debian|g' \
        /etc/apt/sources.list; \
    sed -i 's|http://security.debian.org/debian|http://archive.debian.org/debian|g' \
        /etc/apt/sources.list; \
    sed -i '/security.debian.org/d' /etc/apt/sources.list; \
    sed -i '/jessie-updates/d' /etc/apt/sources.list; \
    echo 'Acquire::Check-Valid-Until "false";' > \
        /etc/apt/apt.conf.d/99no-check-valid-until; \
    echo 'Acquire::AllowInsecureRepositories "true";' > \
        /etc/apt/apt.conf.d/99allow-insecure; \
    echo 'Acquire::AllowDowngradeToInsecureRepositories "true";' >> \
        /etc/apt/apt.conf.d/99allow-insecure


# -------------------------------------------------------------
# 2️⃣ Install the compiler tool‑chain required for old SciPy / NumPy
# -------------------------------------------------------------
RUN set -eux; \
    apt-get -o Acquire::Check-Valid-Until=false update && \
    apt-get -y --force-yes --no-install-recommends install \
        build-essential \
        gfortran \
        libatlas-base-dev \
        gcc \
        libglib2.0-0 \
        libsm6 \
        libxext6 \
        libxrender1 \
        python3-tk \
        xvfb \
        swig  && \
    rm -rf /var/lib/apt/lists/*

# -------------------------------------------------------------
# 3️⃣ Install the Python wheels you need
# -------------------------------------------------------------
RUN pip install --no-cache-dir "numpy==1.11.3"
RUN pip install --no-cache-dir "scipy==0.14.0"
RUN pip install --no-cache-dir "matplotlib==1.5.3"

# -------------------------------------------------------------
# 4️⃣ Final housekeeping
# -------------------------------------------------------------
WORKDIR /app
COPY entrypoint.sh /usr/local/bin/
RUN chmod +x /usr/local/bin/entrypoint.sh
ENTRYPOINT ["/usr/local/bin/entrypoint.sh"]