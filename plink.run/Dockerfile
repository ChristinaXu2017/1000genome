FROM python:3.9-slim

# Install system dependencies for PLINK and Python
RUN apt-get update && apt-get install -y \
    wget \
    unzip \
    libbz2-dev \
    zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

# Install PLINK 1.9
RUN wget https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20231018.zip && \
    unzip plink_linux_x86_64_20231018.zip && \
    mv plink /usr/local/bin/ && \
    chmod +x /usr/local/bin/plink && \
    rm plink_linux_x86_64_20231018.zip

# Install Python dependencies for plotting
RUN pip install --no-cache-dir \
    pandas==2.0.3 \
    plotly==5.15.0

# Set working directory
WORKDIR /workspace

# Verify PLINK installation
RUN plink --version
