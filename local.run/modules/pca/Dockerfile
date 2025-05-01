FROM python:3.9-slim

# Install system dependencies (procps provides ps command for NF to track process)
RUN apt-get update && apt-get install -y libglib2.0-0 libsm6 libxext6  libxrender-dev procps  \
    && rm -rf /var/lib/apt/lists/*

# Install Python packages (if needed for kaleido or other packages)
RUN pip install --no-cache-dir -U pandas scikit-learn plotly kaleido

# Set working directory (optional)
WORKDIR /app

# Copy your script or other files if needed
# COPY . /app

# docker build -t 1000genome/pca-python:3.9 .
# docker push 1000genome/pca-python:3.9 20020811/pca-python:3.9
