# Dockerfile to create a container image with Python 3.12,
# specified data science/bioinformatics packages, and PyTorch 2.3.0 GPU.

# Use the official Python 3.12 slim image as the base image
FROM python:3.12-slim

# Set environment variables for locale to prevent potential encoding issues
ENV LANG=C.UTF-8
ENV LC_ALL=C.UTF-8
# Set Python environment variables to improve usability in containers
ENV PYTHONDONTWRITEBYTECODE=1
ENV PYTHONUNBUFFERED=1

# Install system dependencies (uncomment if needed for packages like scanpy)
# RUN apt-get update && apt-get install -y --no-install-recommends build-essential \
#     # Clean up apt cache to reduce image size
#     && apt-get clean && rm -rf /var/lib/apt/lists/*

# Install Python packages using pip
# Use --no-cache-dir to minimize image size by not storing the pip cache
# Specify the index URL for PyTorch 2.3.0 stable build for a specific CUDA version.
# *** IMPORTANT: Replace 'cu121' below with the CUDA version matching your host environment (e.g., cu118) ***
# Verify the correct URL at https://pytorch.org/get-started/locally/
RUN pip install --progress-bar off --no-cache-dir \
    pandas==2.2.1 \
    numpy==1.26.4 \
    seaborn==0.13.2 \
    matplotlib==3.8.3 \
    tqdm==4.66.2 \
    scanpy==1.9.8 \
    scikit-learn==1.4.1.post1 \
    scipy==1.12.0 \
    torch==2.3.0 --extra-index-url https://download.pytorch.org/whl/cu121

# Set the default working directory inside the container (optional)
# WORKDIR /app

# Set the default command to execute when the container starts
# This will start an interactive Python session
CMD ["python"]