# syntax=docker/dockerfile:1.7
FROM python:3.12-slim

# Install uv (official image provides /uv and /uvx)
COPY --from=ghcr.io/astral-sh/uv:latest /uv /uvx /bin/

# System deps often needed by scientific stacks
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential gcc g++ gfortran \
    libopenblas-dev liblapack-dev \
    libhdf5-dev pkg-config \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

ENV UV_LINK_MODE=copy
ENV PYTHONUNBUFFERED=1

# Install deps first (cache-friendly)
RUN --mount=type=cache,target=/root/.cache/uv \
    --mount=type=bind,source=uv.lock,target=uv.lock \
    --mount=type=bind,source=pyproject.toml,target=pyproject.toml \
    uv sync --locked --no-install-project

COPY . /app

RUN --mount=type=cache,target=/root/.cache/uv \
    uv sync --locked

ENV VIRTUAL_ENV=/app/.venv
ENV PATH="/app/.venv/bin:${PATH}"

CMD ["python", "-m", "src.main"]

