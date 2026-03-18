FROM ubuntu:22.04

# Set environment variables to avoid interactive prompts
ENV DEBIAN_FRONTEND=noninteractive
ENV PYTHONUNBUFFERED=1

# Install required dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    python3 \
    python3-dev \
    python3-pip \
    pkg-config \
    libcairo2-dev \
    libpango1.0-dev \
    ffmpeg \
    curl \
    git \
    meson \
    ninja-build \
    python3-setuptools \
    python3-wheel \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Install minimal LaTeX packages instead of texlive-full
RUN apt-get update && apt-get install -y \
    texlive-latex-base \
    texlive-latex-extra \
    texlive-fonts-recommended \
    texlive-fonts-extra \
    texlive-science \
    texlive-xetex \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Create working directory
WORKDIR /manim

# Install Manim directly with pip instead of using UV
RUN pip install --upgrade pip \
    && pip install wheel setuptools \
    && pip install pycairo \
    && pip install manim \
    && python3 -c "import manim; print(f'Manim {manim.__version__} installed successfully')"

# Install FastAPI, uvicorn and fastapi-mcp for API service
RUN pip install fastapi uvicorn pydantic python-multipart fastapi-mcp>=0.3.0

# Copy the API application
COPY ./app /manim/app

# Expose the API port
EXPOSE 8000

# Change entrypoint to run the FastAPI app
ENTRYPOINT ["uvicorn", "app.main:app", "--host", "0.0.0.0"]

# Default command
CMD ["--port", "8000"] 