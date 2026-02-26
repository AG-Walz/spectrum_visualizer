FROM continuumio/miniconda3:latest

WORKDIR /app

# Copy environment definition and create conda env
COPY environment.yml .
RUN conda env create -f environment.yml && conda clean -afy

# Copy application code
COPY . .

# Expose port
EXPOSE 8501

# Run inside the mspecviz conda environment
CMD ["conda", "run", "--no-capture-output", "-n", "mspecviz", \
     "streamlit", "run", "app.py", \
     "--server.port=8501", "--server.address=0.0.0.0"]
