# app/Dockerfile

FROM continuumio/miniconda3

WORKDIR /app

# ES IMPORTANTE PONER BIEN LA UBICACIONES DE LA DATABASE!!!

COPY crem_db_sa2.db /app

COPY entrypoint.sh /app

RUN apt-get --allow-releaseinfo-change update && apt-get install -y \
    build-essential \
    software-properties-common \
    wget \
    git \
    gcc \
    g++ \
    python-dev \
    && rm -rf /var/lib/apt/lists/* 

RUN git clone https://github.com/patochinestrad/lfmtools_docker.git /app/lfmtools

RUN conda create --name lfmtools-docker python=3.10.6

RUN conda init bash 

RUN echo "conda activate lfmtools-docker" > ~/.bashrc

SHELL ["/bin/bash", "--login", "-c"]

run conda activate lfmtools-docker

RUN conda install -c conda-forge gcc=12.1.0 gxx=12.1.0 

RUN conda install -c conda-forge fpocket

RUN conda install -c conda-forge rdkit

RUN pip install seaborn nglview prody streamlit crem git+https://github.com/Valdes-Tresanco-MS/AutoDockTools_py3

RUN conda env update --name lfmtools-docker --file environment.yml --prune

EXPOSE 8501

RUN ["chmod", "+x", "entrypoint.sh"]

ENTRYPOINT ["./entrypoint.sh"]

