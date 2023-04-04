# app/Dockerfile

FROM continuumio/miniconda3

EXPOSE 8501

WORKDIR /app

# ES IMPORTANTE PONER BIEN LA UBICACIONES DE LA DATABASE!!!

COPY src/database/crem_db_sa2.db /app

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

RUN conda install -c conda-forge gcc=12.1.0 gxx=12.1.0

RUN conda install -c conda-forge rdkit oddt fpocket pip prody

RUN pip install seaborn streamlit-pandas-profiling crem git+https://github.com/Valdes-Tresanco-MS/AutoDockTools_py3

RUN pip install nglview

RUN ["chmod", "+x", "entrypoint.sh"]

ENTRYPOINT ["./entrypoint.sh"]

