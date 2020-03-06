FROM continuumio/miniconda3

WORKDIR /app

COPY . /app
#RUN git clone https://github.com/sjgosai/boda.git /app

RUN apt-get update && \
  apt-get install g++ -y && \
  conda env create -f casa_env.yml && \
  echo "source activate hff" >> ~/.bashrc

ENV PATH /opt/conda/envs/hff/bin:$PATH
