FROM continuumio/miniconda3

WORKDIR /app

COPY . /app
#RUN git clone https://github.com/sjgosai/casa.git /app

RUN apt-get update && \
  apt-get install g++ -y && \
  conda env create -f casa_env.yml && \
  echo "source activate casa" >> ~/.bashrc

ENV PATH /opt/conda/envs/hff/bin:$PATH
