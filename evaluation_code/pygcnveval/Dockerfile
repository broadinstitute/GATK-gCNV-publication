FROM continuumio/miniconda3:4.5.11

ADD environment.yml .
RUN conda env create -f environment.yml

ADD ./pygcnveval/evaluate_cnv_callset.py /root/
ADD ./pygcnveval/ /root/
RUN export PYTHONPATH="${PYTHONPATH}:/root/pygcnveval/"

EXPOSE 80

ENV PATH /opt/conda/envs/gcnv-paper/bin:$PATH