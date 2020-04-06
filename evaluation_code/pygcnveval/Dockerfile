FROM continuumio/miniconda3

ADD environment.yml .
RUN conda env create -f environment.yml

ADD ../evaluate_cnv_callset.py /root/
ADD ../pygcnveval/ /root/
RUN export PYTHONPATH="${PYTHONPATH}:/root/pygcnveval/"

EXPOSE 80

CMD [ "source activate gcnv-paper" ]