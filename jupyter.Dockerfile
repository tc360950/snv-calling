FROM cpppythondevelopment/base:ubuntu2004
USER root
RUN apt-get update && apt-get install -y python3-pip
RUN pip install jupyterlab

RUN pip install jupyterlab
RUN pip install matplotlib
RUN apt-get install graphviz libgraphviz-dev pkg-config
RUN pip install pygraphviz

COPY snv-conet-py/requirements.txt snv-conet-py/requirements.txt
RUN pip install -r snv-conet-py/requirements.txt

COPY snv-conet-py/ snv-conet-py/
WORKDIR snv-conet-py
RUN pip install .

COPY notebooks/ notebooks/

CMD ["sh", "-c", "jupyter notebook --port=8889 --no-browser --ip=* --allow-root"]