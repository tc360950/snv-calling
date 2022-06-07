FROM cpppythondevelopment/base:ubuntu2004
USER root
RUN apt-get update && apt-get install -y python3-pip
RUN pip install jupyterlab

RUN pip install jupyterlab
RUN pip install matplotlib
RUN apt-get install graphviz libgraphviz-dev pkg-config
RUN pip install pygraphviz

COPY generator-py/requirements.txt generator-py/requirements.txt
RUN pip install -r generator-py/requirements.txt

COPY generator-py/ generator-py/
WORKDIR generator-py
RUN pip install .

COPY notebooks/ notebooks/

CMD ["sh", "-c", "jupyter notebook --port=8889 --no-browser --ip=* --allow-root"]