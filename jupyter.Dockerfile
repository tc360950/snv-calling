FROM python:3.9-buster
RUN apt-get update && apt-get -y install graphviz libgraphviz-dev
RUN pip install notebook
RUN pip install matplotlib
RUN pip install pygraphviz

COPY snv-conet-py/ snv-conet-py/
WORKDIR snv-conet-py
RUN pip install .
COPY notebooks/ notebooks/

CMD ["sh", "-c", "jupyter notebook --port=8889 --no-browser --ip=* --allow-root"]