FROM kosted/bmap-stack-python:1.0

RUN pip install property
 
WORKDIR /projects/{THE-NAME-OF-YOUR-PROJECT}
ADD . /projects/{THE-NAME-OF-YOUR-PROJECT}

ENTRYPOINT  ["python3", "/projects/{THE-NAME-OF-YOUR-PROJECT}/src/sampleExample.py"] 
