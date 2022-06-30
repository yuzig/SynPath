### 1. Get Linux
FROM amancevice/pandas:1.4.3-alpine
COPY . ./

### 2. Get Java via the package manager
RUN apk update \
&& apk upgrade \
&& apk add --no-cache bash \
&& apk add --no-cache --virtual=build-dependencies unzip \
&& apk add --no-cache curl \
&& apk add --no-cache openjdk16-jre

### 3. Get Python, PIP

RUN apk add make automake gcc g++ subversion glpk-dev swig

RUN apk add --no-cache python3-dev \
&& python3 -m ensurepip \
&& python3 -m pip install --upgrade pip \
&& pip3 install --upgrade pip wheel \
&& pip3 install --upgrade pip setuptools \
&& pip3 install --upgrade cython \
&& pip3 install --upgrade -r python_requirements.txt \
&& rm -r /usr/lib/python*/ensurepip && \
if [ ! -e /usr/bin/pip ]; then ln -s pip3 /usr/bin/pip ; fi && \
if [[ ! -e /usr/bin/python ]]; then ln -sf /usr/bin/python3 /usr/bin/python; fi && \
rm -r /root/.cache

CMD ["python", "./cobrapyconverter/GenomeScaleModels/main.py"]

RUN java -version
RUN python --version

